"""
Resource Management Module for Scalable KBase Protein Analysis

This module provides comprehensive resource management, monitoring, and optimization
for scalable protein analysis workflows in KBase environments.
"""

import os
import psutil
import gc
import logging
import time
import threading
from typing import Dict, Any, Optional, Callable, List
from dataclasses import dataclass, field
from contextlib import contextmanager
import numpy as np

logger = logging.getLogger(__name__)

@dataclass
class ResourceLimits:
    """
    Server-aware resource limits for KBase DOE servers.
    
    These limits are designed to be respectful of shared server resources
    and ensure efficient operation without impacting other processes.
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/core/resource_manager.py:20
    USED BY: ResourceManager for enforcing resource constraints
    """
    # Percentage-based limits for shared server environments
    max_memory_percent: float = 60.0  # Use max 60% of available memory
    max_cpu_percent: float = 70.0     # Use max 70% of CPU to leave room for other processes
    max_disk_percent: float = 80.0    # Use max 80% of available disk space
    
    # Conservative defaults for DOE server environments
    batch_size_proteins: int = 500     # Smaller batches for shared resources
    max_concurrent_tasks: int = 2      # Limit concurrent tasks on shared servers
    gc_threshold_percent: float = 50.0 # Trigger GC at 50% memory usage
    
    # Dynamic scaling parameters for server efficiency
    scale_factor: float = 0.6          # Conservative scaling for shared environments
    min_batch_size: int = 50           # Smaller minimum for resource efficiency
    max_batch_size: int = 2000         # Reasonable maximum for server environments
    
    # Server detection and adaptation
    auto_detect_server_env: bool = True
    server_safety_margin: float = 0.2  # 20% safety margin for server operations

@dataclass
class ResourceMetrics:
    """Current resource usage metrics."""
    memory_used_gb: float = 0.0
    memory_percent: float = 0.0
    cpu_percent: float = 0.0
    disk_used_gb: float = 0.0
    active_threads: int = 0
    gc_collections: int = 0
    timestamp: float = field(default_factory=time.time)

class ResourceManager:
    """
    Server-aware resource manager for KBase DOE server environments.
    
    This manager is designed for shared server environments and uses percentage-based
    limits to ensure respectful resource usage that doesn't impact other processes.
    
    Features:
    - Server environment detection and adaptation
    - Percentage-based resource limits for shared environments
    - Dynamic memory management with conservative garbage collection
    - CPU utilization monitoring with server-safe throttling
    - Disk space monitoring with safety margins
    - Batch size optimization for server efficiency
    - Thread pool management with server-aware limits
    - Performance profiling optimized for shared resources
    
    CLASS LOCATION: lib/kbase_protein_query_module/src/core/resource_manager.py:65
    USED BY: ParallelProcessor, WorkflowOrchestrator, all pipeline stages
    EXTENDS: None (standalone class)
    """
    
    def __init__(self, limits: Optional[ResourceLimits] = None, max_memory_gb: Optional[float] = None, max_cpu_cores: Optional[int] = None):
        self.limits = limits or ResourceLimits()
        
        # Override limits if specific values provided
        if max_memory_gb is not None:
            self.limits.max_memory_percent = (max_memory_gb / (psutil.virtual_memory().total / (1024**3))) * 100
        if max_cpu_cores is not None:
            self.limits.max_cpu_percent = (max_cpu_cores / psutil.cpu_count()) * 100
        
        self.metrics_history: List[ResourceMetrics] = []
        self.active_tasks: Dict[str, Any] = {}
        self._monitoring = False
        self._monitor_thread: Optional[threading.Thread] = None
        
        # Performance tracking
        self.operation_times: Dict[str, List[float]] = {}
        self.memory_peaks: Dict[str, float] = {}
        
        # Additional attributes for test compatibility
        self.max_memory_gb = max_memory_gb or (psutil.virtual_memory().total / (1024**3)) * 0.6
        self.max_cpu_cores = max_cpu_cores or psutil.cpu_count()
        self.registered_processes: Dict[str, Dict[str, Any]] = {}
        self.allocated_memory: float = 0.0
        self.allocated_cpu: int = 0
        self.active_processes: Dict[str, Dict[str, Any]] = {}
        
        logger.info(f"ResourceManager initialized with limits: {self.limits}")
    
    def get_current_metrics(self) -> ResourceMetrics:
        """Get current system resource metrics."""
        memory = psutil.virtual_memory()
        cpu_percent = psutil.cpu_percent(interval=0.1)
        
        # Disk usage for current working directory
        disk = psutil.disk_usage('.')
        disk_used_gb = (disk.total - disk.free) / (1024**3)
        
        metrics = ResourceMetrics(
            memory_used_gb=memory.used / (1024**3),
            memory_percent=memory.percent,
            cpu_percent=cpu_percent,
            disk_used_gb=disk_used_gb,
            active_threads=threading.active_count(),
            gc_collections=sum(gc.get_stats()[i]['collections'] for i in range(len(gc.get_stats())))
        )
        
        self.metrics_history.append(metrics)
        
        # Keep only recent metrics (last 1000 entries)
        if len(self.metrics_history) > 1000:
            self.metrics_history = self.metrics_history[-1000:]
        
        return metrics
    
    def check_resource_availability(self) -> bool:
        """
        Check if resources are available for processing using server-aware percentage limits.
        
        This method uses percentage-based thresholds designed for shared KBase DOE servers
        to ensure respectful resource usage that doesn't impact other running processes.
        """
        metrics = self.get_current_metrics()
        
        # Check memory usage against percentage limit
        if metrics.memory_percent > self.limits.max_memory_percent:
            logger.warning(f"Memory usage exceeds server limit: {metrics.memory_percent:.1f}% > {self.limits.max_memory_percent}%")
            return False
        
        # Check CPU usage against percentage limit
        if metrics.cpu_percent > self.limits.max_cpu_percent:
            logger.warning(f"CPU usage exceeds server limit: {metrics.cpu_percent:.1f}% > {self.limits.max_cpu_percent}%")
            return False
        
        # Check disk usage against percentage limit
        disk_total_gb = psutil.disk_usage('.').total / (1024**3)
        disk_percent = (metrics.disk_used_gb / disk_total_gb) * 100
        if disk_percent > self.limits.max_disk_percent:
            logger.warning(f"Disk usage exceeds server limit: {disk_percent:.1f}% > {self.limits.max_disk_percent}%")
            return False
        
        return True
    
    def optimize_batch_size(self, base_batch_size: int, data_size: int) -> int:
        """
        Dynamically optimize batch size for server environments using percentage-based limits.
        
        This method ensures batch sizes are appropriate for shared DOE server resources
        by using percentage-based calculations rather than absolute values.
        """
        metrics = self.get_current_metrics()
        
        # Calculate resource availability as percentages
        memory_available_percent = 100 - metrics.memory_percent
        cpu_available_percent = 100 - metrics.cpu_percent
        
        # Calculate scaling factors based on server-safe percentages
        memory_factor = min(memory_available_percent / (100 - self.limits.max_memory_percent), 1.0)
        cpu_factor = min(cpu_available_percent / (100 - self.limits.max_cpu_percent), 1.0)
        
        # Apply server safety margin
        safety_factor = 1.0 - self.limits.server_safety_margin
        
        # Calculate final scaling factor
        scaling_factor = min(memory_factor, cpu_factor) * self.limits.scale_factor * safety_factor
        optimized_batch_size = int(base_batch_size * max(scaling_factor, 0.1))  # Minimum 10% scaling
        
        # Apply server-appropriate bounds
        optimized_batch_size = max(self.limits.min_batch_size, optimized_batch_size)
        optimized_batch_size = min(self.limits.max_batch_size, optimized_batch_size)
        optimized_batch_size = min(optimized_batch_size, data_size)
        
        logger.debug(f"Server-optimized batch size: {optimized_batch_size} (from {base_batch_size}, "
                    f"memory_factor={memory_factor:.2f}, cpu_factor={cpu_factor:.2f})")
        return optimized_batch_size
    
    @contextmanager
    def resource_context(self, operation_name: str):
        """Context manager for resource monitoring during operations."""
        start_time = time.time()
        start_metrics = self.get_current_metrics()
        
        try:
            self.active_tasks[operation_name] = start_time
            logger.debug(f"Starting operation: {operation_name}")
            yield self
            
        finally:
            # Record operation time
            operation_time = time.time() - start_time
            if operation_name not in self.operation_times:
                self.operation_times[operation_name] = []
            self.operation_times[operation_name].append(operation_time)
            
            # Record memory peak
            end_metrics = self.get_current_metrics()
            memory_used = end_metrics.memory_used_gb - start_metrics.memory_used_gb
            if memory_used > 0:
                self.memory_peaks[operation_name] = max(
                    self.memory_peaks.get(operation_name, 0), memory_used
                )
            
            # Cleanup
            if operation_name in self.active_tasks:
                del self.active_tasks[operation_name]
            
            # Trigger garbage collection if memory usage is high
            if end_metrics.memory_used_gb > self.limits.gc_threshold_mb / 1024:
                gc.collect()
            
            logger.debug(f"Completed operation: {operation_name} in {operation_time:.2f}s")
    
    def get_performance_summary(self) -> Dict[str, Any]:
        """Get performance summary for all operations."""
        summary = {}
        
        for operation, times in self.operation_times.items():
            summary[operation] = {
                'count': len(times),
                'avg_time': np.mean(times),
                'max_time': np.max(times),
                'min_time': np.min(times),
                'total_time': np.sum(times),
                'memory_peak_gb': self.memory_peaks.get(operation, 0)
            }
        
        # Overall system metrics
        if self.metrics_history:
            recent_metrics = self.metrics_history[-10:]  # Last 10 measurements
            summary['system'] = {
                'avg_memory_percent': np.mean([m.memory_percent for m in recent_metrics]),
                'avg_cpu_percent': np.mean([m.cpu_percent for m in recent_metrics]),
                'max_memory_gb': np.max([m.memory_used_gb for m in recent_metrics]),
                'active_tasks': len(self.active_tasks)
            }
        
        return summary
    
    def check_memory_availability(self, required_memory_gb: float = None) -> bool:
        """Check if memory is available."""
        if required_memory_gb is not None:
            available_memory = self.max_memory_gb - self.allocated_memory
            return available_memory >= required_memory_gb
        metrics = self.get_current_metrics()
        return metrics.memory_percent < self.limits.max_memory_percent
    
    def check_cpu_availability(self, required_cpu_cores: int = None) -> bool:
        """Check if CPU is available."""
        if required_cpu_cores is not None:
            available_cpu = self.max_cpu_cores - self.allocated_cpu
            return available_cpu >= required_cpu_cores
        metrics = self.get_current_metrics()
        return metrics.cpu_percent < self.limits.max_cpu_percent
    
    def allocate_memory(self, size_gb: float) -> bool:
        """Allocate memory if available."""
        if self.check_memory_availability():
            logger.debug(f"Allocated {size_gb}GB memory")
            return True
        return False
    
    def deallocate_memory(self, size_gb: float):
        """Deallocate memory."""
        logger.debug(f"Deallocated {size_gb}GB memory")
    
    def allocate_cpu(self, cores: int) -> bool:
        """Allocate CPU cores if available."""
        if self.check_cpu_availability():
            logger.debug(f"Allocated {cores} CPU cores")
            return True
        return False
    
    def deallocate_cpu(self, cores: int):
        """Deallocate CPU cores."""
        logger.debug(f"Deallocated {cores} CPU cores")
    
    def register_process(self, process_id: str, memory_gb: float, cpu_cores: int) -> bool:
        """Register a process with resource requirements."""
        if self.allocate_memory(memory_gb) and self.allocate_cpu(cpu_cores):
            self.registered_processes[process_id] = {
                'memory_gb': memory_gb,
                'cpu_cores': cpu_cores,
                'start_time': time.time()
            }
            return True
        return False
    
    def unregister_process(self, process_id: str) -> bool:
        """Unregister a process."""
        if process_id in self.registered_processes:
            process_info = self.registered_processes[process_id]
            self.deallocate_memory(process_info['memory_gb'])
            self.deallocate_cpu(process_info['cpu_cores'])
            del self.registered_processes[process_id]
            return True
        return False
    
    def get_resource_usage(self) -> Dict[str, Any]:
        """Get current resource usage."""
        metrics = self.get_current_metrics()
        return {
            'memory_gb': self.allocated_memory,
            'memory_used_gb': metrics.memory_used_gb,
            'memory_percent': metrics.memory_percent,
            'memory_percentage': (self.allocated_memory / self.max_memory_gb) * 100,
            'cpu_percent': metrics.cpu_percent,
            'cpu_percentage': (self.allocated_cpu / self.max_cpu_cores) * 100,
            'cpu_cores': self.allocated_cpu,
            'active_processes': len(self.registered_processes)
        }
    
    def get_available_resources(self) -> Dict[str, Any]:
        """Get available resources."""
        metrics = self.get_current_metrics()
        available_memory = self.max_memory_gb - self.allocated_memory
        available_cpu = self.max_cpu_cores - self.allocated_cpu
        return {
            'memory_gb': available_memory,
            'available_memory_gb': available_memory,
            'cpu_cores': available_cpu,
            'available_cpu_percent': 100 - metrics.cpu_percent,
            'available_memory_percent': 100 - metrics.memory_percent
        }
    
    def list_active_processes(self) -> List[str]:
        """List active process IDs."""
        return list(self.registered_processes.keys())
    
    def clear_all_processes(self):
        """Clear all registered processes."""
        for process_id in list(self.registered_processes.keys()):
            self.unregister_process(process_id)
    
    def allocate_memory(self, memory_gb: float) -> bool:
        """Allocate memory for a task."""
        if self.check_memory_availability(memory_gb):
            self.allocated_memory += memory_gb
            return True
        return False
    
    def deallocate_memory(self, memory_gb: float) -> None:
        """Deallocate memory from a task."""
        self.allocated_memory = max(0, self.allocated_memory - memory_gb)
    
    def allocate_cpu(self, cpu_cores: int) -> bool:
        """Allocate CPU cores for a task."""
        if self.check_cpu_availability(cpu_cores):
            self.allocated_cpu += cpu_cores
            return True
        return False
    
    def deallocate_cpu(self, cpu_cores: int) -> None:
        """Deallocate CPU cores from a task."""
        self.allocated_cpu = max(0, self.allocated_cpu - cpu_cores)
    
    def register_process(self, process_id: str, memory_gb: float, cpu_cores: int) -> bool:
        """Register a process with resource requirements."""
        if self.allocate_memory(memory_gb) and self.allocate_cpu(cpu_cores):
            self.registered_processes[process_id] = {
                'memory_gb': memory_gb,
                'cpu_cores': cpu_cores,
                'start_time': time.time()
            }
            self.active_processes[process_id] = self.registered_processes[process_id]
            return True
        return False
    
    def unregister_process(self, process_id: str) -> bool:
        """Unregister a process and free its resources."""
        if process_id in self.registered_processes:
            process_info = self.registered_processes[process_id]
            self.deallocate_memory(process_info['memory_gb'])
            self.deallocate_cpu(process_info['cpu_cores'])
            del self.registered_processes[process_id]
            if process_id in self.active_processes:
                del self.active_processes[process_id]
            return True
        return False
    
    def cleanup(self):
        """Cleanup resources and stop monitoring."""
        self._monitoring = False
        if self._monitor_thread and self._monitor_thread.is_alive():
            self._monitor_thread.join(timeout=1.0)
        
        # Clear all processes
        self.clear_all_processes()
        
        # Force garbage collection
        gc.collect()
        
        logger.info("ResourceManager cleanup completed")

class MemoryOptimizer:
    """Memory optimization utilities for large-scale processing."""
    
    @staticmethod
    def optimize_numpy_arrays(*arrays: np.ndarray) -> List[np.ndarray]:
        """Optimize numpy arrays for memory efficiency."""
        optimized = []
        
        for arr in arrays:
            # Convert to most efficient dtype
            if arr.dtype == np.float64:
                # Check if we can safely downcast to float32
                if np.allclose(arr, arr.astype(np.float32), rtol=1e-6):
                    arr = arr.astype(np.float32)
            
            # Make array contiguous for better cache performance
            if not arr.flags.c_contiguous:
                arr = np.ascontiguousarray(arr)
            
            optimized.append(arr)
        
        return optimized
    
    @staticmethod
    def batch_iterator(data: Any, batch_size: int, resource_manager: Optional[ResourceManager] = None):
        """Memory-efficient batch iterator with dynamic sizing."""
        total_size = len(data)
        current_pos = 0
        
        while current_pos < total_size:
            # Dynamically adjust batch size based on resources
            if resource_manager:
                batch_size = resource_manager.optimize_batch_size(batch_size, total_size - current_pos)
            
            end_pos = min(current_pos + batch_size, total_size)
            batch = data[current_pos:end_pos]
            
            yield batch
            current_pos = end_pos
            
            # Cleanup after each batch
            gc.collect()

class ScalabilityConfig:
    """Configuration for scalability parameters."""
    
    def __init__(self, 
                 max_proteins_per_batch: int = 1000,
                 max_memory_gb: float = 8.0,
                 enable_parallel_processing: bool = True,
                 max_workers: int = 4,
                 chunk_size: int = 1000,
                 use_streaming: bool = True,
                 cache_size_mb: int = 512):
        
        self.max_proteins_per_batch = max_proteins_per_batch
        self.max_memory_gb = max_memory_gb
        self.enable_parallel_processing = enable_parallel_processing
        self.max_workers = max_workers
        self.chunk_size = chunk_size
        self.use_streaming = use_streaming
        self.cache_size_mb = cache_size_mb
        
        # Auto-configure based on system resources
        self._auto_configure()
    
    def _auto_configure(self):
        """Auto-configure parameters based on available system resources."""
        memory = psutil.virtual_memory()
        cpu_count = psutil.cpu_count()
        
        # Adjust memory limits based on available RAM
        available_gb = memory.available / (1024**3)
        self.max_memory_gb = min(self.max_memory_gb, available_gb * 0.8)
        
        # Adjust worker count based on CPU cores
        self.max_workers = min(self.max_workers, max(1, cpu_count - 1))
        
        # Adjust batch size based on memory
        memory_factor = self.max_memory_gb / 8.0  # Base factor
        self.max_proteins_per_batch = int(self.max_proteins_per_batch * memory_factor)
        
        logger.info(f"Auto-configured scalability: memory={self.max_memory_gb:.1f}GB, "
                   f"workers={self.max_workers}, batch_size={self.max_proteins_per_batch}")

# Global resource manager instance
_resource_manager: Optional[ResourceManager] = None

def get_resource_manager() -> ResourceManager:
    """Get the global resource manager instance."""
    global _resource_manager
    if _resource_manager is None:
        _resource_manager = ResourceManager()
    return _resource_manager

def cleanup_resources():
    """Cleanup global resources."""
    global _resource_manager
    if _resource_manager is not None:
        _resource_manager.cleanup()
        _resource_manager = None
