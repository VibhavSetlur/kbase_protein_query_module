"""
Parallel Processing Module for Scalable Protein Analysis

This module provides comprehensive parallel processing capabilities optimized
for KBase environments with proper resource management and fault tolerance.
"""

import os
import time
import logging
import threading
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed, Future
from multiprocessing import Manager, Queue, Process
from typing import Dict, Any, List, Callable, Optional, Union, Iterator, Tuple
from dataclasses import dataclass
import numpy as np
from functools import partial
import psutil
import queue

from .resource_manager import ResourceManager, get_resource_manager

logger = logging.getLogger(__name__)

@dataclass
class ProcessingTask:
    """A task for parallel processing."""
    task_id: str
    function: Callable
    args: tuple
    kwargs: dict
    priority: int = 0
    retry_count: int = 0
    max_retries: int = 3

class ParallelProcessor:
    """
    High-performance parallel processor for protein analysis workflows.
    
    Features:
    - Adaptive thread/process pool management
    - Resource-aware task scheduling
    - Fault tolerance with automatic retries
    - Progress tracking and monitoring
    - Memory-efficient batch processing
    - Dynamic load balancing
    """
    
    def __init__(self, 
                 max_workers: Optional[int] = None,
                 use_processes: bool = False,
                 resource_manager: Optional[ResourceManager] = None,
                 enable_monitoring: bool = True):
        
        self.max_workers = max_workers or min(psutil.cpu_count(), 8)
        self.use_processes = use_processes
        self.resource_manager = resource_manager or get_resource_manager()
        self.enable_monitoring = enable_monitoring
        
        # Execution tracking
        self.active_tasks: Dict[str, ProcessingTask] = {}
        self.completed_tasks: Dict[str, Any] = {}
        self.failed_tasks: Dict[str, Exception] = {}
        
        # Performance metrics
        self.task_times: List[float] = []
        self.throughput_history: List[float] = []
        
        # Thread safety
        self._lock = threading.Lock()
        
        logger.info(f"ParallelProcessor initialized: workers={self.max_workers}, "
                   f"processes={use_processes}")
    
    def process_batch(self, 
                     tasks: List[ProcessingTask],
                     progress_callback: Optional[Callable[[int, int], None]] = None) -> Dict[str, Any]:
        """
        Process a batch of tasks in parallel with resource management.
        
        Args:
            tasks: List of processing tasks
            progress_callback: Optional progress callback function
            
        Returns:
            Dictionary of task_id -> result
        """
        if not tasks:
            return {}
        
        results = {}
        
        with self.resource_manager.resource_context("parallel_batch_processing"):
            # Choose executor based on configuration and resource availability
            executor_class = ProcessPoolExecutor if self.use_processes else ThreadPoolExecutor
            
            with executor_class(max_workers=self.max_workers) as executor:
                # Submit all tasks
                future_to_task = {}
                for task in tasks:
                    future = executor.submit(self._execute_task_safely, task)
                    future_to_task[future] = task
                    
                    with self._lock:
                        self.active_tasks[task.task_id] = task
                
                # Collect results as they complete
                completed_count = 0
                for future in as_completed(future_to_task):
                    task = future_to_task[future]
                    
                    try:
                        result = future.result()
                        results[task.task_id] = result
                        
                        with self._lock:
                            self.completed_tasks[task.task_id] = result
                            if task.task_id in self.active_tasks:
                                del self.active_tasks[task.task_id]
                        
                    except Exception as e:
                        logger.error(f"Task {task.task_id} failed: {e}")
                        results[task.task_id] = None
                        
                        with self._lock:
                            self.failed_tasks[task.task_id] = e
                            if task.task_id in self.active_tasks:
                                del self.active_tasks[task.task_id]
                    
                    completed_count += 1
                    
                    # Progress callback
                    if progress_callback:
                        progress_callback(completed_count, len(tasks))
        
        return results
    
    def _execute_task_safely(self, task: ProcessingTask) -> Any:
        """Execute a task with error handling and retries."""
        for attempt in range(task.max_retries + 1):
            try:
                start_time = time.time()
                result = task.function(*task.args, **task.kwargs)
                execution_time = time.time() - start_time
                
                # Track performance
                with self._lock:
                    self.task_times.append(execution_time)
                
                return result
                
            except Exception as e:
                if attempt < task.max_retries:
                    logger.warning(f"Task {task.task_id} attempt {attempt + 1} failed: {e}. Retrying...")
                    time.sleep(0.1 * (2 ** attempt))  # Exponential backoff
                else:
                    logger.error(f"Task {task.task_id} failed after {task.max_retries} retries: {e}")
                    raise
    
    def process_proteins_parallel(self,
                                proteins: List[Any],
                                process_function: Callable[[Any], Any],
                                batch_size: Optional[int] = None,
                                progress_callback: Optional[Callable[[int, int], None]] = None) -> List[Any]:
        """
        Process proteins in parallel with automatic batching and resource management.
        
        Args:
            proteins: List of proteins to process
            process_function: Function to apply to each protein
            batch_size: Optional batch size (auto-calculated if None)
            progress_callback: Optional progress callback
            
        Returns:
            List of processed results
        """
        if not proteins:
            return []
        
        # Optimize batch size based on resources
        if batch_size is None:
            batch_size = self.resource_manager.optimize_batch_size(
                base_batch_size=100, 
                data_size=len(proteins)
            )
        
        # Create tasks for parallel processing
        tasks = []
        for i, protein in enumerate(proteins):
            task = ProcessingTask(
                task_id=f"protein_{i}",
                function=process_function,
                args=(protein,),
                kwargs={}
            )
            tasks.append(task)
        
        # Process in batches to manage memory
        all_results = []
        total_processed = 0
        
        for batch_start in range(0, len(tasks), batch_size):
            batch_end = min(batch_start + batch_size, len(tasks))
            batch_tasks = tasks[batch_start:batch_end]
            
            logger.info(f"Processing batch {batch_start//batch_size + 1}: "
                       f"proteins {batch_start}-{batch_end-1}")
            
            # Process batch
            batch_results = self.process_batch(batch_tasks)
            
            # Collect results in order
            for task in batch_tasks:
                result = batch_results.get(task.task_id)
                all_results.append(result)
            
            total_processed += len(batch_tasks)
            
            # Progress callback for overall progress
            if progress_callback:
                progress_callback(total_processed, len(proteins))
            
            # Check resources and potentially pause
            if not self.resource_manager.check_resource_availability():
                logger.warning("Resource constraints detected, pausing briefly...")
                time.sleep(1.0)
        
        return all_results
    
    def submit_task(self, function: Callable, *args, **kwargs) -> str:
        """Submit a task for execution and return task ID."""
        # Generate a unique task ID
        task_id = f"task_{len(self.active_tasks) + len(self.completed_tasks) + len(self.failed_tasks)}"
        
        task = ProcessingTask(
            task_id=task_id,
            function=function,
            args=args,
            kwargs=kwargs
        )
        
        with self._lock:
            self.active_tasks[task_id] = task
        
        # Execute task immediately in a simple implementation
        try:
            result = self._execute_task_safely(task)
            with self._lock:
                self.completed_tasks[task_id] = result
                # Keep in active_tasks until explicitly removed
        except Exception as e:
            with self._lock:
                self.failed_tasks[task_id] = e
                # Keep in active_tasks until explicitly removed
        
        return task_id
    
    def get_task_result(self, task_id: str) -> Any:
        """Get the result of a completed task."""
        with self._lock:
            if task_id in self.completed_tasks:
                return self.completed_tasks[task_id]
            elif task_id in self.failed_tasks:
                return None
            else:
                return None
    
    def cancel_task(self, task_id: str) -> bool:
        """Cancel a task if it's still active."""
        with self._lock:
            if task_id in self.active_tasks:
                del self.active_tasks[task_id]
                return True
            return False
    
    def get_task_status(self, task_id: str) -> Optional[str]:
        """Get the status of a task."""
        with self._lock:
            if task_id in self.completed_tasks:
                return "completed"
            elif task_id in self.failed_tasks:
                return "failed"
            elif task_id in self.active_tasks:
                return "running"
            else:
                return None
    
    def list_active_tasks(self) -> List[str]:
        """List all active task IDs."""
        with self._lock:
            return list(self.active_tasks.keys())
    
    def list_completed_tasks(self) -> List[str]:
        """List all completed task IDs."""
        with self._lock:
            return list(self.completed_tasks.keys())
    
    def clear_completed_tasks(self) -> int:
        """Clear all completed tasks and return count."""
        with self._lock:
            count = len(self.completed_tasks)
            self.completed_tasks.clear()
            return count
    
    def get_performance_metrics(self) -> Dict[str, Any]:
        """Get performance metrics for the parallel processor."""
        with self._lock:
            if not self.task_times:
                return {}
            
            return {
                'total_tasks': len(self.task_times),
                'completed_tasks': len(self.completed_tasks),
                'failed_tasks': len(self.failed_tasks),
                'active_tasks': len(self.active_tasks),
                'avg_task_time': np.mean(self.task_times),
                'max_task_time': np.max(self.task_times),
                'min_task_time': np.min(self.task_times),
                'total_processing_time': np.sum(self.task_times),
                'tasks_per_second': len(self.task_times) / np.sum(self.task_times) if self.task_times else 0
            }

class StreamingProcessor:
    """
    Streaming processor for handling large datasets that don't fit in memory.
    """
    
    def __init__(self, 
                 chunk_size: int = 1000,
                 max_queue_size: int = 10,
                 resource_manager: Optional[ResourceManager] = None):
        
        self.chunk_size = chunk_size
        self.max_queue_size = max_queue_size
        self.resource_manager = resource_manager or get_resource_manager()
        
        # Streaming state
        self.input_queue = Queue(maxsize=max_queue_size)
        self.output_queue = Queue(maxsize=max_queue_size)
        self.processing_active = False
        
        logger.info(f"StreamingProcessor initialized: chunk_size={chunk_size}")
    
    def stream_process(self,
                      data_iterator: Iterator[Any],
                      process_function: Callable[[Any], Any],
                      output_handler: Callable[[Any], None]) -> None:
        """
        Stream process data with memory-efficient handling.
        
        Args:
            data_iterator: Iterator over input data
            process_function: Function to process each data item
            output_handler: Function to handle processed results
        """
        self.processing_active = True
        
        # Producer thread
        def producer():
            try:
                chunk = []
                for item in data_iterator:
                    chunk.append(item)
                    
                    if len(chunk) >= self.chunk_size:
                        self.input_queue.put(chunk)
                        chunk = []
                    
                    if not self.processing_active:
                        break
                
                # Put remaining items
                if chunk:
                    self.input_queue.put(chunk)
                
            except Exception as e:
                logger.error(f"Producer error: {e}")
            finally:
                self.input_queue.put(None)  # Sentinel
        
        # Consumer thread
        def consumer():
            try:
                while self.processing_active:
                    try:
                        result = self.output_queue.get(timeout=1.0)
                        if result is None:  # Sentinel
                            break
                        
                        output_handler(result)
                        
                    except queue.Empty:
                        continue
                    except Exception as e:
                        logger.error(f"Consumer error: {e}")
                        break
                        
            except Exception as e:
                logger.error(f"Consumer error: {e}")
        
        # Processor thread
        def processor():
            parallel_processor = ParallelProcessor(resource_manager=self.resource_manager)
            
            try:
                while self.processing_active:
                    try:
                        chunk = self.input_queue.get(timeout=1.0)
                        if chunk is None:  # Sentinel
                            break
                        
                        # Process chunk in parallel
                        results = parallel_processor.process_proteins_parallel(
                            chunk, process_function
                        )
                        
                        self.output_queue.put(results)
                        
                    except queue.Empty:
                        continue
                    except Exception as e:
                        logger.error(f"Processor error: {e}")
                        break
                        
            except Exception as e:
                logger.error(f"Processor error: {e}")
            finally:
                self.output_queue.put(None)  # Sentinel
        
        # Start threads
        producer_thread = threading.Thread(target=producer)
        consumer_thread = threading.Thread(target=consumer)
        processor_thread = threading.Thread(target=processor)
        
        producer_thread.start()
        processor_thread.start()
        consumer_thread.start()
        
        # Wait for completion
        producer_thread.join()
        processor_thread.join()
        consumer_thread.join()
        
        self.processing_active = False

class BatchQueue:
    """
    Intelligent batch queue for managing large-scale processing jobs.
    """
    
    def __init__(self, max_batch_size: int = 1000, priority_levels: int = 3):
        self.max_batch_size = max_batch_size
        self.priority_queues = [[] for _ in range(priority_levels)]
        self.batch_history: List[Dict[str, Any]] = []
        self._lock = threading.Lock()
    
    def add_task(self, task: ProcessingTask):
        """Add a task to the appropriate priority queue."""
        with self._lock:
            priority = min(task.priority, len(self.priority_queues) - 1)
            self.priority_queues[priority].append(task)
    
    def get_next_batch(self) -> List[ProcessingTask]:
        """Get the next batch of tasks for processing."""
        with self._lock:
            batch = []
            
            # Process high priority tasks first
            for priority_queue in self.priority_queues:
                while len(batch) < self.max_batch_size and priority_queue:
                    batch.append(priority_queue.pop(0))
            
            if batch:
                batch_info = {
                    'size': len(batch),
                    'timestamp': time.time(),
                    'priority_distribution': {}
                }
                
                for task in batch:
                    priority = task.priority
                    batch_info['priority_distribution'][priority] = \
                        batch_info['priority_distribution'].get(priority, 0) + 1
                
                self.batch_history.append(batch_info)
            
            return batch
    
    def is_empty(self) -> bool:
        """Check if all queues are empty."""
        with self._lock:
            return all(len(q) == 0 for q in self.priority_queues)
    
    def get_queue_status(self) -> Dict[str, Any]:
        """Get current queue status."""
        with self._lock:
            return {
                'queue_lengths': [len(q) for q in self.priority_queues],
                'total_pending': sum(len(q) for q in self.priority_queues),
                'batches_processed': len(self.batch_history)
            }

# Global parallel processor instance
_parallel_processor: Optional[ParallelProcessor] = None

def get_parallel_processor() -> ParallelProcessor:
    """Get the global parallel processor instance."""
    global _parallel_processor
    if _parallel_processor is None:
        _parallel_processor = ParallelProcessor()
    return _parallel_processor
