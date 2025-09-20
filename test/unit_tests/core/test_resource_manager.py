"""
Unit tests for ResourceManager
"""
import pytest
from unittest.mock import Mock, patch
from lib.kbase_protein_query_module.src.core.resource_manager import ResourceManager


class TestResourceManager:
    """Test cases for ResourceManager"""
    
    def test_resource_manager_initialization(self):
        """Test ResourceManager initializes correctly"""
        manager = ResourceManager()
        assert manager.max_memory_gb > 0
        assert manager.max_cpu_cores > 0
        assert manager.allocated_memory == 0
        assert manager.allocated_cpu == 0
        assert manager.active_processes == {}
    
    def test_resource_manager_custom_limits(self):
        """Test ResourceManager with custom limits"""
        manager = ResourceManager(max_memory_gb=16, max_cpu_cores=8)
        assert manager.max_memory_gb == 16
        assert manager.max_cpu_cores == 8
    
    def test_check_memory_availability(self):
        """Test checking memory availability"""
        manager = ResourceManager(max_memory_gb=8)
        
        # Test available memory
        available = manager.check_memory_availability(4)
        assert available is True
        
        # Test insufficient memory
        available = manager.check_memory_availability(10)
        assert available is False
    
    def test_check_cpu_availability(self):
        """Test checking CPU availability"""
        manager = ResourceManager(max_cpu_cores=4)
        
        # Test available CPU
        available = manager.check_cpu_availability(2)
        assert available is True
        
        # Test insufficient CPU
        available = manager.check_cpu_availability(6)
        assert available is False
    
    def test_allocate_memory(self):
        """Test allocating memory"""
        manager = ResourceManager(max_memory_gb=8)
        
        success = manager.allocate_memory(4)
        assert success is True
        assert manager.allocated_memory == 4
        
        # Test allocating more than available
        success = manager.allocate_memory(6)
        assert success is False
        assert manager.allocated_memory == 4  # Should remain unchanged
    
    def test_deallocate_memory(self):
        """Test deallocating memory"""
        manager = ResourceManager(max_memory_gb=8)
        manager.allocated_memory = 4
        
        manager.deallocate_memory(2)
        assert manager.allocated_memory == 2
        
        # Test deallocating more than allocated
        manager.deallocate_memory(4)
        assert manager.allocated_memory == 0  # Should not go negative
    
    def test_allocate_cpu(self):
        """Test allocating CPU cores"""
        manager = ResourceManager(max_cpu_cores=4)
        
        success = manager.allocate_cpu(2)
        assert success is True
        assert manager.allocated_cpu == 2
        
        # Test allocating more than available
        success = manager.allocate_cpu(4)
        assert success is False
        assert manager.allocated_cpu == 2  # Should remain unchanged
    
    def test_deallocate_cpu(self):
        """Test deallocating CPU cores"""
        manager = ResourceManager(max_cpu_cores=4)
        manager.allocated_cpu = 3
        
        manager.deallocate_cpu(1)
        assert manager.allocated_cpu == 2
        
        # Test deallocating more than allocated
        manager.deallocate_cpu(4)
        assert manager.allocated_cpu == 0  # Should not go negative
    
    def test_register_process(self):
        """Test registering a process"""
        manager = ResourceManager()
        
        process_id = "test_process"
        memory_gb = 2
        cpu_cores = 1
        
        success = manager.register_process(process_id, memory_gb, cpu_cores)
        assert success is True
        assert process_id in manager.active_processes
        assert manager.active_processes[process_id]["memory_gb"] == memory_gb
        assert manager.active_processes[process_id]["cpu_cores"] == cpu_cores
    
    def test_register_process_insufficient_resources(self):
        """Test registering process with insufficient resources"""
        manager = ResourceManager(max_memory_gb=1, max_cpu_cores=1)
        
        process_id = "test_process"
        memory_gb = 2  # More than available
        cpu_cores = 1
        
        success = manager.register_process(process_id, memory_gb, cpu_cores)
        assert success is False
        assert process_id not in manager.active_processes
    
    def test_unregister_process(self):
        """Test unregistering a process"""
        manager = ResourceManager()
        
        process_id = "test_process"
        memory_gb = 2
        cpu_cores = 1
        
        manager.register_process(process_id, memory_gb, cpu_cores)
        assert process_id in manager.active_processes
        
        manager.unregister_process(process_id)
        assert process_id not in manager.active_processes
        assert manager.allocated_memory == 0
        assert manager.allocated_cpu == 0
    
    def test_unregister_process_not_found(self):
        """Test unregistering non-existent process"""
        manager = ResourceManager()
        
        # Should not raise exception
        manager.unregister_process("non_existent")
    
    def test_get_resource_usage(self):
        """Test getting resource usage"""
        manager = ResourceManager(max_memory_gb=8, max_cpu_cores=4)
        
        manager.register_process("process1", 2, 1)
        manager.register_process("process2", 3, 2)
        
        usage = manager.get_resource_usage()
        assert usage["memory_gb"] == 5
        assert usage["cpu_cores"] == 3
        assert usage["memory_percentage"] == 62.5  # 5/8 * 100
        assert usage["cpu_percentage"] == 75.0  # 3/4 * 100
    
    def test_get_available_resources(self):
        """Test getting available resources"""
        manager = ResourceManager(max_memory_gb=8, max_cpu_cores=4)
        
        manager.register_process("process1", 2, 1)
        
        available = manager.get_available_resources()
        assert available["memory_gb"] == 6
        assert available["cpu_cores"] == 3
    
    def test_list_active_processes(self):
        """Test listing active processes"""
        manager = ResourceManager()
        
        manager.register_process("process1", 2, 1)
        manager.register_process("process2", 3, 2)
        
        processes = manager.list_active_processes()
        assert len(processes) == 2
        assert "process1" in processes
        assert "process2" in processes
    
    def test_clear_all_processes(self):
        """Test clearing all processes"""
        manager = ResourceManager()
        
        manager.register_process("process1", 2, 1)
        manager.register_process("process2", 3, 2)
        
        assert len(manager.active_processes) == 2
        assert manager.allocated_memory == 5
        assert manager.allocated_cpu == 3
        
        manager.clear_all_processes()
        
        assert len(manager.active_processes) == 0
        assert manager.allocated_memory == 0
        assert manager.allocated_cpu == 0
