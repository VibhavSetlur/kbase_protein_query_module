"""
Unit tests for ParallelProcessor
"""
import pytest
from unittest.mock import Mock, patch
from lib.kbase_protein_query_module.src.core.parallel_processor import ParallelProcessor


class TestParallelProcessor:
    """Test cases for ParallelProcessor"""
    
    def test_parallel_processor_initialization(self):
        """Test ParallelProcessor initializes correctly"""
        processor = ParallelProcessor(max_workers=4)
        assert processor.max_workers == 4
        assert processor.active_tasks == {}
        assert processor.completed_tasks == {}
    
    def test_parallel_processor_default_workers(self):
        """Test ParallelProcessor with default worker count"""
        processor = ParallelProcessor()
        assert processor.max_workers > 0
    
    def test_submit_task(self):
        """Test submitting a task"""
        processor = ParallelProcessor(max_workers=2)
        
        def dummy_task(x):
            return x * 2
        
        task_id = processor.submit_task(dummy_task, 5)
        assert task_id is not None
        assert task_id in processor.active_tasks
    
    def test_get_task_result(self):
        """Test getting task result"""
        processor = ParallelProcessor(max_workers=2)
        
        def dummy_task(x):
            return x * 2
        
        task_id = processor.submit_task(dummy_task, 5)
        
        # Wait a bit for task to complete
        import time
        time.sleep(0.1)
        
        result = processor.get_task_result(task_id)
        assert result == 10
    
    def test_get_task_result_not_found(self):
        """Test getting result for non-existent task"""
        processor = ParallelProcessor()
        result = processor.get_task_result("non_existent")
        assert result is None
    
    def test_cancel_task(self):
        """Test canceling a task"""
        processor = ParallelProcessor(max_workers=2)
        
        def long_task():
            import time
            time.sleep(1)
            return "completed"
        
        task_id = processor.submit_task(long_task)
        success = processor.cancel_task(task_id)
        assert success is True
        assert task_id not in processor.active_tasks
    
    def test_cancel_task_not_found(self):
        """Test canceling non-existent task"""
        processor = ParallelProcessor()
        success = processor.cancel_task("non_existent")
        assert success is False
    
    def test_get_task_status(self):
        """Test getting task status"""
        processor = ParallelProcessor(max_workers=2)
        
        def dummy_task(x):
            return x * 2
        
        task_id = processor.submit_task(dummy_task, 5)
        status = processor.get_task_status(task_id)
        assert status in ["pending", "running", "completed", "failed"]
    
    def test_get_task_status_not_found(self):
        """Test getting status for non-existent task"""
        processor = ParallelProcessor()
        status = processor.get_task_status("non_existent")
        assert status is None
    
    def test_list_active_tasks(self):
        """Test listing active tasks"""
        processor = ParallelProcessor(max_workers=2)
        
        def dummy_task(x):
            return x * 2
        
        task_id1 = processor.submit_task(dummy_task, 5)
        task_id2 = processor.submit_task(dummy_task, 10)
        
        active_tasks = processor.list_active_tasks()
        assert len(active_tasks) == 2
        assert task_id1 in active_tasks
        assert task_id2 in active_tasks
    
    def test_list_completed_tasks(self):
        """Test listing completed tasks"""
        processor = ParallelProcessor(max_workers=2)
        
        def dummy_task(x):
            return x * 2
        
        task_id = processor.submit_task(dummy_task, 5)
        
        # Wait for task to complete
        import time
        time.sleep(0.1)
        
        completed_tasks = processor.list_completed_tasks()
        assert len(completed_tasks) >= 1
        assert task_id in completed_tasks
    
    def test_clear_completed_tasks(self):
        """Test clearing completed tasks"""
        processor = ParallelProcessor(max_workers=2)
        
        def dummy_task(x):
            return x * 2
        
        task_id = processor.submit_task(dummy_task, 5)
        
        # Wait for task to complete
        import time
        time.sleep(0.1)
        
        assert len(processor.completed_tasks) >= 1
        processor.clear_completed_tasks()
        assert len(processor.completed_tasks) == 0
    
    def test_task_with_exception(self):
        """Test task that raises exception"""
        processor = ParallelProcessor(max_workers=2)
        
        def failing_task():
            raise ValueError("Task failed")
        
        task_id = processor.submit_task(failing_task)
        
        # Wait for task to complete
        import time
        time.sleep(0.1)
        
        result = processor.get_task_result(task_id)
        assert result is None
        
        status = processor.get_task_status(task_id)
        assert status == "failed"
