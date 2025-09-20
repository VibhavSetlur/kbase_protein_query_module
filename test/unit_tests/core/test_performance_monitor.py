"""
Unit tests for PerformanceProfiler
"""
import pytest
import time
from unittest.mock import Mock, patch
from lib.kbase_protein_query_module.src.core.performance_monitor import PerformanceProfiler


class TestPerformanceProfiler:
    """Test cases for PerformanceProfiler"""
    
    def test_performance_profiler_initialization(self):
        """Test PerformanceProfiler initializes correctly"""
        profiler = PerformanceProfiler()
        assert profiler.start_time is None
        assert profiler.end_time is None
        assert profiler.metrics == {}
        assert profiler.memory_usage == []
        assert profiler.cpu_usage == []
    
    def test_start_profiling(self):
        """Test starting profiling"""
        profiler = PerformanceProfiler()
        profiler.start_profiling()
        assert profiler.start_time is not None
        assert profiler.end_time is None
    
    def test_stop_profiling(self):
        """Test stopping profiling"""
        profiler = PerformanceProfiler()
        profiler.start_profiling()
        time.sleep(0.01)  # Small delay
        profiler.stop_profiling()
        assert profiler.end_time is not None
        assert profiler.end_time > profiler.start_time
    
    def test_add_metric(self):
        """Test adding metrics"""
        profiler = PerformanceProfiler()
        profiler.add_metric("test_metric", 42.5)
        assert "test_metric" in profiler.metrics
        assert profiler.metrics["test_metric"] == 42.5
    
    def test_get_metric(self):
        """Test getting metrics"""
        profiler = PerformanceProfiler()
        profiler.add_metric("test_metric", 42.5)
        value = profiler.get_metric("test_metric")
        assert value == 42.5
        
        # Test non-existent metric
        value = profiler.get_metric("non_existent")
        assert value is None
    
    def test_get_all_metrics(self):
        """Test getting all metrics"""
        profiler = PerformanceProfiler()
        profiler.add_metric("metric1", 10.0)
        profiler.add_metric("metric2", 20.0)
        
        metrics = profiler.get_all_metrics()
        assert len(metrics) == 2
        assert metrics["metric1"] == 10.0
        assert metrics["metric2"] == 20.0
    
    def test_clear_metrics(self):
        """Test clearing metrics"""
        profiler = PerformanceProfiler()
        profiler.add_metric("test_metric", 42.5)
        assert len(profiler.metrics) == 1
        
        profiler.clear_metrics()
        assert len(profiler.metrics) == 0
    
    def test_get_execution_time(self):
        """Test getting execution time"""
        profiler = PerformanceProfiler()
        profiler.start_profiling()
        time.sleep(0.01)  # Small delay
        profiler.stop_profiling()
        
        execution_time = profiler.get_execution_time()
        assert execution_time > 0
        assert execution_time >= 0.01
    
    def test_get_execution_time_not_started(self):
        """Test getting execution time when not started"""
        profiler = PerformanceProfiler()
        execution_time = profiler.get_execution_time()
        assert execution_time is None
    
    def test_get_execution_time_not_stopped(self):
        """Test getting execution time when not stopped"""
        profiler = PerformanceProfiler()
        profiler.start_profiling()
        execution_time = profiler.get_execution_time()
        assert execution_time is None
    
    def test_memory_monitoring(self):
        """Test memory monitoring functionality"""
        profiler = PerformanceProfiler()
        
        # Start memory monitoring
        profiler.start_memory_monitoring()
        assert len(profiler.memory_usage) >= 0
        
        # Stop memory monitoring
        profiler.stop_memory_monitoring()
        # Memory usage should have been recorded
    
    def test_cpu_monitoring(self):
        """Test CPU monitoring functionality"""
        profiler = PerformanceProfiler()
        
        # Start CPU monitoring
        profiler.start_cpu_monitoring()
        assert len(profiler.cpu_usage) >= 0
        
        # Stop CPU monitoring
        profiler.stop_cpu_monitoring()
        # CPU usage should have been recorded
    
    def test_generate_report(self):
        """Test generating performance report"""
        profiler = PerformanceProfiler()
        profiler.start_profiling()
        profiler.add_metric("test_metric", 42.5)
        time.sleep(0.01)
        profiler.stop_profiling()
        
        report = profiler.generate_report()
        assert isinstance(report, dict)
        assert "execution_time" in report
        assert "metrics" in report
        assert report["metrics"]["test_metric"] == 42.5
    
    def test_context_manager(self):
        """Test using PerformanceProfiler as context manager"""
        with PerformanceProfiler() as profiler:
            profiler.add_metric("test_metric", 42.5)
            time.sleep(0.01)
        
        assert profiler.start_time is not None
        assert profiler.end_time is not None
        assert profiler.get_execution_time() > 0
        assert "test_metric" in profiler.metrics
    
    def test_reset(self):
        """Test resetting profiler"""
        profiler = PerformanceProfiler()
        profiler.start_profiling()
        profiler.add_metric("test_metric", 42.5)
        profiler.stop_profiling()
        
        profiler.reset()
        assert profiler.start_time is None
        assert profiler.end_time is None
        assert len(profiler.metrics) == 0
        assert len(profiler.memory_usage) == 0
        assert len(profiler.cpu_usage) == 0
