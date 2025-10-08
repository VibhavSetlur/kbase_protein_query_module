"""
Unit tests for KBase server module.

Tests the server implementation, request handling, and KBase integration.
"""

import pytest
import sys
import os
from unittest.mock import Mock, patch, MagicMock

# Add lib to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../lib'))

from kbase_protein_query_module.kbase_protein_query_moduleServer import Application, MethodContext, JSONRPCServiceCustom


class TestKBaseServerModule:
    """Test cases for KBase server module."""
    
    def test_application_initialization(self):
        """Test Application class initialization."""
        app = Application()
        
        assert app.userlog is not None
        assert app.serverlog is not None
        assert app.rpc_service is not None
        assert isinstance(app.rpc_service, JSONRPCServiceCustom)
        assert app.method_authentication is not None
        assert app.auth_client is not None
    
    def test_application_method_registration(self):
        """Test that required methods are registered."""
        app = Application()
        
        # Check that required methods are registered
        required_methods = [
            'kbase_protein_query_module.run_protein_query_analysis',
            'kbase_protein_query_module.get_available_analyses',
            'kbase_protein_query_module.status'
        ]
        
        for method in required_methods:
            assert method in app.method_authentication
            assert app.method_authentication[method] == 'required'
    
    def test_method_context_initialization(self):
        """Test MethodContext initialization."""
        mock_logger = Mock()
        ctx = MethodContext(mock_logger)
        
        assert ctx['client_ip'] is None
        assert ctx['user_id'] is None
        assert ctx['authenticated'] is None
        assert ctx['token'] is None
        assert ctx['module'] is None
        assert ctx['method'] is None
        assert ctx['call_id'] is None
        assert ctx['rpc_context'] is None
        assert ctx['provenance'] is None
        assert ctx._logger == mock_logger
    
    def test_method_context_logging(self):
        """Test MethodContext logging methods."""
        mock_logger = Mock()
        ctx = MethodContext(mock_logger)
        
        # Test log_err
        ctx.log_err("Error message")
        mock_logger.log_message.assert_called()
        
        # Test log_info
        ctx.log_info("Info message")
        mock_logger.log_message.assert_called()
        
        # Test log_debug
        ctx.log_debug("Debug message")
        mock_logger.log_message.assert_called()
    
    def test_method_context_provenance(self):
        """Test MethodContext provenance handling."""
        mock_logger = Mock()
        ctx = MethodContext(mock_logger)
        
        # Test with callback URL
        with patch.dict('os.environ', {'SDK_CALLBACK_URL': 'http://callback.test'}):
            with patch('kbase_protein_query_module.kbase_protein_query_moduleServer._requests.post') as mock_post:
                mock_response = Mock()
                mock_response.status_code = 200
                mock_response.ok = True
                mock_response.json.return_value = {'result': [{'test': 'provenance'}]}
                mock_post.return_value = mock_response
                
                provenance = ctx.provenance()
                assert provenance == {'test': 'provenance'}
                mock_post.assert_called_once()
        
        # Test without callback URL
        with patch.dict('os.environ', {}, clear=True):
            ctx['provenance'] = [{'test': 'local_provenance'}]
            provenance = ctx.provenance()
            assert provenance == [{'test': 'local_provenance'}]
    
    def test_jsonrpc_service_call(self):
        """Test JSONRPCServiceCustom call method."""
        service = JSONRPCServiceCustom()
        
        # Mock a simple method
        def test_method(ctx, param):
            return {'result': f'processed_{param}'}
        
        service.add(test_method, name='test.method', types=[str])
        
        # Test call with valid JSON
        ctx = MethodContext(Mock())
        jsondata = {
            'method': 'test.method',
            'params': ['test_param'],
            'id': 1,
            'jsonrpc': '2.0'
        }
        
        result = service.call(ctx, jsondata)
        
        assert result is not None
        import json
        parsed_result = json.loads(result)
        assert parsed_result['result']['result'] == 'processed_test_param'
        assert parsed_result['id'] == 1
    
    def test_jsonrpc_service_call_py(self):
        """Test JSONRPCServiceCustom call_py method."""
        service = JSONRPCServiceCustom()
        
        # Mock a simple method
        def test_method(ctx, param):
            return {'result': f'processed_{param}'}
        
        service.add(test_method, name='test.method', types=[str])
        
        # Test call_py with valid JSON
        ctx = MethodContext(Mock())
        jsondata = {
            'method': 'test.method',
            'params': ['test_param'],
            'id': 1,
            'jsonrpc': '2.0'
        }
        
        result = service.call_py(ctx, jsondata)
        
        assert result is not None
        assert result['result']['result'] == 'processed_test_param'
        assert result['id'] == 1
    
    def test_jsonrpc_service_batch_call(self):
        """Test JSONRPCServiceCustom batch call handling."""
        service = JSONRPCServiceCustom()
        
        # Mock methods
        def test_method1(ctx, param):
            return {'result': f'method1_{param}'}
        
        def test_method2(ctx, param):
            return {'result': f'method2_{param}'}
        
        service.add(test_method1, name='test.method1', types=[str])
        service.add(test_method2, name='test.method2', types=[str])
        
        # Test batch call
        ctx = MethodContext(Mock())
        jsondata = [
            {
                'method': 'test.method1',
                'params': ['param1'],
                'id': 1,
                'jsonrpc': '2.0'
            },
            {
                'method': 'test.method2',
                'params': ['param2'],
                'id': 2,
                'jsonrpc': '2.0'
            }
        ]
        
        result = service.call_py(ctx, jsondata)
        
        assert isinstance(result, list)
        assert len(result) == 2
        assert result[0]['result']['result'] == 'method1_param1'
        assert result[1]['result']['result'] == 'method2_param2'
    
    def test_jsonrpc_service_error_handling(self):
        """Test JSONRPCServiceCustom error handling."""
        service = JSONRPCServiceCustom()
        
        # Mock a method that raises an exception
        def error_method(ctx):
            raise ValueError("Test error")
        
        service.add(error_method, name='test.error')
        
        # Test error handling
        ctx = MethodContext(Mock())
        jsondata = {
            'method': 'test.error',
            'params': [],
            'id': 1,
            'jsonrpc': '2.0'
        }
        
        result = service.call_py(ctx, jsondata)
        
        assert 'error' in result
        assert result['error']['code'] == 0  # Unexpected Server Error
        assert result['error']['name'] == 'Unexpected Server Error'
        assert 'Test error' in result['error']['message']
    
    def test_application_call_method(self):
        """Test Application __call__ method (WSGI interface)."""
        app = Application()
        
        # Mock WSGI environment
        environ = {
            'REQUEST_METHOD': 'POST',
            'CONTENT_LENGTH': '100',
            'HTTP_AUTHORIZATION': 'Bearer test_token'
        }
        
        # Mock start_response
        from unittest.mock import Mock as _Mock
        start_response = _Mock()
        
        # Mock wsgi.input
        mock_input = _Mock()
        mock_input.read.return_value = b'{"method": "kbase_protein_query_module.status", "params": [], "id": 1, "jsonrpc": "2.0"}'
        
        with patch.dict('os.environ', {}, clear=True):
            # Mock the auth client
            app.auth_client.get_user.return_value = 'test_user'

            # Mock the RPC service
            from unittest.mock import Mock
            app.rpc_service.call = Mock(return_value='{"result": {"state": "OK"}, "id": 1}')

            # Test the call
            result = app(environ, start_response)

            # Verify start_response was called with success status and headers list
            args, _ = start_response.call_args
            assert args[0] == '200 OK'
            assert isinstance(args[1], list)
            assert all(isinstance(h, tuple) and len(h) == 2 for h in args[1])

            # Verify result is returned
            assert isinstance(result, list)
            assert len(result) == 1
    
    def test_application_authentication_required(self):
        """Test Application authentication requirement handling."""
        app = Application()
        
        # Mock WSGI environment without authorization
        environ = {
            'REQUEST_METHOD': 'POST',
            'CONTENT_LENGTH': '100'
            # No HTTP_AUTHORIZATION header
        }
        
        start_response = Mock()
        
        # Mock wsgi.input for a method that requires authentication
        mock_input = Mock()
        mock_input.read.return_value = b'{"method": "kbase_protein_query_module.run_protein_query_analysis", "params": [{}], "id": 1, "jsonrpc": "2.0"}'
        
        with patch.dict('os.environ', {}, clear=True):
            result = app(environ, start_response)

            # Should return error response; validate status and headers shape
            args, _ = start_response.call_args
            assert args[0] == '500 Internal Server Error'
            assert isinstance(args[1], list)
    
    def test_application_options_request(self):
        """Test Application handling of OPTIONS request."""
        app = Application()
        
        # Mock WSGI environment for OPTIONS request
        environ = {
            'REQUEST_METHOD': 'OPTIONS'
        }
        
        start_response = Mock()
        
        result = app(environ, start_response)

        # Should return 200 OK for OPTIONS
        args, _ = start_response.call_args
        assert args[0] == '200 OK'
        assert isinstance(args[1], list)
    
    def test_application_error_handling(self):
        """Test Application error handling."""
        app = Application()
        
        # Mock WSGI environment
        environ = {
            'REQUEST_METHOD': 'POST',
            'CONTENT_LENGTH': '100'
        }
        
        start_response = Mock()
        
        # Mock wsgi.input with invalid JSON
        mock_input = Mock()
        mock_input.read.return_value = b'invalid json'
        
        with patch.dict('os.environ', {}, clear=True):
            result = app(environ, start_response)

            # Should handle JSON parsing error
            args, _ = start_response.call_args
            assert args[0] == '500 Internal Server Error'
            assert isinstance(args[1], list)
    
    def test_get_ip_address(self):
        """Test IP address extraction."""
        from kbase_protein_query_module.kbase_protein_query_moduleServer import getIPAddress
        
        # Test with X-Forwarded-For header
        environ = {
            'HTTP_X_FORWARDED_FOR': '192.168.1.1, 10.0.0.1',
            'REMOTE_ADDR': '127.0.0.1'
        }
        
        with patch('kbase_protein_query_module.kbase_protein_query_moduleServer.config', None):
            ip = getIPAddress(environ)
            assert ip == '192.168.1.1'  # Should take first IP
        
        # Test with X-Real-IP header
        environ = {
            'HTTP_X_REAL_IP': '203.0.113.1',
            'REMOTE_ADDR': '127.0.0.1'
        }
        
        with patch('kbase_protein_query_module.kbase_protein_query_moduleServer.config', None):
            ip = getIPAddress(environ)
            assert ip == '203.0.113.1'
        
        # Test with only REMOTE_ADDR
        environ = {
            'REMOTE_ADDR': '127.0.0.1'
        }
        
        with patch('kbase_protein_query_module.kbase_protein_query_moduleServer.config', None):
            ip = getIPAddress(environ)
            assert ip == '127.0.0.1'
    
    def test_server_error_class(self):
        """Test ServerError class."""
        from kbase_protein_query_module.kbase_protein_query_moduleServer import ServerError
        
        # Test ServerError initialization
        error = ServerError('TestError', 100, 'Test error message', 'Stack trace')
        
        assert error.name == 'TestError'
        assert error.code == 100
        assert error.message == 'Test error message'
        assert error.data == 'Stack trace'
        
        # Test string representation
        error_str = str(error)
        assert 'TestError' in error_str
        assert '100' in error_str
        assert 'Test error message' in error_str
        assert 'Stack trace' in error_str
    
    def test_json_object_encoder(self):
        """Test JSONObjectEncoder class."""
        from kbase_protein_query_module.kbase_protein_query_moduleServer import JSONObjectEncoder
        import json
        
        encoder = JSONObjectEncoder()
        
        # Test with set
        test_set = {1, 2, 3}
        result = encoder.default(test_set)
        assert isinstance(result, list)
        assert set(result) == test_set
        
        # Test with frozenset
        test_frozenset = frozenset([1, 2, 3])
        result = encoder.default(test_frozenset)
        assert isinstance(result, list)
        assert set(result) == test_frozenset
        
        # Test with object having toJSONable method
        class MockObject:
            def toJSONable(self):
                return {'mock': 'object'}
        
        mock_obj = MockObject()
        result = encoder.default(mock_obj)
        assert result == {'mock': 'object'}
        
        # Test with regular object (should use default encoder)
        regular_obj = "test string"
        result = encoder.default(regular_obj)
        assert result == "test string"

