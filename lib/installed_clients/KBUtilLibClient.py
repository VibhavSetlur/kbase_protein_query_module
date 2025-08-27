#!/usr/bin/env python
"""
KBUtilLib Client for KBase SDK Integration

This client provides access to KBUtilLib utilities within the KBase SDK framework.
"""

import os
import sys
import logging
from typing import Dict, Any, Optional, List

# Add kbutillib to path
kbutillib_path = os.path.join(os.path.dirname(__file__), 'kbutillib')
if kbutillib_path not in sys.path:
    sys.path.insert(0, kbutillib_path)

# Also add parent directory to ensure imports work
parent_path = os.path.dirname(__file__)
if parent_path not in sys.path:
    sys.path.insert(0, parent_path)

try:
    from kbutillib import BaseUtils, KBWSUtils, KBSDKUtils
except ImportError as e:
    logging.warning(f"Could not import KBUtilLib modules: {e}")
    BaseUtils = None
    KBWSUtils = None
    KBSDKUtils = None

logger = logging.getLogger(__name__)


class KBUtilLib:
    """
    KBUtilLib client for KBase SDK integration.
    
    Provides access to KBUtilLib utilities including:
    - Workspace operations (KBWSUtils)
    - SDK utilities (KBSDKUtils) 
    - Base utilities (BaseUtils)
    - File operations and data handling
    """
    
    def __init__(self, callback_url: str, token: str = None, **kwargs):
        """
        Initialize KBUtilLib client.
        
        Args:
            callback_url: KBase callback URL
            token: Authentication token
            **kwargs: Additional parameters
        """
        self.callback_url = callback_url
        self.token = token or os.environ.get('KB_AUTH_TOKEN')
        
        # Initialize utility classes if available
        self.base_utils = None
        self.ws_utils = None
        self.sdk_utils = None
        
        self._initialize_utils(**kwargs)
    
    def _initialize_utils(self, **kwargs):
        """Initialize available utility classes."""
        try:
            if BaseUtils:
                self.base_utils = BaseUtils(**kwargs)
            
            if KBWSUtils:
                # Set up workspace utilities with KBase environment
                ws_kwargs = {
                    'kb_version': 'prod',
                    'max_retry': 3,
                    **kwargs
                }
                self.ws_utils = KBWSUtils(**ws_kwargs)
            
            if KBSDKUtils:
                # Set up SDK utilities 
                sdk_kwargs = {
                    'callback_url': self.callback_url,
                    'working_dir': kwargs.get('scratch', '/tmp'),
                    **kwargs
                }
                self.sdk_utils = KBSDKUtils(**sdk_kwargs)
                
            logger.info("KBUtilLib utilities initialized successfully")
            
        except Exception as e:
            logger.warning(f"Failed to initialize some KBUtilLib utilities: {e}")
    
    # Workspace utility methods
    def get_workspace_client(self):
        """Get workspace client from KBUtilLib."""
        if self.ws_utils:
            return self.ws_utils.ws_client()
        return None
    
    def set_workspace(self, workspace):
        """Set workspace using KBUtilLib utilities."""
        if self.ws_utils:
            return self.ws_utils.set_ws(workspace)
        return None
    
    def save_workspace_object(self, workspace_name: str, object_name: str, 
                            object_type: str, object_data: Dict[str, Any]) -> str:
        """
        Save object to workspace using KBUtilLib patterns.
        
        Args:
            workspace_name: Name of workspace
            object_name: Name for the object
            object_type: KBase object type
            object_data: Object data to save
            
        Returns:
            Object reference string
        """
        if not self.ws_utils:
            raise RuntimeError("Workspace utilities not available")
        
        try:
            self.ws_utils.set_ws(workspace_name)
            ws_client = self.ws_utils.ws_client()
            
            result = ws_client.save_objects({
                'id': workspace_name,
                'objects': [{
                    'name': object_name,
                    'type': object_type,
                    'data': object_data
                }]
            })
            
            return f"{workspace_name}/{object_name}"
            
        except Exception as e:
            logger.error(f"Failed to save workspace object: {e}")
            raise
    
    def get_workspace_object(self, object_ref: str) -> Dict[str, Any]:
        """
        Get object from workspace using KBUtilLib patterns.
        
        Args:
            object_ref: Workspace object reference
            
        Returns:
            Object data
        """
        if not self.ws_utils:
            raise RuntimeError("Workspace utilities not available")
        
        try:
            ws_client = self.ws_utils.ws_client()
            result = ws_client.get_objects2({
                'objects': [{'ref': object_ref}]
            })
            
            return result['data'][0]['data']
            
        except Exception as e:
            logger.error(f"Failed to get workspace object: {e}")
            raise
    
    # File utility methods
    def create_report(self, report_data: Dict[str, Any], workspace_name: str, 
                     report_name: str = None) -> Dict[str, Any]:
        """
        Create KBase report using KBUtilLib patterns.
        
        Args:
            report_data: Report data including message, objects_created, etc.
            workspace_name: Workspace for the report
            report_name: Name for the report object
            
        Returns:
            Report information
        """
        try:
            # Use KBase report patterns from KBUtilLib
            from installed_clients.KBaseReportClient import KBaseReport
            
            report_client = KBaseReport(self.callback_url)
            
            report_params = {
                'message': report_data.get('message', 'Analysis completed'),
                'workspace_name': workspace_name,
                'report_object_name': report_name or 'analysis_report'
            }
            
            # Add optional parameters
            if 'objects_created' in report_data:
                report_params['objects_created'] = report_data['objects_created']
            if 'html_links' in report_data:
                report_params['html_links'] = report_data['html_links']
            if 'file_links' in report_data:
                report_params['file_links'] = report_data['file_links']
            
            return report_client.create_extended_report(report_params)
            
        except Exception as e:
            logger.error(f"Failed to create report: {e}")
            raise
    
    def validate_input_data(self, params: Dict[str, Any], required_fields: List[str]) -> bool:
        """
        Validate input data using KBUtilLib patterns.
        
        Args:
            params: Input parameters to validate
            required_fields: List of required field names
            
        Returns:
            True if validation passes
        """
        try:
            for field in required_fields:
                if field not in params or params[field] is None:
                    raise ValueError(f"Required parameter '{field}' is missing or None")
            
            return True
            
        except Exception as e:
            logger.error(f"Input validation failed: {e}")
            raise
    
    def log_info(self, message: str):
        """Log info message using KBUtilLib patterns."""
        if self.base_utils:
            self.base_utils.log_info(message)
        else:
            logger.info(message)
    
    def log_error(self, message: str):
        """Log error message using KBUtilLib patterns."""
        if self.base_utils:
            self.base_utils.log_error(message)
        else:
            logger.error(message)
    
    def log_warning(self, message: str):
        """Log warning message using KBUtilLib patterns.""" 
        if self.base_utils:
            self.base_utils.log_warning(message)
        else:
            logger.warning(message)
