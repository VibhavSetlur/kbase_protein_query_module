"""
Mock biokbase.log module for testing purposes.
"""

import logging
import sys
from typing import Optional, Callable

# Log levels
ERR = 3
WARN = 4
INFO = 6
DEBUG = 7

class MockLog:
    """Mock log class that mimics biokbase.log functionality."""
    
    def __init__(self, name: str, **kwargs):
        self.name = name
        self.logger = logging.getLogger(name)
        self.logger.setLevel(logging.DEBUG)
        
        # Set up console handler if not already present
        if not self.logger.handlers:
            handler = logging.StreamHandler(sys.stdout)
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
    
    def get_log_file(self) -> str:
        """Return a mock log file path."""
        return "/tmp/mock_kbase.log"
    
    def set_log_file(self, log_file: str):
        """Mock method to set log file."""
        pass
    
    def set_log_level(self, level: int):
        """Set the log level."""
        if level <= ERR:
            self.logger.setLevel(logging.ERROR)
        elif level <= WARN:
            self.logger.setLevel(logging.WARNING)
        elif level <= INFO:
            self.logger.setLevel(logging.INFO)
        else:
            self.logger.setLevel(logging.DEBUG)
    
    def log_message(self, level: int, message: str, client_ip: str = None, 
                   user_id: str = None, module: str = None, 
                   method: str = None, call_id: str = None):
        """Log a message with context."""
        if level <= ERR:
            self.logger.error(message)
        elif level <= WARN:
            self.logger.warning(message)
        elif level <= INFO:
            self.logger.info(message)
        else:
            self.logger.debug(message)

def log(name: str, ip_address: bool = False, authuser: bool = False, 
        module: bool = False, method: bool = False, call_id: bool = False,
        changecallback: Optional[Callable] = None, config: Optional[str] = None,
        logfile: Optional[str] = None) -> MockLog:
    """Create a mock log instance."""
    return MockLog(name)
