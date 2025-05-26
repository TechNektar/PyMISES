"""
Logging utilities for PyMISES.

This module provides logging configuration and utilities for
consistent logging throughout the PyMISES software.
"""

import os
import sys
import logging
from typing import Optional, Union, Dict, Any, TextIO

# Define log levels with descriptive names
LEVELS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

class LoggerManager:
    """
    Manager for configuring and accessing loggers in PyMISES.
    
    This class provides a centralized way to configure logging for
    all modules in PyMISES, ensuring consistent log formatting and
    behavior.
    """
    
    def __init__(self):
        """Initialize the logger manager."""
        self.root_logger = logging.getLogger('pymises')
        self.root_logger.setLevel(logging.INFO)
        self.handlers = {}
        
        # By default, don't propagate to the root logger
        self.root_logger.propagate = False
        
        # Create a console handler by default
        self.add_console_handler()
    
    def configure(self, 
                  level: Union[str, int] = 'info',
                  log_file: Optional[str] = None,
                  console: bool = True,
                  format_string: Optional[str] = None) -> None:
        """
        Configure the logger with specified settings.
        
        Args:
            level: Log level ('debug', 'info', 'warning', 'error', 'critical') or logging constant.
            log_file: Path to log file. If provided, logs will be written to this file.
            console: Whether to log to console.
            format_string: Custom format string for log messages. If None, use default format.
        """
        # Set log level
        if isinstance(level, str):
            level = LEVELS.get(level.lower(), logging.INFO)
        self.root_logger.setLevel(level)
        
        # Clear existing handlers
        for handler in list(self.root_logger.handlers):
            self.root_logger.removeHandler(handler)
        self.handlers = {}
        
        # Create a console handler if requested
        if console:
            self.add_console_handler(level, format_string)
        
        # Create a file handler if log_file is provided
        if log_file:
            self.add_file_handler(log_file, level, format_string)
    
    def add_console_handler(self, 
                           level: Union[str, int] = None,
                           format_string: Optional[str] = None) -> None:
        """
        Add a console handler to the logger.
        
        Args:
            level: Log level for this handler. If None, use the root logger's level.
            format_string: Custom format string for log messages. If None, use default format.
        """
        if 'console' in self.handlers:
            return
        
        handler = logging.StreamHandler(sys.stdout)
        
        if level is not None:
            if isinstance(level, str):
                level = LEVELS.get(level.lower(), logging.INFO)
            handler.setLevel(level)
        
        formatter = logging.Formatter(
            format_string or '[%(levelname)s] %(asctime)s - %(name)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)
        
        self.root_logger.addHandler(handler)
        self.handlers['console'] = handler
    
    def add_file_handler(self, 
                        log_file: str,
                        level: Union[str, int] = None,
                        format_string: Optional[str] = None) -> None:
        """
        Add a file handler to the logger.
        
        Args:
            log_file: Path to log file.
            level: Log level for this handler. If None, use the root logger's level.
            format_string: Custom format string for log messages. If None, use default format.
        """
        if 'file' in self.handlers:
            self.root_logger.removeHandler(self.handlers['file'])
        
        # Ensure the directory exists
        os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)
        
        handler = logging.FileHandler(log_file)
        
        if level is not None:
            if isinstance(level, str):
                level = LEVELS.get(level.lower(), logging.INFO)
            handler.setLevel(level)
        
        formatter = logging.Formatter(
            format_string or '[%(levelname)s] %(asctime)s - %(name)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)
        
        self.root_logger.addHandler(handler)
        self.handlers['file'] = handler
    
    def get_logger(self, name: str) -> logging.Logger:
        """
        Get a logger for a specific module.
        
        Args:
            name: Name of the module (typically __name__).
            
        Returns:
            A configured logger for the specified module.
        """
        if name.startswith('pymises.'):
            return logging.getLogger(name)
        else:
            return logging.getLogger(f'pymises.{name}')

# Create a global logger manager instance
logger_manager = LoggerManager()

def setup_logging(level: Union[str, int] = 'info',
                 log_file: Optional[str] = None,
                 console: bool = True,
                 format_string: Optional[str] = None) -> None:
    """
    Configure logging for PyMISES.
    
    Args:
        level: Log level ('debug', 'info', 'warning', 'error', 'critical') or logging constant.
        log_file: Path to log file. If provided, logs will be written to this file.
        console: Whether to log to console.
        format_string: Custom format string for log messages. If None, use default format.
    """
    logger_manager.configure(level, log_file, console, format_string)

def get_logger(name: str) -> logging.Logger:
    """
    Get a logger for a specific module.
    
    Args:
        name: Name of the module (typically __name__).
        
    Returns:
        A configured logger for the specified module.
    """
    return logger_manager.get_logger(name)
