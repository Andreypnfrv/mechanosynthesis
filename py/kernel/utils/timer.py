"""
Timer utility functions for mechanosynthesis research pipeline.
"""

import time
from typing import Tuple


def format_elapsed_time(elapsed_seconds: float) -> str:
    """
    Format elapsed time in seconds to a human-readable Russian format.
    
    Args:
        elapsed_seconds (float): Time elapsed in seconds
        
    Returns:
        str: Formatted time string (e.g., "1ч 23м 45с", "15м 30с", "45с")
    """
    hours, remainder = divmod(int(elapsed_seconds), 3600)
    minutes, seconds = divmod(remainder, 60)
    
    if hours > 0:
        return f"{hours}ч {minutes}м {seconds}с"
    elif minutes > 0:
        return f"{minutes}м {seconds}с"
    else:
        return f"{seconds}с"


def calculate_and_format_time(start_time: float) -> Tuple[float, str]:
    """
    Calculate elapsed time from start_time and return both elapsed time and formatted string.
    
    Args:
        start_time (float): Start time from time.time()
        
    Returns:
        tuple: (elapsed_seconds, formatted_time_string)
    """
    end_time = time.time()
    elapsed_time = end_time - start_time
    formatted_time = format_elapsed_time(elapsed_time)
    return elapsed_time, formatted_time


class Timer:
    """
    Context manager for timing operations with automatic formatting.
    
    Usage:
        with Timer() as timer:
            # do some work
            pass
        print(f"Operation completed {timer.formatted_time}")
        
    Or manually:
        timer = Timer()
        timer.start()
        # do some work
        timer.stop()
        print(f"Operation took {timer.formatted_time}")
    """
    
    def __init__(self):
        self.start_time = None
        self.end_time = None
        self.elapsed_time = None
        self.formatted_time = None
    
    def start(self):
        """Start the timer."""
        self.start_time = time.time()
        return self
    
    def stop(self):
        """Stop the timer and calculate elapsed time."""
        if self.start_time is None:
            raise RuntimeError("Timer not started")
        
        self.end_time = time.time()
        self.elapsed_time = self.end_time - self.start_time
        self.formatted_time = format_elapsed_time(self.elapsed_time)
        return self.elapsed_time
    
    def __enter__(self):
        """Context manager entry."""
        self.start()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.stop()
        return False  # Don't suppress exceptions


def print_completion_message(operation: str, filename: str, elapsed_time: float):
    """
    Print a standardized completion message with timing.
    
    Args:
        operation (str): Name of the operation (e.g., "Orca optimization", "Hessian calculation")
        filename (str): Name of the file being processed
        elapsed_time (float): Time elapsed in seconds
    """
    formatted_time = format_elapsed_time(elapsed_time)
    print(f"{operation} completed for {filename} за {formatted_time}")


def print_failure_message(operation: str, filename: str, elapsed_time: float, error_msg: str = ""):
    """
    Print a standardized failure message with timing.
    
    Args:
        operation (str): Name of the operation (e.g., "Orca optimization", "Hessian calculation")
        filename (str): Name of the file being processed
        elapsed_time (float): Time elapsed in seconds
        error_msg (str): Optional error message
    """
    formatted_time = format_elapsed_time(elapsed_time)
    if error_msg:
        print(f"{operation} failed for {filename} после {formatted_time}: {error_msg}")
    else:
        print(f"{operation} failed for {filename} после {formatted_time}")


def print_timeout_message(operation: str, filename: str, elapsed_time: float):
    """
    Print a standardized timeout message with timing.
    
    Args:
        operation (str): Name of the operation (e.g., "Orca optimization", "Hessian calculation")
        filename (str): Name of the file being processed
        elapsed_time (float): Time elapsed in seconds
    """
    formatted_time = format_elapsed_time(elapsed_time)
    print(f"{operation} timed out for {filename} после {formatted_time}")