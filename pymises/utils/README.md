# PyMISES Utilities Module

This directory contains utility functions and helper classes used throughout the PyMISES solver. These utilities provide common functionality needed by multiple modules in the codebase.

## Overview

The utils module provides:

1. **Configuration Management**: Tools for handling configuration settings
2. **Logging**: Customized logging facilities for the solver
3. **Validation**: Input validation and verification utilities

## Modules

### `config.py`

This module implements configuration management for PyMISES.

#### Key Functions and Classes

- **`ConfigManager`**: Manages solver configuration settings.
  - Loads configuration from various sources (files, dictionaries, etc.)
  - Provides hierarchical access to settings
  - Validates configuration against schemas
  - Handles default values and overrides

- **`load_config`**: Loads configuration from a file.
  - Supports JSON, YAML, and INI formats
  - Handles relative paths and environment variables
  - Validates against schema if available

- **`save_config`**: Saves configuration to a file.
  - Preserves comments when possible
  - Supports multiple formats
  - Can generate documentation for settings

#### Implementation Details

- Configuration validation using JSON Schema
- Hierarchical access to nested configuration using dot notation
- Environment variable substitution for deployment flexibility
- Cross-validation between related settings
- Type conversion and validation for numeric values

### `logger.py`

This module implements customized logging facilities for PyMISES.

#### Key Functions and Classes

- **`get_logger`**: Creates and returns a logger instance.
  - Configured with appropriate handlers and formatters
  - Supports different verbosity levels
  - Integrates with module hierarchy
  - Handles both console and file output

- **`LogFormatter`**: Custom formatter for log messages.
  - Adds timestamp, level, and module information
  - Color coding for different log levels
  - Supports different formats for console and file output

- **`configure_logging`**: Sets up global logging configuration.
  - Sets log level based on configuration
  - Configures log file location and rotation
  - Sets up exception handling and reporting

#### Implementation Details

- Uses Python's built-in logging module
- Supports hierarchical loggers for component-specific logging
- Thread-safe implementation for parallel processing
- Performance optimizations for high-volume logging
- Log rotation for managing log file sizes
- Special handling for progress indicators and real-time updates

### `validation.py`

This module implements input validation and verification utilities.

#### Key Functions and Classes

- **`validate_array`**: Validates numpy arrays against specified criteria.
  - Shape validation
  - Value range checking
  - NaN and infinity detection
  - Type checking

- **`validate_airfoil_coordinates`**: Validates airfoil geometry.
  - Checks for closed contour
  - Identifies self-intersections
  - Verifies point distribution quality
  - Normalizes coordinates if requested

- **`validate_grid`**: Validates computational grid properties.
  - Checks for negative cell areas
  - Verifies grid quality metrics
  - Examines node distribution
  - Identifies potential issues for solver stability

- **`validate_boundary_layer_inputs`**: Validates boundary layer solver inputs.
  - Checks for valid edge velocity distribution
  - Verifies Reynolds number ranges
  - Validates surface coordinate monotonicity
  - Checks for missing or inconsistent data

- **`validate_jacobian`**: Validates Jacobian matrix against finite differences.
  - Compares analytical and numerical Jacobians
  - Reports maximum discrepancies
  - Identifies problematic entries
  - Useful for debugging derivative calculations

#### Implementation Details

- Comprehensive error messages for debugging
- Validation functions return True/False for conditional execution
- Type hints for better IDE integration
- Optional validation levels for different stages of development
- Integration with logging system for traceability

## Dependencies

- NumPy: For array operations and validation
- Python logging: For logging functionality

## Example Usage

### Configuration Management

```python
from pymises.utils.config import ConfigManager, load_config

# Load configuration from file
config = load_config('pymises_config.json')

# Create a configuration manager
config_manager = ConfigManager(config)

# Access configuration options
far_field_distance = config_manager.get('grid.far_field_distance', default=20.0)
artificial_dissipation = config_manager.get('euler.artificial_dissipation.enabled', default=True)

# Set a configuration option
config_manager.set('boundary_layer.transition.turbulence_level', 0.01)

# Save updated configuration
config_manager.save('updated_config.json')
```

### Logging

```python
from pymises.utils.logger import get_logger, configure_logging

# Configure global logging
configure_logging(level='INFO', log_file='pymises.log')

# Get a module-specific logger
logger = get_logger(__name__)

# Use the logger
logger.info("Starting grid generation")
logger.debug("Grid parameters: %s", grid_params)

try:
    # Some operation
    result = complex_operation()
except Exception as e:
    logger.error("Operation failed: %s", str(e), exc_info=True)
    raise

logger.info("Grid generation completed successfully")
```

### Validation

```python
from pymises.utils.validation import validate_array, validate_grid

# Validate a numpy array
def process_density(density):
    try:
        validate_array(
            density,
            min_val=0.0,  # Density must be positive
            name="density"
        )
    except ValueError as e:
        logger.error(str(e))
        return False
        
    # Process the validated density
    return True

# Validate a computational grid
def check_grid_quality(x_grid, y_grid):
    try:
        validate_grid(x_grid, y_grid)
        return True
    except ValueError as e:
        logger.warning("Grid validation failed: %s", str(e))
        logger.info("Attempting to improve grid quality...")
        return False
```
