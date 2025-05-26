# PyMISES Test Suite

This directory contains tests for the PyMISES (Multiple-blade Interacting Streamtube Euler Solver) implementation.

## Test Organization

The tests are organized into the following directories:

- `core/`: Tests for core functionality (grid, euler, boundary_layer, coupling, etc.)
- `boundary_conditions/`: Tests for boundary condition implementations
- `numerics/`: Tests for numerical methods (jacobian, linear solver, etc.)
- `physics/`: Tests for physical models (transition, turbulence, etc.)
- `integration/`: End-to-end tests that validate the complete workflow

## Running Tests

You can run the tests using the `run_tests.py` script in the root directory:

```bash
# Run all tests
python run_tests.py

# Run specific module tests
python run_tests.py grid euler

# Run integration tests only
python run_tests.py -i

# Run unit tests only
python run_tests.py -u

# Run tests with verbose output
python run_tests.py -v

# List available test modules
python run_tests.py -l
```

Alternatively, you can use pytest directly:

```bash
# Run all tests
pytest

# Run specific test module
pytest tests/core/test_grid.py

# Run specific test class
pytest tests/core/test_grid.py::TestStreamlineGrid

# Run specific test
pytest tests/core/test_grid.py::TestStreamlineGrid::test_initialization
```

## Writing New Tests

When adding new functionality, please also add corresponding tests. Follow these guidelines:

1. Place unit tests in the appropriate directory based on the module being tested
2. Name test files with the pattern `test_*.py`
3. Name test classes with the pattern `Test*`
4. Name test methods with the pattern `test_*`
5. Use descriptive names that indicate what is being tested
6. Include docstrings explaining the purpose of the test
7. Keep tests focused on testing a single aspect of functionality
8. Use appropriate assertions to validate results

For integration tests:
1. Place integration tests in the `integration/` directory
2. Name files based on the workflow being tested (e.g., `test_airfoil_analysis.py`)
3. Minimize dependency on external files or resources
4. Use simplified test cases that run relatively quickly

## Test Data

The `data/` directory contains test data used by the tests, such as:
- Reference airfoil geometries
- Validation data from experiments
- Pre-computed reference solutions

Please ensure that any test data added to the repository is small and well-documented.
