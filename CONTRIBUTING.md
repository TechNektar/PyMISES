# Contributing to PyMISES

Thank you for your interest in contributing to PyMISES! This document provides guidelines and information for contributors.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Making Changes](#making-changes)
- [Submitting Changes](#submitting-changes)
- [Code Style](#code-style)
- [Testing](#testing)
- [Documentation](#documentation)

## Code of Conduct

This project follows a code of conduct to ensure a welcoming environment for all contributors. Please be respectful and professional in all interactions.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/yourusername/PyMISES.git
   cd PyMISES
   ```
3. **Add the upstream repository**:
   ```bash
   git remote add upstream https://github.com/originalowner/PyMISES.git
   ```

## Development Setup

### Prerequisites

- Python 3.8 or higher
- Git
- Virtual environment tool (venv, conda, etc.)

### Environment Setup

1. **Create a virtual environment**:
   ```bash
   python -m venv venv
   # On Windows:
   venv\Scripts\activate
   # On Unix/macOS:
   source venv/bin/activate
   ```

2. **Install development dependencies**:
   ```bash
   pip install -e ".[dev]"
   ```

3. **Install pre-commit hooks** (optional but recommended):
   ```bash
   pre-commit install
   ```

### Verify Installation

Run the tests to ensure everything is working:
```bash
pytest
```

## Making Changes

### Before You Start

1. **Check existing issues** and pull requests
2. **Create an issue** if one doesn't exist for your change
3. **Discuss major changes** before implementing

### Workflow

1. **Create a new branch** for your feature/fix:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following the guidelines below

3. **Add tests** for new functionality

4. **Update documentation** if needed

5. **Commit your changes**:
   ```bash
   git add .
   git commit -m "Add: brief description of your changes"
   ```

### Commit Messages

Use clear, descriptive commit messages:
- **Add**: for new features
- **Fix**: for bug fixes
- **Update**: for modifications to existing features
- **Remove**: for deleted features
- **Docs**: for documentation changes

Example:
```
Add: boundary layer transition prediction model

- Implement Abu-Ghannam Shaw transition criterion
- Add tests for transition detection
- Update documentation with transition model details
```

## Submitting Changes

1. **Push your branch**:
   ```bash
   git push origin feature/your-feature-name
   ```

2. **Create a Pull Request** on GitHub with:
   - Clear description of changes
   - Reference to related issues
   - Screenshots/plots if applicable

3. **Respond to feedback** and make requested changes

4. **Ensure CI passes** before requesting final review

## Code Style

### Python Style

- Follow **PEP 8** guidelines
- Use **Black** for code formatting (line length: 88)
- Use **type hints** where appropriate
- Write **docstrings** for all public functions/classes

### Code Formatting

Format your code with Black:
```bash
black pymises/
```

Check style with flake8:
```bash
flake8 pymises/
```

### Docstring Style

Use NumPy-style docstrings:

```python
def calculate_lift_coefficient(pressure_dist, chord_length):
    """
    Calculate lift coefficient from pressure distribution.

    Parameters
    ----------
    pressure_dist : np.ndarray
        Pressure coefficient distribution along airfoil surface.
    chord_length : float
        Airfoil chord length in meters.

    Returns
    -------
    float
        Lift coefficient (dimensionless).

    Examples
    --------
    >>> cp = np.array([-1.2, 0.5, 0.8])
    >>> cl = calculate_lift_coefficient(cp, 1.0)
    >>> print(f"Cl = {cl:.4f}")
    """
```

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=pymises

# Run specific test file
pytest tests/test_core.py

# Run tests matching pattern
pytest -k "test_boundary_layer"
```

### Writing Tests

- Write tests for all new functionality
- Use **pytest** framework
- Place tests in `tests/` directory
- Name test files `test_*.py`
- Use descriptive test names

Example test:
```python
import pytest
import numpy as np
from pymises.core import MISESSolver

class TestMISESSolver:
    def test_solver_initialization(self):
        """Test that solver initializes with default parameters."""
        solver = MISESSolver()
        assert solver.mach == 0.0
        assert solver.reynolds > 0
        
    def test_invalid_mach_number(self):
        """Test that invalid Mach number raises error."""
        solver = MISESSolver()
        with pytest.raises(ValueError):
            solver.set_mach(-0.5)
```

### Test Categories

Mark tests appropriately:
- `@pytest.mark.unit`: Unit tests (fast)
- `@pytest.mark.integration`: Integration tests (slower)
- `@pytest.mark.slow`: Computationally expensive tests

## Documentation

### Building Documentation

```bash
cd docs/
make html
```

### Documentation Style

- Use **Sphinx** with **NumPy docstring** style
- Include **examples** in docstrings
- Add **mathematical equations** using LaTeX when needed
- Keep documentation **up-to-date** with code changes

## Types of Contributions

### Bug Reports

When reporting bugs, include:
- Operating system and Python version
- Steps to reproduce
- Expected vs actual behavior
- Minimal code example
- Relevant error messages

### Feature Requests

For new features:
- Describe the use case
- Explain why it would be valuable
- Consider implementation complexity
- Provide examples if possible

### Code Contributions

We welcome:
- Bug fixes
- New features
- Performance improvements
- Documentation improvements
- Test coverage improvements

## Review Process

1. **Automated checks** must pass (CI/CD)
2. **Code review** by maintainers
3. **Testing** on different platforms
4. **Documentation** review if applicable
5. **Final approval** and merge

## Getting Help

- **GitHub Issues**: For bugs and feature requests
- **GitHub Discussions**: For questions and general discussion
- **Email**: For private inquiries

## Recognition

Contributors will be acknowledged in:
- `CONTRIBUTORS.md` file
- Release notes
- Documentation credits

Thank you for contributing to PyMISES!
