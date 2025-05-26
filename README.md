# PyMISES

**Python Multiple-Blade Interacting Streamtube Euler Solver**

A Python implementation of the MISES solver originally developed by Mark Drela at MIT. This package provides tools for aerodynamic analysis and design of airfoils and turbine cascades using a coupled Euler/boundary-layer approach.

## Features

- **Streamline-based Euler solver** with accurate shock capturing
- **Integral boundary layer solver** with transition prediction
- **Viscous-inviscid coupling** for separated flows
- **Inverse design capabilities** for airfoil optimization
- **Comprehensive post-processing tools** for analysis and visualization
- **Python-native implementation** with modern software engineering practices

## Installation

### Prerequisites

- Python 3.8 or higher
- NumPy
- SciPy
- Matplotlib
- Pandas
- h5py

### Install from source

```bash
git clone https://github.com/yourusername/PyMISES.git
cd PyMISES
pip install -e .
```

### Install dependencies

```bash
pip install -r requirements.txt
```

## Quick Start

```python
import pymises
from pymises.core import MISESSolver

# Create a solver instance
solver = MISESSolver()

# Load an airfoil geometry
airfoil = pymises.load_airfoil('naca0012')

# Set flow conditions
solver.set_conditions(
    mach=0.5,
    reynolds=1e6,
    alpha=2.0  # angle of attack in degrees
)

# Run the analysis
results = solver.solve(airfoil)

# Post-process results
results.plot_pressure_distribution()
results.plot_velocity_field()

print(f"Lift coefficient: {results.cl:.4f}")
print(f"Drag coefficient: {results.cd:.4f}")
```

## Documentation

Comprehensive documentation is available in the `docs/` directory. Key documents include:

- **Installation Guide**: Detailed setup instructions
- **User Manual**: Complete usage examples and tutorials
- **API Reference**: Detailed function and class documentation
- **Theory Guide**: Mathematical background and implementation details

## Project Structure

```
PyMISES/
├── pymises/                    # Main package
│   ├── core/                   # Core solver components
│   ├── boundary_conditions/    # Boundary condition implementations
│   ├── numerics/              # Numerical methods
│   ├── physics/               # Physical models
│   ├── postprocessing/        # Analysis and visualization tools
│   ├── ui/                    # User interface components
│   └── utils/                 # Utility functions
├── tests/                     # Test suite
├── examples/                  # Example scripts and tutorials
├── docs/                      # Documentation
└── requirements.txt           # Python dependencies
```

## Testing

Run the test suite using pytest:

```bash
# Run all tests
pytest

# Run tests with coverage
pytest --cov=pymises

# Run specific test module
pytest tests/test_core.py
```

## Examples

The `examples/` directory contains various usage examples:

- **Basic Analysis**: Simple airfoil analysis
- **Design Optimization**: Inverse design examples
- **Cascade Analysis**: Multi-blade configurations
- **Validation Cases**: Comparison with experimental data

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details on:

- Setting up the development environment
- Code style and standards
- Submitting pull requests
- Reporting issues

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Mark Drela** and the MIT team for the original MISES solver
- The Python scientific computing community for excellent tools and libraries
- Contributors and users who help improve this project

## Citation

If you use PyMISES in your research, please cite:

```bibtex
@software{pymises2024,
  title={PyMISES: Python Multiple-Blade Interacting Streamtube Euler Solver},
  author={PyMISES Development Team},
  year={2024},
  url={https://github.com/yourusername/PyMISES}
}
```

## References

1. Drela, M., "MISES Implementation of Modified Abu-Ghannam Shaw Transition Criterion", MIT Aero-Astro, 1998
2. Drela, M., Youngren, H., "A User's Guide to MISES 2.1", MIT Aero & Astro, 1995

## Status

🚧 **Development Status**: Active development - Beta version

- ✅ Core Euler solver implementation
- ✅ Boundary layer coupling
- ✅ Basic post-processing
- 🔄 Advanced design tools (in progress)
- 📋 Comprehensive validation (planned)

## Support

- **Issues**: Report bugs and request features via [GitHub Issues](https://github.com/yourusername/PyMISES/issues)
- **Discussions**: Join the community discussion on [GitHub Discussions](https://github.com/yourusername/PyMISES/discussions)
- **Email**: For private inquiries, contact [your-email@example.com]

---

**Note**: This is an independent Python implementation inspired by the original MISES solver. It is not officially affiliated with MIT or the original authors.
