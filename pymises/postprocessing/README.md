# PyMISES Postprocessing Module

This directory contains tools for postprocessing and visualizing PyMISES solutions. These modules provide functionality for extracting meaningful information from numerical results and presenting it in useful formats.

## Overview

The postprocessing module provides:

1. **Visualization Tools**: Functions for plotting results and flow fields
2. **Performance Analysis**: Tools for calculating and interpreting aerodynamic performance
3. **Export Utilities**: Functions for exporting results to various formats

## Modules

### `visualize.py`

This module implements visualization tools for PyMISES solutions.

#### Key Functions

- **`plot_pressure`**: Creates pressure distribution plots.
  - Surface pressure coefficient (Cp) vs. x/c
  - Options for upper/lower surface identification
  - Comparison with target or experimental data
  - Shock and separation indicators

- **`plot_grid`**: Visualizes computational grids.
  - Full domain or zoomed views
  - Grid quality visualization
  - Streamline overlay options
  - Boundary condition identification

- **`plot_mach_contours`**: Creates Mach number contour plots.
  - Contour levels with customizable ranges
  - Shock detection and highlighting
  - Streamline overlays
  - Boundary layer edge indicators

- **`plot_boundary_layer`**: Visualizes boundary layer properties.
  - Displacement thickness and shape factor
  - Skin friction distribution
  - Transition and separation locations
  - Comparison with experimental data

- **`plot_convergence_history`**: Displays convergence history.
  - Residual reduction over iterations
  - Force coefficient convergence
  - Multiple convergence metrics
  - Semi-log and linear scaling options

#### Implementation Details

- Modern Matplotlib-based plotting with customizable styling
- Support for interactive visualization when used in notebooks
- SVG export for publication-quality figures
- Multiple data series for comparison studies
- Intelligent default parameters for common use cases
- Consistent color schemes and styling

### `performance.py`

This module implements tools for calculating and analyzing aerodynamic performance.

#### Key Functions

- **`calculate_forces`**: Computes aerodynamic forces and moments.
  - Lift, drag, and moment coefficients
  - Pressure and viscous contributions
  - Integration along specified surfaces
  - Moment reference point customization

- **`analyze_boundary_layer`**: Analyzes boundary layer properties.
  - Integrated boundary layer parameters
  - Transition and separation detection
  - Skin friction and dissipation analysis
  - Loss generation mechanisms

- **`cascade_performance`**: Calculates performance for turbomachinery cascades.
  - Flow turning and deflection
  - Total pressure loss coefficient
  - Diffusion factor
  - Zweifel loading coefficient

- **`off_design_analysis`**: Evaluates performance at off-design conditions.
  - Angle of attack or incidence sweeps
  - Mach number sensitivity
  - Drag divergence detection
  - Stall margin estimation

#### Implementation Details

- Force computation through surface pressure integration
- Viscous drag calculation from momentum thickness at trailing edge
- Wake analysis for profile loss quantification
- Vectorized implementation for efficiency with large datasets
- Robust handling of irregular geometries and grids
- Accurate treatment of leading and trailing edge regions

### `export.py`

This module implements functions for exporting PyMISES results to various formats.

#### Key Functions

- **`export_to_vtk`**: Exports solution to VTK format.
  - Compatible with ParaView and other visualization tools
  - Full flow field and boundary layer properties
  - Structured and unstructured grid support
  - Multi-block format for complex geometries

- **`export_to_tecplot`**: Exports solution to Tecplot format.
  - Zone-based organization for different regions
  - ASCII and binary format options
  - Support for multiple solution variables
  - Auxiliary data for boundary conditions

- **`export_pressure_distribution`**: Exports pressure distribution data.
  - CSV format for easy plotting in other tools
  - Upper and lower surface identification
  - Includes x/c, Cp, Mach, and boundary layer data
  - Optional interpolation to uniform spacing

- **`export_performance_report`**: Generates a comprehensive performance report.
  - Text-based summary of key performance metrics
  - HTML report with embedded figures
  - LaTeX report template for documentation
  - Comparison with reference data

- **`export_geometry`**: Exports geometric information.
  - Airfoil coordinates in Selig or Lednicer format
  - Structured grid in PLOT3D format
  - STL export for 3D visualization
  - CAD-compatible formats

#### Implementation Details

- Consistent file naming conventions
- Options for compression of large datasets
- Metadata inclusion for provenance tracking
- Batch processing capabilities for multiple cases
- Built-in validation of exported data
- Support for both absolute and normalized coordinates

## Dependencies

- NumPy: For array operations
- Matplotlib: For visualization
- VTK: For VTK file export (optional)
- Pandas: For data manipulation and CSV export

## Example Usage

### Plotting Pressure Distribution

```python
from pymises.postprocessing.visualize import plot_pressure

# After running a simulation and obtaining the solution
solution = euler_solver.get_solution()
airfoil = airfoil_geometry

# Create pressure distribution plot
fig, ax = plot_pressure(
    solution, 
    airfoil,
    show_upper_lower=True,
    include_mach=True,
    show_separation=True,
    comparison_data=experimental_data
)

# Add additional information
ax.set_title(f'NACA 0012, M={mach}, α={alpha}°, Re={reynolds:.1e}')

# Save the figure
fig.savefig('pressure_distribution.png', dpi=300)
```

### Calculating Performance Metrics

```python
from pymises.postprocessing.performance import calculate_forces, analyze_boundary_layer

# After running a simulation
solution = coupled_solver.get_solution()
airfoil = airfoil_geometry

# Calculate forces
forces = calculate_forces(
    solution, 
    airfoil,
    reference_chord=1.0,
    reference_area=1.0,
    moment_reference=[0.25, 0.0]  # Quarter-chord
)

# Analyze boundary layer
bl_analysis = analyze_boundary_layer(
    solution,
    airfoil,
    detect_transition=True,
    detect_separation=True
)

# Print results
print(f"Lift coefficient: {forces['cl']:.4f}")
print(f"Drag coefficient: {forces['cd']:.6f}")
print(f"Moment coefficient: {forces['cm']:.4f}")

if 'transition_location' in bl_analysis:
    print(f"Transition (x/c): Upper={bl_analysis['transition_location']['upper']:.3f}, " 
          f"Lower={bl_analysis['transition_location']['lower']:.3f}")

if 'separation_location' in bl_analysis:
    if bl_analysis['separation_location']['upper'] < 1.0:
        print(f"Separation on upper surface at x/c={bl_analysis['separation_location']['upper']:.3f}")
```

### Exporting Results

```python
from pymises.postprocessing.export import export_to_vtk, export_pressure_distribution, export_performance_report

# After solving and obtaining the solution
solution = coupled_solver.get_solution()
airfoil = airfoil_geometry
grid = computational_grid

# Export to VTK for visualization in ParaView
export_to_vtk(
    solution,
    grid,
    filename='airfoil_flow.vtk',
    include_boundary_layer=True,
    include_vorticity=True
)

# Export pressure distribution for plotting in Excel/etc.
export_pressure_distribution(
    solution,
    airfoil,
    filename='pressure_distribution.csv',
    variables=['x/c', 'cp', 'mach', 'cf', 'h']
)

# Generate a comprehensive performance report
export_performance_report(
    solution,
    filename='performance_report.html',
    format='html',
    include_figures=True,
    simulation_parameters={
        'mach': mach,
        'alpha': alpha,
        'reynolds': reynolds,
        'airfoil': 'NACA 0012'
    }
)
```
