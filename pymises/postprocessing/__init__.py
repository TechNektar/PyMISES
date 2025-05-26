"""
PyMISES - Postprocessing Package

This package provides functionality for post-processing, visualization,
and exporting simulation results.
"""

# Import visualization functions
from pymises.postprocessing.visualize import (
    plot_pressure,
    plot_grid,
    plot_mach_contours,
    plot_streamlines,
    plot_boundary_layer,
    plot_geometry_comparison,
    plot_convergence_history,
    create_performance_report
)

# Import performance metrics functions
from pymises.postprocessing.performance import (
    calculate_airfoil_forces,
    calculate_cascade_performance,
    calculate_boundary_layer_parameters,
    calculate_energy_dissipation,
    analyze_wake
)

# Import export functions
from pymises.postprocessing.export import (
    export_solution_to_csv,
    export_solution_to_vtk,
    export_solution_to_tecplot,
    export_pressure_distribution,
    export_forces_to_json,
    export_boundary_layer_data,
    export_wake_data,
    export_grid,
    export_performance_report,
    export_geometry,
    export_convergence_history
)

# Define package exports
__all__ = [
    # Visualization functions
    'plot_pressure',
    'plot_grid',
    'plot_mach_contours',
    'plot_streamlines',
    'plot_boundary_layer',
    'plot_geometry_comparison',
    'plot_convergence_history',
    'create_performance_report',
    
    # Performance metrics functions
    'calculate_airfoil_forces',
    'calculate_cascade_performance',
    'calculate_boundary_layer_parameters',
    'calculate_energy_dissipation',
    'analyze_wake',
    
    # Export functions
    'export_solution_to_csv',
    'export_solution_to_vtk',
    'export_solution_to_tecplot',
    'export_pressure_distribution',
    'export_forces_to_json',
    'export_boundary_layer_data',
    'export_wake_data',
    'export_grid',
    'export_performance_report',
    'export_geometry',
    'export_convergence_history'
]