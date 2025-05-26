"""
PyMISES - Export Module

This module provides functions for exporting simulation results to various formats,
including CSV, VTK, and Tecplot.
"""

import numpy as np
import os
import csv
import json
from typing import Dict, List, Tuple, Optional, Union

def export_solution_to_csv(solution: Dict, filename: str):
    """
    Export solution data to CSV file.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
    filename : str
        Output filename
    """
    # Check if required data is present
    required_fields = ['grid_x', 'grid_y', 'density', 'velocity_x', 'velocity_y', 'pressure']
    for field in required_fields:
        if field not in solution:
            raise ValueError(f"Solution missing required field: {field}")
    
    # Extract data
    x = solution['grid_x']
    y = solution['grid_y']
    rho = solution['density']
    vx = solution['velocity_x']
    vy = solution['velocity_y']
    p = solution['pressure']
    
    # Combine into a single array
    data = np.column_stack((x, y, rho, vx, vy, p))
    
    # Additional fields if available
    if 'mach' in solution:
        mach = solution['mach']
        data = np.column_stack((data, mach))
        header = ['x', 'y', 'density', 'velocity_x', 'velocity_y', 'pressure', 'mach']
    else:
        header = ['x', 'y', 'density', 'velocity_x', 'velocity_y', 'pressure']
    
    # Add boundary layer data if available
    if 'displacement_thickness' in solution:
        delta_star = solution['displacement_thickness']
        data = np.column_stack((data, delta_star))
        header.append('displacement_thickness')
    
    if 'momentum_thickness' in solution:
        theta = solution['momentum_thickness']
        data = np.column_stack((data, theta))
        header.append('momentum_thickness')
    
    # Write to CSV
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerows(data)

def export_solution_to_vtk(solution: Dict, filename: str):
    """
    Export solution data to VTK file for visualization in ParaView or VisIt.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
    filename : str
        Output filename
    """
    # Check if required data is present
    required_fields = ['grid_x', 'grid_y', 'density', 'velocity_x', 'velocity_y', 'pressure']
    for field in required_fields:
        if field not in solution:
            raise ValueError(f"Solution missing required field: {field}")
    
    # Extract data
    x = solution['grid_x']
    y = solution['grid_y']
    rho = solution['density']
    vx = solution['velocity_x']
    vy = solution['velocity_y']
    p = solution['pressure']
    
    # Get grid dimensions
    if 'grid_ni' in solution and 'grid_nj' in solution:
        ni = solution['grid_ni']
        nj = solution['grid_nj']
    else:
        # Try to infer grid dimensions
        # This assumes the grid is structured and ordered correctly
        ni = 0
        nj = 0
        x_diff = np.diff(x)
        y_diff = np.diff(y)
        
        # Find the first large jump in coordinates
        for i in range(1, len(x)):
            if abs(x_diff[i-1]) > 100 * abs(x_diff[i-2]) or abs(y_diff[i-1]) > 100 * abs(y_diff[i-2]):
                nj = i
                break
        
        if nj > 0:
            ni = len(x) // nj
            if ni * nj != len(x):
                raise ValueError("Cannot infer grid dimensions from data")
        else:
            raise ValueError("Cannot infer grid dimensions from data")
    
    # Reshape data to grid dimensions
    x_grid = x.reshape(ni, nj)
    y_grid = y.reshape(ni, nj)
    rho_grid = rho.reshape(ni, nj)
    vx_grid = vx.reshape(ni, nj)
    vy_grid = vy.reshape(ni, nj)
    p_grid = p.reshape(ni, nj)
    
    # Add z-coordinate for 2D grid (all zeros)
    z_grid = np.zeros_like(x_grid)
    
    # Open VTK file for writing
    with open(filename, 'w') as f:
        # Write VTK header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("PyMISES Solution\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {nj} {ni} 1\n")
        f.write(f"POINTS {ni*nj} float\n")
        
        # Write grid points
        for i in range(ni):
            for j in range(nj):
                f.write(f"{x_grid[i,j]} {y_grid[i,j]} {z_grid[i,j]}\n")
        
        # Write scalar and vector data
        f.write(f"POINT_DATA {ni*nj}\n")
        
        # Density
        f.write("SCALARS density float\n")
        f.write("LOOKUP_TABLE default\n")
        for i in range(ni):
            for j in range(nj):
                f.write(f"{rho_grid[i,j]}\n")
        
        # Pressure
        f.write("SCALARS pressure float\n")
        f.write("LOOKUP_TABLE default\n")
        for i in range(ni):
            for j in range(nj):
                f.write(f"{p_grid[i,j]}\n")
        
        # Velocity vector
        f.write("VECTORS velocity float\n")
        for i in range(ni):
            for j in range(nj):
                f.write(f"{vx_grid[i,j]} {vy_grid[i,j]} 0.0\n")
        
        # Mach number if available
        if 'mach' in solution:
            mach = solution['mach']
            mach_grid = mach.reshape(ni, nj)
            
            f.write("SCALARS mach float\n")
            f.write("LOOKUP_TABLE default\n")
            for i in range(ni):
                for j in range(nj):
                    f.write(f"{mach_grid[i,j]}\n")

def export_solution_to_tecplot(solution: Dict, filename: str):
    """
    Export solution data to Tecplot format.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
    filename : str
        Output filename
    """
    # Check if required data is present
    required_fields = ['grid_x', 'grid_y', 'density', 'velocity_x', 'velocity_y', 'pressure']
    for field in required_fields:
        if field not in solution:
            raise ValueError(f"Solution missing required field: {field}")
    
    # Extract data
    x = solution['grid_x']
    y = solution['grid_y']
    rho = solution['density']
    vx = solution['velocity_x']
    vy = solution['velocity_y']
    p = solution['pressure']
    
    # Get grid dimensions
    if 'grid_ni' in solution and 'grid_nj' in solution:
        ni = solution['grid_ni']
        nj = solution['grid_nj']
    else:
        # Try to infer grid dimensions (same as in VTK export)
        ni = 0
        nj = 0
        x_diff = np.diff(x)
        y_diff = np.diff(y)
        
        # Find the first large jump in coordinates
        for i in range(1, len(x)):
            if abs(x_diff[i-1]) > 100 * abs(x_diff[i-2]) or abs(y_diff[i-1]) > 100 * abs(y_diff[i-2]):
                nj = i
                break
        
        if nj > 0:
            ni = len(x) // nj
            if ni * nj != len(x):
                raise ValueError("Cannot infer grid dimensions from data")
        else:
            raise ValueError("Cannot infer grid dimensions from data")
    
    # Reshape data to grid dimensions
    x_grid = x.reshape(ni, nj)
    y_grid = y.reshape(ni, nj)
    rho_grid = rho.reshape(ni, nj)
    vx_grid = vx.reshape(ni, nj)
    vy_grid = vy.reshape(ni, nj)
    p_grid = p.reshape(ni, nj)
    
    # Define variable names
    variables = ['x', 'y', 'density', 'velocity_x', 'velocity_y', 'pressure']
    
    # Add Mach number if available
    if 'mach' in solution:
        mach = solution['mach']
        mach_grid = mach.reshape(ni, nj)
        variables.append('mach')
    
    # Open Tecplot file for writing
    with open(filename, 'w') as f:
        # Write header
        f.write('TITLE = "PyMISES Solution"\n')
        f.write('VARIABLES = ' + ', '.join(f'"{var}"' for var in variables) + '\n')
        f.write(f'ZONE T="Zone 1", I={nj}, J={ni}, F=POINT\n')
        
        # Write data
        for i in range(ni):
            for j in range(nj):
                line = [x_grid[i,j], y_grid[i,j], rho_grid[i,j], vx_grid[i,j], vy_grid[i,j], p_grid[i,j]]
                
                # Add Mach number if available
                if 'mach' in solution:
                    line.append(mach_grid[i,j])
                
                f.write(' '.join([f"{val}" for val in line]) + '\n')

def export_grid(grid, filename: str, format: str = 'plot3d'):
    """
    Export computational grid to file.
    
    Parameters
    ----------
    grid : Grid
        Computational grid object
    filename : str
        Output filename
    format : str, optional
        Grid format ('plot3d', 'cgns', 'vtk'), defaults to 'plot3d'
    """
    # Extract grid coordinates
    x = grid.x
    y = grid.y
    ni = grid.ni
    nj = grid.nj
    
    if format.lower() == 'plot3d':
        # Export to PLOT3D format
        with open(filename, 'w') as f:
            # Write number of blocks
            f.write("1\n")
            
            # Write grid dimensions
            f.write(f"{ni} {nj}\n")
            
            # Write x-coordinates
            for i in range(ni):
                for j in range(nj):
                    f.write(f"{x[i, j]}\n")
            
            # Write y-coordinates
            for i in range(ni):
                for j in range(nj):
                    f.write(f"{y[i, j]}\n")
            
            # Write z-coordinates (all zeros for 2D grid)
            for i in range(ni):
                for j in range(nj):
                    f.write("0.0\n")
    
    elif format.lower() == 'vtk':
        # Export to VTK format
        with open(filename, 'w') as f:
            # Write VTK header
            f.write("# vtk DataFile Version 3.0\n")
            f.write("PyMISES Grid\n")
            f.write("ASCII\n")
            f.write("DATASET STRUCTURED_GRID\n")
            f.write(f"DIMENSIONS {nj} {ni} 1\n")
            f.write(f"POINTS {ni*nj} float\n")
            
            # Write grid points
            for i in range(ni):
                for j in range(nj):
                    f.write(f"{x[i,j]} {y[i,j]} 0.0\n")
    
    else:
        raise ValueError(f"Unsupported grid format: {format}")

def export_performance_report(solution: Dict, filename: str):
    """
    Export performance report to a text file.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing performance data
    filename : str
        Output filename
    """
    # Generate report text
    report = "PyMISES Performance Report\n"
    report += "========================\n\n"
    
    # Flow conditions
    report += "Flow Conditions\n"
    report += "--------------\n"
    if 'mach' in solution:
        report += f"Mach Number: {solution['mach']:.3f}\n"
    if 'reynolds' in solution:
        report += f"Reynolds Number: {solution['reynolds']:.2e}\n"
    if 'alpha' in solution:
        report += f"Angle of Attack: {np.degrees(solution['alpha']):.2f}°\n"
    report += "\n"
    
    # Aerodynamic coefficients
    report += "Aerodynamic Coefficients\n"
    report += "-----------------------\n"
    if 'cl' in solution:
        report += f"Lift Coefficient (CL): {solution['cl']:.4f}\n"
    if 'cd' in solution:
        report += f"Drag Coefficient (CD): {solution['cd']:.6f}\n"
    if 'cm' in solution:
        report += f"Moment Coefficient (CM): {solution['cm']:.4f}\n"
    if 'cl' in solution and 'cd' in solution:
        report += f"Lift-to-Drag Ratio (L/D): {solution['cl']/solution['cd']:.2f}\n"
    report += "\n"
    
    # Boundary layer properties
    report += "Boundary Layer Properties\n"
    report += "------------------------\n"
    if 'transition_location' in solution:
        trans_loc = solution['transition_location']
        if isinstance(trans_loc, dict):
            # For airfoil with upper/lower surfaces
            report += f"Transition Location (Upper): {trans_loc.get('upper', 'N/A'):.4f} x/c\n"
            report += f"Transition Location (Lower): {trans_loc.get('lower', 'N/A'):.4f} x/c\n"
        else:
            # For cascade with single surface
            report += f"Transition Location: {trans_loc:.4f} x/c\n"
    
    if 'separation_location' in solution:
        sep_loc = solution['separation_location']
        if isinstance(sep_loc, dict):
            # For airfoil with upper/lower surfaces
            upper_sep = sep_loc.get('upper', -1)
            lower_sep = sep_loc.get('lower', -1)
            if upper_sep > 0:
                report += f"Separation Location (Upper): {upper_sep:.4f} x/c\n"
            else:
                report += "Separation Location (Upper): None\n"
            if lower_sep > 0:
                report += f"Separation Location (Lower): {lower_sep:.4f} x/c\n"
            else:
                report += "Separation Location (Lower): None\n"
        else:
            # For cascade with single surface
            if sep_loc > 0:
                report += f"Separation Location: {sep_loc:.4f} x/c\n"
            else:
                report += "Separation Location: None\n"
    report += "\n"
    
    # Additional cascade performance metrics
    if 'loss_coefficient' in solution:
        report += "Cascade Performance\n"
        report += "-------------------\n"
        report += f"Total Pressure Loss Coefficient: {solution['loss_coefficient']:.6f}\n"
    if 'diffusion_factor' in solution:
        report += f"Diffusion Factor: {solution['diffusion_factor']:.4f}\n"
    if 'outlet_angle' in solution:
        report += f"Outlet Flow Angle: {np.degrees(solution['outlet_angle']):.2f}°\n"
    
    # Write to file
    with open(filename, 'w') as f:
        f.write(report)

def export_geometry(geometry, filename: str, format: str = 'dat'):
    """
    Export geometry to file.
    
    Parameters
    ----------
    geometry : AirfoilGeometry or BladeGeometry
        Geometry object
    filename : str
        Output filename
    format : str, optional
        Geometry format ('dat', 'igs', 'stp'), defaults to 'dat'
    """
    # Extract coordinates
    x = geometry.x
    y = geometry.y
    
    if format.lower() == 'dat':
        # Export to simple text format (e.g., Selig format for airfoils)
        with open(filename, 'w') as f:
            # Write header with geometry name
            name = os.path.splitext(os.path.basename(filename))[0]
            f.write(f"{name}\n")
            
            # Write coordinates
            for i in range(len(x)):
                f.write(f"{x[i]:.6f} {y[i]:.6f}\n")
    
    elif format.lower() in ['igs', 'iges']:
        # This would require a proper IGES writer
        # For simplicity, we'll just export as a simple text file
        with open(filename, 'w') as f:
            f.write("IGES export not implemented yet.\n")
            f.write("Coordinates in simple text format:\n")
            for i in range(len(x)):
                f.write(f"{x[i]:.6f} {y[i]:.6f}\n")
    
    elif format.lower() in ['stp', 'step']:
        # This would require a proper STEP writer
        # For simplicity, we'll just export as a simple text file
        with open(filename, 'w') as f:
            f.write("STEP export not implemented yet.\n")
            f.write("Coordinates in simple text format:\n")
            for i in range(len(x)):
                f.write(f"{x[i]:.6f} {y[i]:.6f}\n")
    
    else:
        raise ValueError(f"Unsupported geometry format: {format}")
    
def export_convergence_history(convergence_history: dict, filename: str):
    """
    Export convergence history to CSV file.
    
    Parameters
    ----------
    convergence_history : dict
        Dictionary containing convergence history data
    filename : str
        Output filename
    """
    # Check if convergence history is a dict with iteration data
    if isinstance(convergence_history, dict):
        # Get iterations
        iterations = convergence_history.get('iterations', np.arange(1, len(next(iter(convergence_history.values()))) + 1))
        
        # Get residual types
        residual_types = [key for key in convergence_history.keys() if key != 'iterations']
        
        # Prepare data for export
        data = np.column_stack([iterations] + [convergence_history[key] for key in residual_types])
        
        # Write to CSV
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['iteration'] + residual_types)
            writer.writerows(data)
    
    # Check if convergence history is a simple array or list
    elif isinstance(convergence_history, (np.ndarray, list)):
        # Prepare data with iteration numbers
        iterations = np.arange(1, len(convergence_history) + 1)
        data = np.column_stack((iterations, convergence_history))
        
        # Write to CSV
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['iteration', 'residual'])
            writer.writerows(data)
    
    else:
        raise ValueError("Convergence history must be a dictionary or array-like object")

def export_solution_to_hdf5(solution: Dict, filename: str):
    """
    Export solution data to HDF5 file for efficient storage and retrieval.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
    filename : str
        Output filename
    """
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py package is required for HDF5 export.")
    
    # Check if required data is present
    required_fields = ['grid_x', 'grid_y', 'density', 'velocity_x', 'velocity_y', 'pressure']
    for field in required_fields:
        if field not in solution:
            raise ValueError(f"Solution missing required field: {field}")
    
    # Create HDF5 file
    with h5py.File(filename, 'w') as f:
        # Create groups
        grid_group = f.create_group('grid')
        flow_group = f.create_group('flow')
        meta_group = f.create_group('metadata')
        
        # Add grid data
        grid_group.create_dataset('x', data=solution['grid_x'])
        grid_group.create_dataset('y', data=solution['grid_y'])
        
        if 'grid_ni' in solution and 'grid_nj' in solution:
            grid_group.attrs['ni'] = solution['grid_ni']
            grid_group.attrs['nj'] = solution['grid_nj']
        
        # Add flow data
        flow_group.create_dataset('density', data=solution['density'])
        flow_group.create_dataset('velocity_x', data=solution['velocity_x'])
        flow_group.create_dataset('velocity_y', data=solution['velocity_y'])
        flow_group.create_dataset('pressure', data=solution['pressure'])
        
        # Add additional flow variables if available
        if 'mach' in solution:
            flow_group.create_dataset('mach', data=solution['mach'])
        
        if 'temperature' in solution:
            flow_group.create_dataset('temperature', data=solution['temperature'])
        
        # Add boundary layer data if available
        if 'displacement_thickness' in solution:
            bl_group = f.create_group('boundary_layer')
            bl_group.create_dataset('displacement_thickness', data=solution['displacement_thickness'])
            
            if 'momentum_thickness' in solution:
                bl_group.create_dataset('momentum_thickness', data=solution['momentum_thickness'])
            
            if 'cf' in solution:
                bl_group.create_dataset('skin_friction', data=solution['cf'])
            
            if 'transition_location' in solution:
                trans_loc = solution['transition_location']
                if isinstance(trans_loc, dict):
                    for key, value in trans_loc.items():
                        bl_group.attrs[f'transition_{key}'] = value
                else:
                    bl_group.attrs['transition'] = trans_loc
        
        # Add performance coefficients if available
        perf_group = f.create_group('performance')
        for coef in ['cl', 'cd', 'cm']:
            if coef in solution:
                perf_group.attrs[coef] = solution[coef]
        
        # Add flow conditions as metadata
        for param in ['mach', 'reynolds', 'alpha', 'p_inf', 'T_inf']:
            if param in solution:
                meta_group.attrs[param] = solution[param]
        
        # Add convergence history if available
        if 'convergence_history' in solution:
            conv_history = solution['convergence_history']
            if isinstance(conv_history, dict):
                conv_group = f.create_group('convergence')
                for key, values in conv_history.items():
                    conv_group.create_dataset(key, data=values)
            elif isinstance(conv_history, (np.ndarray, list)):
                f.create_dataset('convergence/residual', data=conv_history)

def export_mesh_to_gmsh(grid, filename: str):
    """
    Export grid to Gmsh format.
    
    Parameters
    ----------
    grid : Grid
        Computational grid object
    filename : str
        Output filename
    """
    # Extract grid coordinates
    x = grid.x
    y = grid.y
    ni = grid.ni
    nj = grid.nj
    
    # Open file for writing
    with open(filename, 'w') as f:
        # Write Gmsh header
        f.write('$MeshFormat\n')
        f.write('2.2 0 8\n')
        f.write('$EndMeshFormat\n')
        
        # Write nodes
        f.write('$Nodes\n')
        f.write(f'{ni*nj}\n')
        
        node_id = 1
        for i in range(ni):
            for j in range(nj):
                f.write(f'{node_id} {x[i,j]} {y[i,j]} 0.0\n')
                node_id += 1
        
        f.write('$EndNodes\n')
        
        # Write elements (quads)
        f.write('$Elements\n')
        
        # Calculate number of elements
        n_elements = (ni-1) * (nj-1)
        f.write(f'{n_elements}\n')
        
        # Element type 3 is quad
        element_id = 1
        for i in range(ni-1):
            for j in range(nj-1):
                # Calculate node IDs for this quad
                node1 = i*nj + j + 1
                node2 = i*nj + (j+1) + 1
                node3 = (i+1)*nj + (j+1) + 1
                node4 = (i+1)*nj + j + 1
                
                # Write element
                f.write(f'{element_id} 3 2 0 0 {node1} {node2} {node3} {node4}\n')
                element_id += 1
        
        f.write('$EndElements\n')

def export_solution_to_json(solution: Dict, filename: str):
    """
    Export solution data to JSON format.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
    filename : str
        Output filename
    """
    import json
    
    # Create a serializable version of the solution
    json_data = {}
    
    # Add scalar values
    for key, value in solution.items():
        if isinstance(value, (int, float, str, bool)) or value is None:
            json_data[key] = value
        elif isinstance(value, dict):
            # Handle nested dictionaries
            json_data[key] = {}
            for sub_key, sub_value in value.items():
                if isinstance(sub_value, (int, float, str, bool)) or sub_value is None:
                    json_data[key][sub_key] = sub_value
                elif isinstance(sub_value, np.ndarray):
                    json_data[key][sub_key] = sub_value.tolist()
        elif isinstance(value, np.ndarray):
            # Convert numpy arrays to lists
            json_data[key] = value.tolist()
    
    # Write to file
    with open(filename, 'w') as f:
        json.dump(json_data, f, indent=4)

def export_surface_geometry(geometry, filename: str, format: str = 'stl'):
    """
    Export surface geometry to file formats suitable for CAD or meshing.
    
    Parameters
    ----------
    geometry : AirfoilGeometry or BladeGeometry
        Geometry object
    filename : str
        Output filename
    format : str, optional
        Geometry format ('stl', 'obj', 'dxf'), defaults to 'stl'
    """
    # Extract coordinates
    x = geometry.x
    y = geometry.y
    
    if format.lower() == 'stl':
        # Create a simple STL representing the 2D geometry as a thin 3D object
        with open(filename, 'w') as f:
            f.write("solid PyMISES_Geometry\n")
            
            # Create triangles with a small z thickness
            z_thickness = 0.01
            
            # Front face triangles
            for i in range(len(x) - 2):
                f.write("  facet normal 0 0 1\n")
                f.write("    outer loop\n")
                f.write(f"      vertex {x[0]} {y[0]} {z_thickness/2}\n")
                f.write(f"      vertex {x[i+1]} {y[i+1]} {z_thickness/2}\n")
                f.write(f"      vertex {x[i+2]} {y[i+2]} {z_thickness/2}\n")
                f.write("    endloop\n")
                f.write("  endfacet\n")
            
            # Back face triangles
            for i in range(len(x) - 2):
                f.write("  facet normal 0 0 -1\n")
                f.write("    outer loop\n")
                f.write(f"      vertex {x[0]} {y[0]} {-z_thickness/2}\n")
                f.write(f"      vertex {x[i+2]} {y[i+2]} {-z_thickness/2}\n")
                f.write(f"      vertex {x[i+1]} {y[i+1]} {-z_thickness/2}\n")
                f.write("    endloop\n")
                f.write("  endfacet\n")
            
            # Side faces (connecting front and back)
            for i in range(len(x) - 1):
                f.write("  facet normal 0 0 0\n")  # Proper normal calculation omitted for simplicity
                f.write("    outer loop\n")
                f.write(f"      vertex {x[i]} {y[i]} {z_thickness/2}\n")
                f.write(f"      vertex {x[i+1]} {y[i+1]} {z_thickness/2}\n")
                f.write(f"      vertex {x[i+1]} {y[i+1]} {-z_thickness/2}\n")
                f.write("    endloop\n")
                f.write("  endfacet\n")
                
                f.write("  facet normal 0 0 0\n")
                f.write("    outer loop\n")
                f.write(f"      vertex {x[i]} {y[i]} {z_thickness/2}\n")
                f.write(f"      vertex {x[i+1]} {y[i+1]} {-z_thickness/2}\n")
                f.write(f"      vertex {x[i]} {y[i]} {-z_thickness/2}\n")
                f.write("    endloop\n")
                f.write("  endfacet\n")
            
            f.write("endsolid PyMISES_Geometry\n")
    
    elif format.lower() == 'obj':
        # Export as Wavefront OBJ format
        with open(filename, 'w') as f:
            f.write("# PyMISES Geometry\n")
            
            # Write vertices
            for i in range(len(x)):
                f.write(f"v {x[i]} {y[i]} 0.0\n")
            
            # Write another set of vertices with small z offset for thickness
            z_thickness = 0.01
            for i in range(len(x)):
                f.write(f"v {x[i]} {y[i]} {z_thickness}\n")
            
            # Write faces (quads connecting front and back)
            for i in range(len(x) - 1):
                v1 = i + 1  # OBJ indices start at 1
                v2 = i + 2
                v3 = i + 2 + len(x)
                v4 = i + 1 + len(x)
                f.write(f"f {v1} {v2} {v3} {v4}\n")
    
    elif format.lower() == 'dxf':
        try:
            import ezdxf
        except ImportError:
            raise ImportError("ezdxf package is required for DXF export.")
        
        # Create DXF document
        doc = ezdxf.new('R2010')
        msp = doc.modelspace()
        
        # Add polyline for the geometry
        points = [(x[i], y[i]) for i in range(len(x))]
        msp.add_lwpolyline(points, close=True)
        
        # Save DXF file
        doc.saveas(filename)
    
    else:
        raise ValueError(f"Unsupported geometry format: {format}")
    
def export_pressure_distribution(solution: Dict, geometry, filename: str):
    """
    Export pressure distribution around an airfoil or blade.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing pressure and grid data
    geometry : AirfoilGeometry or BladeGeometry
        Geometry object containing airfoil or blade information
    filename : str
        Output filename
    """
    # Extract data
    pressure = solution['pressure']
    grid_x = solution['grid_x']
    grid_y = solution['grid_y']
    
    # Get airfoil surface points
    surface_indices = geometry.get_surface_indices()
    
    # Compute pressure coefficient
    if 'mach' in solution and 'p_inf' in solution:
        mach = solution['mach']
        p_inf = solution['p_inf']
        cp = (pressure[surface_indices] - p_inf) / (0.5 * p_inf * 1.4 * mach**2)
    else:
        # Normalize pressure for visualization if no reference values available
        p_min = np.min(pressure[surface_indices])
        p_max = np.max(pressure[surface_indices])
        cp = (pressure[surface_indices] - p_min) / (p_max - p_min)
    
    # Get surface coordinates
    x = grid_x[surface_indices]
    y = grid_y[surface_indices]
    
    # Calculate arc length along surface
    s = np.zeros_like(x)
    for i in range(1, len(x)):
        s[i] = s[i-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2)
    
    # Normalize arc length
    s = s / s[-1]
    
    # Prepare data for export
    data = np.column_stack((x, y, s, pressure[surface_indices], cp))
    
    # Write to CSV
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['x', 'y', 's', 'pressure', 'Cp'])
        writer.writerows(data)

def export_boundary_layer_data(solution: Dict, geometry, filename: str):
    """
    Export boundary layer data along the surface.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing boundary layer data
    geometry : AirfoilGeometry or BladeGeometry
        Geometry object containing airfoil or blade information
    filename : str
        Output filename
    """
    # Check if boundary layer data is present
    if 'displacement_thickness' not in solution or 'momentum_thickness' not in solution:
        raise ValueError("Solution missing boundary layer data")
    
    # Extract boundary layer data
    delta_star = solution['displacement_thickness']
    theta = solution['momentum_thickness']
    H = delta_star / theta  # Shape factor
    
    # Get surface coordinates
    surface_indices = geometry.get_surface_indices()
    x = solution['grid_x'][surface_indices]
    y = solution['grid_y'][surface_indices]
    
    # Calculate arc length along surface
    s = np.zeros_like(x)
    for i in range(1, len(x)):
        s[i] = s[i-1] + np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2)
    
    # Normalize arc length
    s = s / s[-1]
    
    # Prepare data for export
    data = np.column_stack((x, y, s, delta_star, theta, H))
    
    # Additional boundary layer data if available
    header = ['x', 'y', 's', 'delta_star', 'theta', 'H']
    
    if 'cf' in solution:
        cf = solution['cf']
        data = np.column_stack((data, cf))
        header.append('cf')
    
    if 'transition' in solution:
        transition = solution['transition']
        data = np.column_stack((data, transition))
        header.append('transition')
    
    # Write to CSV
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerows(data)

def export_wake_data(solution: Dict, wake_location: float, filename: str):
    """
    Export wake data at a specific downstream location.
    
    Parameters
    ----------
    solution : Dict
        Solution dictionary containing flow data
    wake_location : float
        Downstream location (x/c) for wake analysis
    filename : str
        Output filename
    """
    # Extract grid information
    if 'grid_x' not in solution or 'grid_y' not in solution:
        raise ValueError("Solution missing grid data")
    
    x = solution['grid_x']
    y = solution['grid_y']
    
    # Find points at the specified downstream location
    tolerance = 0.01  # Tolerance for finding points at wake_location
    wake_indices = np.where(np.abs(x - wake_location) < tolerance)[0]
    
    if len(wake_indices) == 0:
        raise ValueError(f"No points found at downstream location x/c = {wake_location}")
    
    # Extract flow variables at wake
    rho_wake = solution['density'][wake_indices]
    vx_wake = solution['velocity_x'][wake_indices]
    vy_wake = solution['velocity_y'][wake_indices]
    p_wake = solution['pressure'][wake_indices]
    
    # Calculate velocity magnitude at wake
    v_wake = np.sqrt(vx_wake**2 + vy_wake**2)
    
    # Get y-coordinates at wake
    y_wake = y[wake_indices]
    
    # Sort points by y-coordinate
    sort_idx = np.argsort(y_wake)
    y_wake = y_wake[sort_idx]
    v_wake = v_wake[sort_idx]
    vx_wake = vx_wake[sort_idx]
    vy_wake = vy_wake[sort_idx]
    p_wake = p_wake[sort_idx]
    rho_wake = rho_wake[sort_idx]
    
    # Calculate normalized velocity deficit
    v_max = np.max(v_wake)
    v_deficit = 1.0 - v_wake / v_max
    
    # Prepare data for export
    data = np.column_stack((y_wake, rho_wake, vx_wake, vy_wake, v_wake, p_wake, v_deficit))
    
    # Write to CSV
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['y', 'density', 'velocity_x', 'velocity_y', 'velocity_magnitude', 'pressure', 'deficit'])
        writer.writerows(data)
    
