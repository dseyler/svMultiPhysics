#!/usr/bin/env python3
"""
Analyze center of mass displacement for block_eb test case.

Compares simulations with and without energy balance to show the effect
on rigid body translation.
"""

import vtk
import numpy as np
import os
import glob
import matplotlib.pyplot as plt


def read_vtk_unstructured_grid(filename):
    """Read a VTU file and return the reader output."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def compute_center_of_mass(data):
    """Compute center of mass (centroid) of the mesh points."""
    n_points = data.GetNumberOfPoints()
    com = np.zeros(3)
    for i in range(n_points):
        pt = data.GetPoint(i)
        com += np.array(pt)
    return com / n_points


def compute_displaced_com(data):
    """Compute center of mass in deformed configuration."""
    n_points = data.GetNumberOfPoints()
    disp_array = data.GetPointData().GetArray('Displacement')
    
    if disp_array is None:
        raise ValueError("Displacement array not found in VTK file")
    
    com = np.zeros(3)
    for i in range(n_points):
        pt = np.array(data.GetPoint(i))
        disp = np.array(disp_array.GetTuple3(i))
        com += pt + disp
    return com / n_points


def compute_com_displacement(data):
    """Compute the displacement of the center of mass."""
    n_points = data.GetNumberOfPoints()
    disp_array = data.GetPointData().GetArray('Displacement')
    
    if disp_array is None:
        raise ValueError("Displacement array not found in VTK file")
    
    # Average displacement over all points (= COM displacement for uniform density)
    avg_disp = np.zeros(3)
    for i in range(n_points):
        disp = np.array(disp_array.GetTuple3(i))
        avg_disp += disp
    return avg_disp / n_points


def analyze_results(results_dir, label):
    """Analyze all result files in a directory."""
    result_files = sorted(glob.glob(os.path.join(results_dir, 'result_*.vtu')))
    
    if not result_files:
        print(f"No result files found in {results_dir}")
        return None
    
    print(f"\n{'='*60}")
    print(f"Analysis: {label}")
    print(f"{'='*60}")
    print(f"Found {len(result_files)} result files")
    
    # Get initial mesh (reference configuration)
    data0 = read_vtk_unstructured_grid(result_files[0])
    com0 = compute_center_of_mass(data0)
    print(f"\nInitial center of mass: ({com0[0]:.6e}, {com0[1]:.6e}, {com0[2]:.6e})")
    
    timesteps = []
    com_displacements = []
    com_disp_magnitudes = []
    
    for i, f in enumerate(result_files):
        data = read_vtk_unstructured_grid(f)
        com_disp = compute_com_displacement(data)
        com_disp_mag = np.linalg.norm(com_disp)
        
        timesteps.append(i + 1)
        com_displacements.append(com_disp)
        com_disp_magnitudes.append(com_disp_mag)
        
        print(f"  Step {i+1:3d}: COM displacement = ({com_disp[0]:+.6e}, {com_disp[1]:+.6e}, {com_disp[2]:+.6e}), |d| = {com_disp_mag:.6e}")
    
    # Final summary
    final_disp = com_displacements[-1]
    final_mag = com_disp_magnitudes[-1]
    print(f"\nFinal COM displacement: ({final_disp[0]:+.6e}, {final_disp[1]:+.6e}, {final_disp[2]:+.6e})")
    print(f"Final COM displacement magnitude: {final_mag:.6e}")
    
    return {
        'timesteps': timesteps,
        'com_displacements': np.array(com_displacements),
        'com_disp_magnitudes': np.array(com_disp_magnitudes),
        'label': label
    }


def plot_comparison(results_eb, results_no_eb, output_file='com_displacement_comparison.png'):
    """Plot comparison of COM displacement with and without energy balance."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot X, Y, Z components
    components = ['X', 'Y', 'Z']
    for i, comp in enumerate(components):
        ax = axes[i // 2, i % 2]
        if results_eb is not None:
            ax.plot(results_eb['timesteps'], results_eb['com_displacements'][:, i], 
                   'b-o', label='With Energy Balance', markersize=4)
        if results_no_eb is not None:
            ax.plot(results_no_eb['timesteps'], results_no_eb['com_displacements'][:, i], 
                   'r-s', label='Without Energy Balance', markersize=4)
        ax.set_xlabel('Time Step')
        ax.set_ylabel(f'{comp} Displacement')
        ax.set_title(f'COM {comp}-Displacement')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Plot magnitude
    ax = axes[1, 1]
    if results_eb is not None:
        ax.plot(results_eb['timesteps'], results_eb['com_disp_magnitudes'], 
               'b-o', label='With Energy Balance', markersize=4)
    if results_no_eb is not None:
        ax.plot(results_no_eb['timesteps'], results_no_eb['com_disp_magnitudes'], 
               'r-s', label='Without Energy Balance', markersize=4)
    ax.set_xlabel('Time Step')
    ax.set_ylabel('Displacement Magnitude')
    ax.set_title('COM Displacement Magnitude')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"\nPlot saved to: {output_file}")
    plt.close()


def main():
    """Main analysis function."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    results_eb_dir = os.path.join(script_dir, 'results')
    results_no_eb_dir = os.path.join(script_dir, 'results_no_eb')
    
    results_eb = None
    results_no_eb = None
    
    # Analyze with energy balance
    if os.path.exists(results_eb_dir):
        results_eb = analyze_results(results_eb_dir, "WITH Energy Balance")
    else:
        print(f"\nResults directory not found: {results_eb_dir}")
        print("Run: mpirun -np 1 svmultiphysics solver.xml")
    
    # Analyze without energy balance
    if os.path.exists(results_no_eb_dir):
        results_no_eb = analyze_results(results_no_eb_dir, "WITHOUT Energy Balance")
    else:
        print(f"\nResults directory not found: {results_no_eb_dir}")
        print("Run: mpirun -np 1 svmultiphysics solver_no_eb.xml")
    
    # Compare results
    if results_eb is not None and results_no_eb is not None:
        print(f"\n{'='*60}")
        print("COMPARISON SUMMARY")
        print(f"{'='*60}")
        
        final_eb = results_eb['com_displacements'][-1]
        final_no_eb = results_no_eb['com_displacements'][-1]
        
        print(f"\nFinal COM displacement WITH energy balance:    ({final_eb[0]:+.6e}, {final_eb[1]:+.6e}, {final_eb[2]:+.6e})")
        print(f"Final COM displacement WITHOUT energy balance: ({final_no_eb[0]:+.6e}, {final_no_eb[1]:+.6e}, {final_no_eb[2]:+.6e})")
        print(f"\nDifference: ({final_eb[0]-final_no_eb[0]:+.6e}, {final_eb[1]-final_no_eb[1]:+.6e}, {final_eb[2]-final_no_eb[2]:+.6e})")
        
        ratio_x = abs(final_eb[0] / final_no_eb[0]) if abs(final_no_eb[0]) > 1e-15 else float('inf')
        print(f"\nX-displacement ratio (EB/No-EB): {ratio_x:.4f}")
        
        # Create comparison plot
        plot_comparison(results_eb, results_no_eb, os.path.join(script_dir, 'com_displacement_comparison.png'))
    elif results_eb is not None or results_no_eb is not None:
        # Create plot with available data
        plot_comparison(results_eb, results_no_eb, os.path.join(script_dir, 'com_displacement_comparison.png'))


if __name__ == '__main__':
    main()
