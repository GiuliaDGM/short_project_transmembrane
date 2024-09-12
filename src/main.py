"""
File: main.py
Author: <Giulia Di Gennaro>
Date: <September 12th 2024>
Description: 
    This script is the main entry point for the membrane region detection tool. It executes
    the workflow for analyzing the hydrophobicity and solvent accessibility of protein structures.
    The tool identifies the most likely membrane-exposed regions by slicing the protein along
    different axes, calculating hydrophobicity ratios, and generating membrane planes for
    visualization.

    Workflow:
        1. Parse the protein structure from a PDB file using Biopython.
        2. Run DSSP to compute solvent accessibility and retain accessible residues.
        3. Build a dictionary of residue information and calculate the protein's center of mass.
        4. Generate points on a hemisphere for slicing the protein along different axes.
        5. Use parallel processing to analyze hydrophobicity along each axis.
        6. Identify the optimal membrane axis and generate two membrane planes (above and below 
           the center of mass) for visualization.
        7. Write the results to a new PDB file and generate a PyMOL script for 3D visualization.

Dependencies:
    - Python
    - Biopython: For parsing PDB files and running the DSSP algorithm.
    - NumPy: For numerical operations, array handling, and meshgrid creation.
    - PyMOL: For visualization of the protein and membrane planes.
    - Multiprocessing: For parallel processing of protein slices along different axes.
    - Custom modules:
        - vector3d.py: Handles vector operations in 3D space.
        - fibonnacci_hemisphere_points.py: Generates Fibonacci lattice points on a hemisphere.
        - dssp_analysis.py: Runs DSSP and filters solvent-accessible residues.
        - hydrophobicity_analysis.py: Calculates hydrophobicity for protein slices.
        - protein_struct_com.py: Builds a dictionary of accessible residues and calculates center of mass.
        - best_axis_membrane.py: Identifies the optimal membrane axis and generates membrane planes.
        - visualization.py: Handles PDB file writing and PyMOL script generation.

Usage:
    This script can be run as a standalone program. The user should provide a PDB file, and the script
    will output a modified PDB file with membrane planes and a PyMOL script for visualization.

License: MIT License
"""

# Imports
import numpy as np
from datetime import datetime
from multiprocessing import Pool, cpu_count
from functools import partial
from Bio.PDB import PDBParser

from vector3d import  *
from fibonnacci_hemisphere_points import *
from dssp_analysis import *
from hydrophobicity_analysis import *
from protein_struct_com import *
from best_axis_membrane import *
from visualization import *

if __name__ == "__main__":
    startTime = datetime.now()

    pdb_file = "../data/2moz.pdb"
    thickness = 15
    resolution = 1
    nb_points = 200

    # Load protein structure
    pdb_struct = PDBParser(QUIET=True)
    struct = pdb_struct.get_structure(pdb_file[:4].upper(), pdb_file)
    model = struct[0]

    # Run DSSP and keep only accessible residues
    accessible_residues = run_dssp_and_filter_accessible_residues(pdb_file)
    
    # Build the protein dictionary
    prot_dict, center_of_mass = build_protein_structure_dict_and_com(model, accessible_residues)

    # Generate points on the sphere
    sphere_points = generate_fibonacci_hemisphere(nb_points)

    # Parallel processing of lines
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(partial(process_slicing_plane, prot_dict=prot_dict, thickness=thickness, resolution=resolution), sphere_points)

    # Find best results
    best_results = get_best_results(results)
    
    # best_results[0] is the best_line dictionary, and its "axis_average_hydrophobicity" key contains the coordinates
    best_line = best_results[0]
    coordinates = best_line["axis_average_hydrophobicity"][0]  # This is the [x, y, z] list of coordinates
    hydrophobicity_factor = best_line["axis_average_hydrophobicity"][1]

    display_best_result(coordinates, center_of_mass, hydrophobicity_factor, startTime)

    # Generate and write output files
    pts_mb_1, pts_mb_2 = generate_membranes(results, best_results, resolution)
    new_pdb_file = write_pdb(pdb_file, pts_mb_1, pts_mb_2)
    write_pml_script(coordinates, new_pdb_file, center_of_mass)
