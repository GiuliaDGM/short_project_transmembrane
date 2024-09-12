from datetime import datetime
import os
import numpy as np
from vector3d import * 

def write_pdb(pdb_file, points_membrane_1, points_membrane_2):
    """
    Copy the original PDB file into a new file and write the coordinates of DUM
    atoms to visualize the membranes.

    Args:
        pdb_file (str): Path to the PDB file.
        points_membrane_1 (numpy.ndarray): 3D coordinates of the points representing the first membrane.
        points_membrane_2 (numpy.ndarray): 3D coordinates of the points representing the second membrane.

    Returns:
        str: Path to the new PDB file with DUM atoms simulating the membranes.
    """
    # Ensure the results directory exists
    output_dir = "../results/"
    os.makedirs(output_dir, exist_ok=True)

    # Define the new PDB file name in the ../results directory
    base_name = os.path.basename(pdb_file)  # Get only the filename from the pdb_file path
    base, ext = os.path.splitext(base_name)
    new_pdb_file = os.path.join(output_dir, base + "_new" + ext)

    # Copy the original PDB file
    with open(pdb_file, "r") as original_file, open(new_pdb_file, "w") as new_file:
        for line in original_file:
            new_file.write(line)

    # Helper function to append points as atoms in the new file
    def write_points(points, start_serial, chain_id, start_resseq):
        """Helper function to write membrane points to the new PDB file."""
        with open(new_pdb_file, "a") as file_in:  # Open in append mode
            for i, point in enumerate(points, start=start_serial):
                atom_type = "DUM"  # Use a consistent atom name for visualization
                residue_seq_number = start_resseq + i
                line = "{:<6}{:>5} {:^4}{:1} {:>3} {:>4}    {:>8.3f} {:>8.3f} {:>8.3f}\n".format(
                    "HETATM", i, atom_type, " ", "MEM", residue_seq_number, float(point[0]), float(point[1]), float(point[2]))
                file_in.write(line)

    # Assign starting serial numbers and residue sequence numbers for membranes
    # Ensure serial numbers do not conflict with original PDB atoms
    write_points(points_membrane_1, start_serial=10000, chain_id='D', start_resseq=10000)
    write_points(points_membrane_2, start_serial=20000, chain_id='D', start_resseq=20000)

    print(f"New PDB file saved at: {new_pdb_file}")
    return new_pdb_file


def write_pml_script(coordinates, pdb_file, center_of_mass):
    """
    Write a PyMol script to visualize the best axis along with the protein structure.
    
    Args:
        sphere_point (tuple or list): A point on the hemisphere representing the best direction.
        pdb_file (str): Path to the PDB file of the protein structure.
    
    Returns:
        str: The full path of the created PyMol script.
    """

    if not isinstance(coordinates, Vector3D):
        coordinates = np.array(coordinates)
        if coordinates.shape != (3,):
            raise ValueError(f"Coordinates should be a 1D array of shape (3,), but got shape {coordinates.shape}")
        coordinates = Vector3D(*coordinates)  # Convert to Vector3D if it's not already
    if not isinstance(center_of_mass, Vector3D):
        center_of_mass = np.array(center_of_mass)
        if center_of_mass.shape != (3,):
            raise ValueError(f"Center of Mass should be a 1D array of shape (3,), but got shape {center_of_mass.shape}")
        center_of_mass = Vector3D(*center_of_mass)  # Convert to Vector3D if it's not already

    point_on_sphere = coordinates + center_of_mass

    # Define the output directory and filename
    output_dir = "../results/"
    os.makedirs(output_dir, exist_ok=True)  # Ensure the results directory exists
    
    # Create the output filename
    pml_filename = os.path.join(output_dir, "visualization.pml")

    # Write the PyMol script to the file
    with open(pml_filename, "w") as file_in:
        file_in.write(f"# PyMol script to visualize the best axis\n")
        file_in.write(f"cmd.load('{pdb_file}')\n\n")
        file_in.write(f"# Create two pseudoatoms representing the best axis\n")
        file_in.write(f"cmd.pseudoatom('pt1', pos=[{-point_on_sphere.x}, {-point_on_sphere.y}, {-point_on_sphere.z}])\n")
        file_in.write(f"cmd.pseudoatom('pt2', pos=[{point_on_sphere.x}, {point_on_sphere.y}, {point_on_sphere.z}])\n\n")
        file_in.write(f"# Draw a dashed line between the two points\n")
        file_in.write(f"cmd.distance('axis', 'pt1', 'pt2')\n\n")
        file_in.write(f"# Customize the appearance of the line\n")
        file_in.write(f"cmd.set('dash_gap', '0')\n")
        file_in.write(f"cmd.set('dash_radius', '0.3')\n")
        file_in.write(f"cmd.set('dash_round_ends', '0')\n")
        file_in.write(f"cmd.set('dash_color', 'yellow', 'axis')\n")
        file_in.write(f"cmd.hide('labels', 'axis')\n\n")
        file_in.write(f"# Show spheres at both points representing the line endpoints\n")
        file_in.write(f"cmd.show('spheres', 'pt1')\n")
        file_in.write(f"cmd.show('spheres', 'pt2')\n\n")
        file_in.write(f"# Set the sphere size for both points\n")
        file_in.write(f"cmd.set('sphere_scale', '0.1', 'pt1')\n")
        file_in.write(f"cmd.set('sphere_scale', '0.1', 'pt2')\n")

    # Print and return the path of the generated PyMol script
    print(f"PyMol script saved at: {pml_filename}\n")
    return pml_filename


def display_best_result(coordinates, center_of_mass, hydrophobicity_factor, startTime):
    """
    Function to display the best axis, hydrophobicity factor, and runtime in a clean format.
    
    Args:
        coordinates (numpy.ndarray or Vector3D): The best direction coordinates (point on the sphere).
        center_of_mass (numpy.ndarray or Vector3D): The center of mass of the protein structure.
        hydrophobicity_factor (float): The calculated hydrophobicity factor for the best direction.
        startTime (datetime): The time when the program started (to calculate runtime).
    
    Outputs:
        A clean, formatted output displaying the best direction, center of mass, hydrophobicity factor, 
        and the program's runtime in a readable format.
    """
    # Ensure that coordinates and center_of_mass are properly formatted
    if not isinstance(coordinates, Vector3D):
        coordinates = np.array(coordinates)
        if coordinates.shape != (3,):
            raise ValueError(f"Coordinates should be a 1D array of shape (3,), but got shape {coordinates.shape}")
        coordinates = Vector3D(*coordinates)  # Convert to Vector3D if it's not already
    if not isinstance(center_of_mass, Vector3D):
        center_of_mass = np.array(center_of_mass)
        if center_of_mass.shape != (3,):
            raise ValueError(f"Center of Mass should be a 1D array of shape (3,), but got shape {center_of_mass.shape}")
        center_of_mass = Vector3D(*center_of_mass)  # Convert to Vector3D if it's not already
    
    # Calculate the point of sphere (shift the best direction coordinates by the center of mass)
    point_of_sphere = coordinates + center_of_mass

    # Output formatting
    print("\n" + "="*50)
    print("{:^50}".format("Best Axis Information"))
    print("="*50)
    print(f"{'Best Direction (Point of Sphere)':<30}: {point_of_sphere}")
    print(f"{'Center of Mass':<30}: {center_of_mass}")
    print(f"{'Highest Hydrophobicity Factor':<30}: {hydrophobicity_factor:.4f}")
    print(f"{'Program Runtime':<30}: {datetime.now() - startTime}")
    print("="*50 + "\n")
