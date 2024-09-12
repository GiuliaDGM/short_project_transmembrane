import numpy as np
import math

from vector3d import *

def calculate_slice_hydrophobicity(residues_in_slice, total_residues_in_slice):
    """
    Calculate the relative hydrophobicity of a slice of residues.
    
    Args:
        residues_in_slice (list): List of residue names in the current slice.
        total_residues_in_slice (int): The total number of residues in the current slice.

    Returns:
        float: The relative hydrophobicity of the slice.
    """
    hydrophobic_residues = ["GLY", "PHE", "ILE", "LEU", "MET", "TRP", "TYR", "VAL"]  # Removed GLY and added LEU

    if total_residues_in_slice == 0:
        # Optionally, return 0.0 or raise an exception
        return 0.0  # Alternatively: raise ValueError("No residues in slice.")

    hydrophobic_count = sum(1 for res in residues_in_slice if res in hydrophobic_residues)
    hydrophobicity_ratio = hydrophobic_count / total_residues_in_slice

    assert isinstance(hydrophobicity_ratio, (float, int)), "Error: rel_hydro should be an int or a float."
    
    return hydrophobicity_ratio

def process_slicing_plane(point_on_hemisphere, prot_dict, thickness, resolution, center_of_mass=(0.0, 0.0, 0.0)):
    """
    Process a single slicing plane through the protein structure and calculate hydrophobicity.

    Args:
        point_on_hemisphere (tuple or list): Coordinates of a point on the hemisphere (x, y, z).
        prot_dict (dict): Dictionary of protein data with residue identifiers as keys and 
                          dictionaries containing 'Coordinates', 'Accessibility', and 'Residue_name' as values.
        thickness (float): Thickness of the slices in angstroms.
        resolution (float): Step size of the sliding plane in angstroms.
        center_of_mass (tuple or list, optional): Coordinates of the protein's center of mass. Defaults to (0.0, 0.0, 0.0).

    Returns:
        dict: A dictionary with results for this slicing plane, including:
              - "slice_hydrophobicity": NumPy array of hydrophobicity ratios per slice.
              - "axis_average_hydrophobicity": Tuple containing the plane normal (Vector3D) and the average hydrophobicity.
              - "total_slices": Total number of slices processed.
              - "min_distance": Minimum distance of residues from the plane.
    """
    # Validate inputs
    if not isinstance(point_on_hemisphere, (list, tuple)) or len(point_on_hemisphere) != 3:
        raise ValueError("point_on_hemisphere must be a list or tuple of three numerical values.")
    if not isinstance(prot_dict, dict):
        raise TypeError("prot_dict must be a dictionary.")
    if not isinstance(thickness, (int, float)) or thickness <= 0:
        raise ValueError("thickness must be a positive number.")
    if not isinstance(resolution, (int, float)) or resolution <= 0:
        raise ValueError("resolution must be a positive number.")
    
    # Normalize the plane normal vector
    try:
        plane_normal = Vector3D(*point_on_hemisphere).normalized()
    except ValueError as e:
        raise

    # Calculate plane offset (d) to pass through the center of mass
    d = - (plane_normal.x * center_of_mass[0] + plane_normal.y * center_of_mass[1] + plane_normal.z * center_of_mass[2])

    # Calculate distances from C-alpha atoms to the plane
    distance_ca_to_plane = {}
    for res_id, infos in prot_dict.items():
        try:
            distance = infos['Coordinates'].dist_to_plane(plane_normal, d)
            distance_ca_to_plane[res_id] = distance
        except Exception as e:
            continue

    if not distance_ca_to_plane:
        return {
            "slice_hydrophobicity": np.array([]),
            "axis_average_hydrophobicity": (plane_normal, 0.0),
            "total_slices": 0,
            "min_distance": 0.0
        }

    # Find the minimum and maximum distances to determine slicing bounds
    min_distance = min(distance_ca_to_plane.values())
    max_distance = max(distance_ca_to_plane.values())

    # Adjust distances by subtracting the shortest distance to start slicing from zero
    adjusted_distances = {res_id: dist - min_distance for res_id, dist in distance_ca_to_plane.items()}

    # Determine the number of steps (slices)
    total_slices = math.ceil((max_distance - min_distance) / resolution)

    # Prepare the result dictionary
    slicing_results = {
        "slice_hydrophobicity": np.zeros(total_slices),  # Array to store hydrophobicity of each slice
        "axis_average_hydrophobicity": None,
        "total_slices": total_slices,
        "min_distance": min_distance
    }

    # Iterate over each slice
    for slice_index in range(total_slices):
        # Define current slice boundaries
        lower_bound = slice_index * resolution
        upper_bound = lower_bound + thickness

        # Collect residues within this slice
        residues_in_slice = [
            prot_dict[res_id]["Residue_name"]
            for res_id, dist in adjusted_distances.items()
            if lower_bound <= dist < upper_bound
        ]

        # Calculate and store hydrophobicity for the slice
        hydrophobicity = calculate_slice_hydrophobicity(residues_in_slice, len(residues_in_slice))
        slicing_results["slice_hydrophobicity"][slice_index] = hydrophobicity

    # Compute the average hydrophobicity for the entire slicing plane
    if total_slices > 0:
        line_average_hydro = np.mean(slicing_results["slice_hydrophobicity"])
    else:
        line_average_hydro = 0.0

    slicing_results["axis_average_hydrophobicity"] = (plane_normal, line_average_hydro)

    return slicing_results