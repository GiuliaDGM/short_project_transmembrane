"""
This code is NOT NECESSARY to the project's input. 
It is merely used to understand and visulise the uniformity of the Fibonacci sphere/hemisphere genrated.
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import PDBParser

# Ensure results directory exists
output_dir = "../results/"
os.makedirs(output_dir, exist_ok=True)

# Function to calculate the center of gravity of the protein
def calculate_center_of_gravity(pdb_file):
    """
    Calculate the center of gravity (centroid) of all atoms in a PDB file.
    
    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        numpy.ndarray: The (x, y, z) coordinates of the center of gravity.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    atoms = [atom for atom in structure.get_atoms() if atom.get_coord() is not None]
    coords = np.array([atom.get_coord() for atom in atoms])
    center_of_gravity = coords.mean(axis=0)
    return center_of_gravity

# Function to generate points on a Fibonacci sphere
def fibonacci_sphere(samples=1000):
    """
    Generate points on a sphere using Fibonacci sampling.

    Args:
        samples (int): Number of points to generate on the sphere.

    Returns:
        numpy.ndarray: Array of shape (samples, 3) containing (x, y, z) coordinates.
    """
    points = []
    phi = math.pi * (math.sqrt(5.) - 1.)  # Golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))

    return np.array(points)

# Function to plot the sphere and center of gravity
def plot_sphere_with_center(center_of_gravity, sphere_points):
    """
    Plot a 3D visualization of Fibonacci sphere points and the protein's center of gravity.

    Args:
        center_of_gravity (numpy.ndarray): The (x, y, z) coordinates of the center of gravity.
        sphere_points (numpy.ndarray): Array of shape (samples, 3) containing the sphere points.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Shift sphere points to match the center of gravity
    shifted_sphere_points = sphere_points + center_of_gravity
    ax.scatter(shifted_sphere_points[:, 0], shifted_sphere_points[:, 1], shifted_sphere_points[:, 2], 
               c=shifted_sphere_points[:, 2], cmap='viridis', alpha=0.6, label='Fibonacci sphere points', s=2)
    
    # Plot the center of gravity
    ax.scatter(center_of_gravity[0], center_of_gravity[1], center_of_gravity[2], color='red', s=50, label='Center of mass')
    
    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Fibonacci Sphere with Center of Mass')
    ax.legend()

    # Save the plot
    output_path = os.path.join(output_dir, "fibonacci_sphere_with_center.png")
    plt.savefig(output_path)
    plt.show()
    print(f"\nPlot saved to {output_path}\n")

# Load C-alpha coordinates from PDB file
def load_pdb_ca_coordinates(pdb_filename):
    """
    Load atomic coordinates of C-alpha atoms from a PDB file.

    Args:
        pdb_filename (str): Path to the PDB file.

    Returns:
        numpy.ndarray: An array of shape (num_ca_atoms, 3) containing the (x, y, z) coordinates of the C-alpha atoms.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_filename)
    
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_coords.append(residue['CA'].get_coord())  # Extract C-alpha coordinates
    
    return np.array(ca_coords)

# Visualize points on a hemisphere and the protein structure
def visualize_fibonacci_hemisphere_with_protein(scale_up, points, protein_coords):
    """
    Visualize Fibonacci hemisphere points and a protein structure, aligning the protein's
    center of gravity with the hemisphere center. Only points in the upper hemisphere are displayed.

    Args:
        points (numpy.ndarray): Array of shape (num_points, 3) containing sphere points.
        protein_coords (numpy.ndarray): Array of shape (num_atoms, 3) containing protein C-alpha coordinates.
    """
    # Filter points to only include the upper hemisphere (z >= 0)
    hemisphere_points = points[points[:, 2] >= 0] * scale_up
    
    # Calculate the center of gravity directly from the protein coordinates
    protein_center_of_gravity = np.mean(protein_coords, axis=0)
    
    # Translate protein coordinates so that its center of gravity matches the sphere's center
    translation_vector = -protein_center_of_gravity
    translated_protein_coords = protein_coords + translation_vector

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot hemisphere points (z >= 0)
    ax.scatter(hemisphere_points[:, 0], hemisphere_points[:, 1], hemisphere_points[:, 2], 
               c=hemisphere_points[:, 2], cmap='viridis', alpha=0.6, label='Fibonacci hemisphere points')
    
    # Plot protein atoms
    ax.scatter(translated_protein_coords[:, 0], translated_protein_coords[:, 1], translated_protein_coords[:, 2], 
               c='red', marker='o', alpha=0.7, label='Protein structure')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Fibonacci Hemisphere with Protein Structure')
    ax.legend()

    output_path = os.path.join(output_dir, "fibonacci_hemisphere_with_protein.png")
    plt.savefig(output_path)
    plt.show()
    print(f"\nPlot saved to {output_path}\n")
    print("\nFor our project, we will pass different axis through different points of the sphere and the center of mass.")
    print("We will then calculate the hydrophobicity ratio along each axis. \nThe axis with the highest hydrophobicity ratio is selected as the 'Best Axis'.\n")

# Main execution
if __name__ == "__main__":
    pdb_file = '../data/2moz.pdb'  # Replace with your actual PDB file path
    n_points = 1000
    scale_up = 100
    
    # Calculate center of gravity and generate Fibonacci sphere points
    center_of_gravity = calculate_center_of_gravity(pdb_file)
    sphere_points = fibonacci_sphere(samples=n_points)
    
    # Plot the Fibonacci sphere with the protein's center of gravity
    plot_sphere_with_center(center_of_gravity, sphere_points)
    
    # Load protein C-alpha coordinates and visualize them with Fibonacci points
    protein_ca_coords = load_pdb_ca_coordinates(pdb_file)
    visualize_fibonacci_hemisphere_with_protein(scale_up, sphere_points, protein_ca_coords)
