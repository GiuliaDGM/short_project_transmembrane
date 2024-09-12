from vector3d import *

def build_protein_structure_dict_and_com(model, accessible_residues):
    """
    Build a dictionary with protein residue information and calculate the center of mass using a PDB model from Biopython.
    
    Args:
        model (Bio.PDB.Model.Model): The PDB model object obtained from the parsed structure.
        accessible_residues (dict): Dictionary with residue identifiers as keys (e.g., "A_100") and 
                                    their solvent accessibility values as values.
    
    Returns:
        dict: A dictionary with residue identifiers as keys and their associated 
              information (coordinates, accessibility, residue name).
        tuple: The center of mass coordinates (x_com, y_com, z_com) based on all Cα atoms.
    """
    
    protein_structure_dict = {}  # Dictionary to hold residue information
    ca_atom_coords = []  # List to store Cα atom coordinates for center of mass calculation

    # Iterate over all chains and residues in the model
    for chain in model:
        chain_id = chain.id  # Get the chain identifier
        for residue in chain:
            # Skip hetero residues and water molecules
            if not residue.id[0] == ' ':
                continue

            # Check if residue has a Cα atom
            if residue.has_id('CA'):
                ca_atom = residue['CA']
                x, y, z = ca_atom.coord  # Get the coordinates of the Cα atom
                
                # Create a unique residue identifier (e.g., "A_100")
                resseq = residue.get_id()[1]  # Residue sequence number
                residue_identifier = f"{chain_id}_{resseq}"
                
                # If the residue is in the accessible_residues, add to protein_structure_dict
                if residue_identifier in accessible_residues:
                    protein_structure_dict[residue_identifier] = {
                        'Coordinates': Vector3D(x, y, z),
                        'Accessibility': accessible_residues[residue_identifier],
                        'Residue_name': residue.get_resname()
                    }
                
                # Append Cα coordinates for center of mass calculation
                ca_atom_coords.append([x, y, z])

    # Convert the list of coordinates to a NumPy array for vectorized calculation
    if ca_atom_coords:
        coords_array = np.array(ca_atom_coords)
        center_of_mass = tuple(np.mean(coords_array, axis=0))  # Calculate the arithmetic mean
    else:
        center_of_mass = (0.0, 0.0, 0.0)  # Default center of mass if no Cα atoms are found
    
    return protein_structure_dict, center_of_mass
