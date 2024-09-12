from Bio.PDB import PDBParser, DSSP
import os

def run_dssp_and_filter_accessible_residues(pdb_file, accessibility_threshold=0.3):
    """
    Run DSSP on the provided PDB file and filter residues based on their solvent accessibility.
    
    Args:
        pdb_file (str): The path to the PDB file of the protein.
        accessibility_threshold (float): The minimum solvent accessibility percentage to keep a residue (default is 30%).
    
    Returns:
        dict: A dictionary with residue sequence numbers as keys and their solvent accessibility values as values,
              filtered for those with accessibility >= accessibility_threshold.
    """
    
    # Step 1: Parse the PDB structure
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_file)
    except Exception as e:
        raise ValueError(f"Error parsing PDB file '{pdb_file}': {e}")
    
    model = structure[0]  # DSSP works on the first model (index 0)
    
    # Step 2: Run DSSP to compute solvent accessibility and secondary structure
    try:
        dssp = DSSP(model, pdb_file)
    except Exception as e:
        raise RuntimeError(f"Error running DSSP on '{pdb_file}': {e}")
    
    # Step 3: Filter residues based on solvent accessibility
    accessible_residues = {}
    
    # DSSP property dictionary contains tuples with structure:
    # (AA, secondary_structure, relative_ASA, phi, psi, ...)
    dssp_data = dssp.property_dict
    
    for key, values in dssp_data.items():
        # Extract Relative Solvent Accessibility (Relative ASA) from index 3
        accessibility = values[3]
        
        # Ensure accessibility is a valid float
        if isinstance(accessibility, (float, int)):
            if accessibility >= accessibility_threshold:
                # key is a tuple: (chain_id, (' ', resseq, ' '))
                chain_id, res_id = key
                resseq = res_id[1]  # Residue sequence number
                residue_identifier = f"{chain_id}_{resseq}"
                accessible_residues[residue_identifier] = accessibility

        else:
            # Handle cases where accessibility might be 'None' or invalid
            continue
    
    return accessible_residues
