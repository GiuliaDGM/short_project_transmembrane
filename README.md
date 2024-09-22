# SHORT PROJECT - ASSIGNATION AND DETECTION OF PROTEIN TRANSMEMBRANE REGIONS

This project focuses on detecting and assigning the transmembrane regions of a protein using a computational tool. The tool takes a protein's PDB structure and analyzes it to identify membrane-exposed regions based on hydrophobicity and solvent accessibility. The results include membrane planes and a PyMOL visualization script.

## Features
- **Transmembrane region detection**: Detects membrane boundaries using hydrophobicity analysis of protein slices.
- **Multiprocessing**: Efficient parallel computation of different slicing planes.
- **Visualization**: Outputs PyMOL scripts to visualize membrane planes along with the protein structure.

## Download the Repository
```bash
git clone https://github.com/GiuliaDGM/short_project_transmembrane.com
```

## Setup Conda Environment

1. Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
2. Create the Conda environment from the `environment.yml` file:
   ```bash
   conda env create -f environment.yml
   ```
3. Activate the environment:
   ```bash
   conda activate projet_court_transmembranaire
   ```

## Running the Project
The main entry point for the project is the `main.py` file. This script takes a PDB file as input and generates a new PDB file with membrane planes and a PyMOL visualization script.

1. Run the main script with the following command:
   ```bash
   cd src
   python main.py
   ```
2. The output includes:
   - A modified PDB file with the calculated membrane planes.
   - A `.pml` file that can be opened in PyMOL for 3D visualization.

## How to Open the `.pml` File with PyMOL

1. **Install PyMOL**: If you don't have PyMOL installed, you can download it from [here](https://pymol.org/2/).
2. **Open the `.pml` file**: After generating the PyMOL script (`visualization.pml`), follow these steps to visualize it:
   - Open PyMOL.
   - In the PyMOL console, use the command:
     ```bash
     load ../results/visualization.pml
     ```
3. The visualization will display:
   - The protein structure.
   - The optimal membrane planes based on hydrophobicity. (Doesn't do this yet.)
   - The best axis identified for membrane placement.

## Fibonacci Sphere Visualization

To better understand how the protein is analyzed using the Fibonacci sampling method, you can visualize the Fibonacci sphere or hemisphere points in combination with the protein structure.

1. Open the `fibonacci_sphere_visualization.py` file.
2. You can change the PDB file path in this script to visualize a different protein:
   ```python
   pdb_file = '../data/2moz.pdb'  # Change this line to your desired PDB file
   ```
3. Run the script with the following command:
   ```bash
   python fibonacci_sphere_visualization.py
   ```
4. The script will generate:
   - A 3D plot of the Fibonacci sphere or hemisphere aligned with the protein's center of mass.
   - Visualizations are saved in the `results` folder as PNG images.

## Dependencies

- **Python**: The main language for the tool.
- **Biopython**: For parsing PDB files and running DSSP.
- **NumPy**: For numerical operations and array handling..
- **Multiprocessing**: For efficient parallel processing.

All dependencies are managed via Conda, and the full environment specification can be found in the `environment.yml` file.

