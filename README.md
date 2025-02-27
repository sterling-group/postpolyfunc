# Polymer Upcycling

This project processes and functionalizes a polymer supercell using the Atomic Simulation Environment (ASE) library. The script reads an initial structure, creates a supercell, removes edge atoms and their attached hydrogen atoms, identifies extreme carbon atoms, adds hydrogen atoms to these carbons, assigns residues, and functionalizes the carbon atoms.

## Requirements

- Python 3.6+
- ASE (Atomic Simulation Environment)
- NumPy

You can install the required packages using pip:

```bash
pip install ase numpy

Usage
The script can be run from the command line with the following arguments:

Arguments
input_file: Path to the input file containing the initial structure.
output_file: Path to the output file to save the modified structure.
--transformation_matrix: Transformation matrix for the supercell (default: [3, 0, 0, 0, 3, 0, 0, 0, 10]).
--tolerance: Tolerance value for edge detection (default: 1.0).
--functionalization_ratio: Ratio of carbon atoms to be functionalized (default: 0.1).