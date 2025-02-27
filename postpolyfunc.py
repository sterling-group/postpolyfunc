import ase
from ase.visualize import view
from ase.io import read,write
from ase.neighborlist import neighbor_list

def remove_edge_atoms_and_hydrogens(atoms, tolerance=1.0):
    """
    Removes atoms that are on the edges of the supercell, where their coordinates are within a tolerance of 0 in either x, y, or z,
    and also removes hydrogen atoms that are attached to these edge atoms.
    
    Parameters:
    atoms (ase.Atoms): The supercell as an ASE Atoms object.
    tolerance (float): The tolerance value for edge detection.
    
    Returns:
    ase.Atoms: The modified supercell with edge atoms and their attached hydrogen atoms removed.
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    # Identify edge atoms
    edge_indices = []
    for i, pos in enumerate(positions):
        if (abs(pos[0]) < tolerance or abs(pos[0] - atoms.get_cell()[0, 0]) < tolerance or
            abs(pos[1]) < tolerance or abs(pos[1] - atoms.get_cell()[1, 1]) < tolerance or
            abs(pos[2]) < tolerance or abs(pos[2] - atoms.get_cell()[2, 2]) < tolerance):
            edge_indices.append(i)
    
    # Use neighbor list to find hydrogen atoms attached to edge atoms
    i_list, j_list = neighbor_list('ij', atoms, 1.5, self_interaction=False)
    hydrogen_indices = set()
    for i, j in zip(i_list, j_list):
        if i in edge_indices and symbols[j] == 'H':
            hydrogen_indices.add(j)
        elif j in edge_indices and symbols[i] == 'H':
            hydrogen_indices.add(i)
    
    # Combine edge and hydrogen indices to remove
    indices_to_remove = set(edge_indices).union(hydrogen_indices)
    
    # Filter out the atoms to be removed
    filtered_positions = [pos for i, pos in enumerate(positions) if i not in indices_to_remove]
    filtered_symbols = [sym for i, sym in enumerate(symbols) if i not in indices_to_remove]
    
    # Create a new Atoms object with the filtered positions and symbols
    new_atoms = Atoms(symbols=filtered_symbols, positions=filtered_positions, cell=atoms.get_cell(), pbc=atoms.get_pbc())
    
    return new_atoms

# Define the transformation matrix for the supercell
transformation_matrix = [[3, 0, 0], [0, 3, 0], [0, 0, 10]]

# Create the supercell
supercell = make_supercell(atom, transformation_matrix)

# Remove edge atoms and their attached hydrogen atoms with a tolerance of 1.0
clean_supercell = remove_edge_atoms_and_hydrogens(supercell, tolerance=1.0)

def identify_extreme_carbons(atoms):
    """
    Identifies the carbon atoms with the smallest and largest z-coordinate values.
    
    Parameters:
    atoms (ase.Atoms): The supercell as an ASE Atoms object.
    
    Returns:
    tuple: Lists of indices of the carbon atoms with the smallest and largest z-coordinate values.
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    min_z = float('inf')
    max_z = float('-inf')
    
    # Find the minimum and maximum z-coordinate values
    for pos, sym in zip(positions, symbols):
        if sym == 'C':
            if pos[2] < min_z:
                min_z = pos[2]
            if pos[2] > max_z:
                max_z = pos[2]
    
    # Find all carbon atoms with the minimum and maximum z-coordinate values
    min_z_indices = [i for i, (pos, sym) in enumerate(zip(positions, symbols)) if sym == 'C' and pos[2] == min_z]
    max_z_indices = [i for i, (pos, sym) in enumerate(zip(positions, symbols)) if sym == 'C' and pos[2] == max_z]
    
    return min_z_indices, max_z_indices

def identify_extreme_carbons(atoms):
    """
    Identifies the carbon atoms with the smallest and largest z-coordinate values.
    
    Parameters:
    atoms (ase.Atoms): The supercell as an ASE Atoms object.
    
    Returns:
    tuple: Lists of indices of the carbon atoms with the smallest and largest z-coordinate values.
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    min_z = float('inf')
    max_z = float('-inf')
    
    # Find the minimum and maximum z-coordinate values
    for pos, sym in zip(positions, symbols):
        if sym == 'C':
            if pos[2] < min_z:
                min_z = pos[2]
            if pos[2] > max_z:
                max_z = pos[2]
    
    # Find all carbon atoms with the minimum and maximum z-coordinate values
    min_z_indices = [i for i, (pos, sym) in enumerate(zip(positions, symbols)) if sym == 'C' and pos[2] == min_z]
    max_z_indices = [i for i, (pos, sym) in enumerate(zip(positions, symbols)) if sym == 'C' and pos[2] == max_z]
    
    return min_z_indices, max_z_indices

def add_hydrogens_to_extreme_carbons(atoms, min_z_indices, max_z_indices, bond_length=1.09):
    """
    Adds hydrogen atoms to the carbon atoms with the smallest and largest z-coordinate values to convert them from CH2 to CH3.
    
    Parameters:
    atoms (ase.Atoms): The supercell as an ASE Atoms object.
    min_z_indices (list): Indices of the carbon atoms with the smallest z-coordinate values.
    max_z_indices (list): Indices of the carbon atoms with the largest z-coordinate values.
    bond_length (float): The bond length for the new hydrogen atoms.
    
    Returns:
    ase.Atoms: The modified supercell with added hydrogen atoms.
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    new_positions = positions.tolist()
    new_symbols = symbols[:]
    
    # Use neighbor list to find existing hydrogen atoms attached to the extreme carbon atoms
    i_list, j_list = neighbor_list('ij', atoms, 1.5, self_interaction=False)
    
    def add_hydrogen(carbon_index, direction, mirror=False):
        carbon_pos = positions[carbon_index]
        attached_hydrogens = [positions[j] for i, j in zip(i_list, j_list) if i == carbon_index and symbols[j] == 'H']
        
        if len(attached_hydrogens) == 2:
            # Calculate the position for the new hydrogen atom
            vec1 = attached_hydrogens[0] - carbon_pos
            vec2 = attached_hydrogens[1] - carbon_pos
            if mirror:
                new_h_pos = carbon_pos - direction * bond_length * (vec1 + vec2) / np.linalg.norm(vec1 + vec2)
            else:
                new_h_pos = carbon_pos + direction * bond_length * (vec1 + vec2) / np.linalg.norm(vec1 + vec2)
            new_positions.append(new_h_pos)
            new_symbols.append('H')
    
    # Add hydrogen atoms to the carbon atoms with the smallest z-coordinate values
    for index in min_z_indices:
        add_hydrogen(index, -1)
    
    # Add hydrogen atoms to the carbon atoms with the largest z-coordinate values (mirror image)
    for index in max_z_indices:
        add_hydrogen(index, 1, mirror=True)
    
    # Create a new Atoms object with the added hydrogen atoms
    new_atoms = Atoms(symbols=new_symbols, positions=new_positions, cell=atoms.get_cell(), pbc=atoms.get_pbc())
    
    return new_atoms

def assign_residues(atoms):
    """
    Assigns different residue names or numbers to each strand in the supercell.
    
    Parameters:
    atoms (ase.Atoms): The supercell as an ASE Atoms object.
    
    Returns:
    ase.Atoms: The modified supercell with assigned residue information.
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    residues = []
    
    # Assign residue numbers based on the z-coordinate
    residue_number = 1
    current_z = positions[0][2]
    for pos in positions:
        if abs(pos[2] - current_z) > 1.0:  # Adjust the threshold as needed
            residue_number += 1
            current_z = pos[2]
        residues.append(residue_number)
    
    # Create a new Atoms object with the residue information
    new_atoms = Atoms(symbols=symbols, positions=positions, cell=atoms.get_cell(), pbc=atoms.get_pbc())
    new_atoms.set_array('residues', np.array(residues))
    
    return new_atoms

def functionalize_carbons(atoms, functionalization_ratio=0.1):
    """
    Selects random carbon atoms, deletes one hydrogen and replaces the other with an oxygen.
    
    Parameters:
    atoms (ase.Atoms): The supercell as an ASE Atoms object.
    functionalization_ratio (float): The ratio of carbon atoms to be functionalized.
    
    Returns:
    ase.Atoms: The modified supercell with functionalized carbon atoms.
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    # Identify all carbon atoms
    carbon_indices = [i for i, sym in enumerate(symbols) if sym == 'C']
    
    # Calculate the number of carbon atoms to be functionalized based on the ratio
    num_carbons_to_functionalize = int(len(carbon_indices) * functionalization_ratio)
    
    # Randomly select carbon atoms to be functionalized
    selected_carbons = random.sample(carbon_indices, num_carbons_to_functionalize)
    
    new_positions = positions.tolist()
    new_symbols = symbols[:]
    
    # Use neighbor list to find hydrogen atoms attached to the selected carbon atoms
    i_list, j_list = neighbor_list('ij', atoms, 1.5, self_interaction=False)
    
    for carbon_index in selected_carbons:
        attached_hydrogens = [j for i, j in zip(i_list, j_list) if i == carbon_index and symbols[j] == 'H']
        
        if len(attached_hydrogens) >= 2:
            # Delete one hydrogen atom
            hydrogen_to_delete = attached_hydrogens[0]
            new_positions[hydrogen_to_delete] = None
            new_symbols[hydrogen_to_delete] = None
            
            # Replace the other hydrogen with an oxygen atom
            hydrogen_to_replace = attached_hydrogens[1]
            new_symbols[hydrogen_to_replace] = 'O'
    
    # Remove None entries from positions and symbols
    new_positions = [pos for pos in new_positions if pos is not None]
    new_symbols = [sym for sym in new_symbols if sym is not None]
    
    # Create a new Atoms object with the modified structure
    new_atoms = Atoms(symbols=new_symbols, positions=new_positions, cell=atoms.get_cell(), pbc=atoms.get_pbc())
    
    return new_atoms

def main():
    parser = argparse.ArgumentParser(description='Process and functionalize a polymer supercell.')
    parser.add_argument('input_file', type=str, help='Path to the input file containing the initial structure.')
    parser.add_argument('output_file', type=str, help='Path to the output file to save the modified structure.')
    parser.add_argument('--transformation_matrix', type=int, nargs=9, default=[3, 0, 0, 0, 3, 0, 0, 0, 10],
                        help='Transformation matrix for the supercell (default: [3, 0, 0, 0, 3, 0, 0, 0, 10]).')
    parser.add_argument('--tolerance', type=float, default=1.0, help='Tolerance value for edge detection (default: 1.0).')
    parser.add_argument('--functionalization_ratio', type=float, default=0.1, help='Ratio of carbon atoms to be functionalized (default: 0.1).')
    
    args = parser.parse_args()
    
    # Read the initial structure
    atom = read(args.input_file)
    
    # Define the transformation matrix for the supercell
    transformation_matrix = np.array(args.transformation_matrix).reshape(3, 3)
    
    # Create the supercell
    supercell = make_supercell(atom, transformation_matrix)
    
    # Remove edge atoms and their attached hydrogen atoms
    clean_supercell = remove_edge_atoms_and_hydrogens(supercell, tolerance=args.tolerance)
    
    # Identify the carbon atoms with the smallest and largest z-coordinate values
    min_z_indices, max_z_indices = identify_extreme_carbons(clean_supercell)
    
    # Add hydrogen atoms to the extreme carbon atoms to convert them from CH2 to CH3
    modified_supercell = add_hydrogens_to_extreme_carbons(clean_supercell, min_z_indices, max_z_indices)
    
    # Assign residues to the modified supercell
    modified_supercell_with_residues = assign_residues(modified_supercell)
    
    # Functionalize the carbon atoms
    functionalized_supercell = functionalize_carbons(modified_supercell_with_residues, functionalization_ratio=args.functionalization_ratio)
    
    # Save the modified structure
    write(args.output_file, functionalized_supercell)

if __name__ == '__main__':
    main()