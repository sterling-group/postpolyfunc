# func.py
import random
from ase import Atoms
from ase.neighborlist import neighbor_list
import math


class PolymerFunctionalizer:
    """
    Handles polymer functionalization — randomly converts selected CH2 sites
    into CH–OH by removing one H and replacing another with O.

    Usage:
        >>> from func import PolymerFunctionalizer
        >>> from ase.io import read, write
        >>> atoms = read("polymer.pdb")
        >>> func = PolymerFunctionalizer(functionalization_ratio=0.1, seed=42)
        >>> new_atoms = func.functionalize_carbons(atoms)
        >>> write("polymer_func.pdb", new_atoms)
    """

    def __init__(self, functionalization_ratio=0.1, seed=42,mode="carbonyl"):
        self.functionalization_ratio = functionalization_ratio
        self.seed = seed
        self.mode = mode
        
    def functionalize_carbons(self, atoms: Atoms) -> Atoms:
        # ✅ Skip functionalization entirely if ratio is 0.0
        if self.functionalization_ratio == 0.0:
            return atoms.copy()

        if self.mode == "carbonyl":
            return self._carbonyl(atoms)
        # Future hooks:
        # elif self.mode == "hydroxyl":
        #     return self._hydroxyl(atoms)
        # elif self.mode == "epoxide":
        #     return self._epoxide(atoms)
        else:
            raise NotImplementedError(f"Functionalization mode not implemented: {self.mode}")

        
    def _carbonyl(self, atoms):
        """
        Selects random carbon atoms, deletes one hydrogen and replaces another with an oxygen.

        Parameters:
        atoms (ase.Atoms): The supercell as an ASE Atoms object.

        Returns:
        ase.Atoms: The modified supercell with functionalized carbon atoms.
        """
        random.seed(self.seed)

        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()

        # Identify all carbon atoms
        carbon_indices = [i for i, sym in enumerate(symbols) if sym == 'C']

        # Calculate number of carbons to functionalize
        num_carbons_to_functionalize = math.ceil(len(carbon_indices) * self.functionalization_ratio)
        if num_carbons_to_functionalize == 0:
            return atoms.copy()

        # Randomly select carbons
        selected_carbons = random.sample(carbon_indices, num_carbons_to_functionalize)

        new_positions = positions.tolist()
        new_symbols = symbols[:]

        # Build neighbor list for H detection
        i_list, j_list = neighbor_list('ij', atoms, 1.5, self_interaction=False)

        for carbon_index in selected_carbons:
            attached_hydrogens = [j for i, j in zip(i_list, j_list)
                                  if i == carbon_index and symbols[j] == 'H']

            if len(attached_hydrogens) >= 2:
                # Delete one hydrogen
                hydrogen_to_delete = attached_hydrogens[0]
                new_positions[hydrogen_to_delete] = None
                new_symbols[hydrogen_to_delete] = None

                # Replace the second with an oxygen atom
                hydrogen_to_replace = attached_hydrogens[1]
                new_symbols[hydrogen_to_replace] = 'O'

        # Remove deleted atoms
        new_positions = [p for p in new_positions if p is not None]
        new_symbols = [s for s in new_symbols if s is not None]

        # Create new Atoms object
        new_atoms = Atoms(symbols=new_symbols,
                          positions=new_positions,
                          cell=atoms.get_cell(),
                          pbc=atoms.get_pbc())

        return new_atoms
