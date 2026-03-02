# func.py
import random
from ase import Atoms
from ase.neighborlist import neighbor_list
import math
import time
import numpy as np


class PolymerFunctionalizer:
    """
    Handles polymer functionalization — randomly converts selected CH2 sites
    into functional groups.

    Modes implemented:
      - "carbonyl": CH2 -> C(=O)H (approx via delete one H, replace one H with O)
      - "hydroxyl": CH2 -> CH–OH  (delete one H, replace one H with O, add new H on O)

    Usage:
        >>> from func import PolymerFunctionalizer
        >>> from ase.io import read, write
        >>> atoms = read("polymer.pdb")
        >>> func = PolymerFunctionalizer(functionalization_ratio=0.1, seed=42, mode="hydroxyl")
        >>> new_atoms = func.functionalize_carbons(atoms)
        >>> write("polymer_func.pdb", new_atoms)
    """

    def __init__(self, functionalization_ratio=0.1, seed=42, mode="carbonyl"):
        self.functionalization_ratio = functionalization_ratio
        self.seed = seed
        self.mode = mode

    def functionalize_carbons(self, atoms: Atoms) -> Atoms:
        # ✅ Skip functionalization entirely if ratio is 0.0
        if self.functionalization_ratio == 0.0:
            return atoms.copy()

        if self.mode == "carbonyl":
            return self._carbonyl(atoms)
        elif self.mode == "hydroxyl":
            return self._hydroxyl(atoms)
        else:
            raise NotImplementedError(f"Functionalization mode not implemented: {self.mode}")

    def _select_carbons(self, symbols):
        carbon_indices = [i for i, sym in enumerate(symbols) if sym == "C"]
        num = math.ceil(len(carbon_indices) * self.functionalization_ratio)
        if num <= 0:
            return []
        return random.sample(carbon_indices, min(num, len(carbon_indices)))

    def _carbonyl(self, atoms: Atoms) -> Atoms:
        """
        Selects random carbon atoms, deletes one hydrogen and replaces another with an oxygen.
        """
        random.seed(self.seed if self.seed is not None else int(time.time()))

        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()

        selected_carbons = self._select_carbons(symbols)
        if not selected_carbons:
            return atoms.copy()

        new_positions = positions.tolist()
        new_symbols = symbols[:]

        # neighbor list for H detection
        i_list, j_list = neighbor_list("ij", atoms, 1.5, self_interaction=False)

        for carbon_index in selected_carbons:
            attached_hydrogens = [
                j for i, j in zip(i_list, j_list)
                if i == carbon_index and symbols[j] == "H"
            ]

            if len(attached_hydrogens) >= 2:
                # Delete one hydrogen
                h_del = attached_hydrogens[0]
                new_positions[h_del] = None
                new_symbols[h_del] = None

                # Replace the second with O (at same coordinates)
                h_rep = attached_hydrogens[1]
                new_symbols[h_rep] = "O"

        # Remove deleted atoms
        new_positions = [p for p in new_positions if p is not None]
        new_symbols = [s for s in new_symbols if s is not None]

        return Atoms(
            symbols=new_symbols,
            positions=new_positions,
            cell=atoms.get_cell(),
            pbc=atoms.get_pbc(),
        )

    def _hydroxyl(self, atoms: Atoms) -> Atoms:
        """
        Hydroxyl functionalization:
          CH2 -> CH–OH

        Implementation mirrors carbonyl selection:
          - pick carbon atoms
          - require >=2 attached H
          - delete one H
          - replace one H with O (O placed at the replaced-H position)
          - add a new H bonded to O along the C->O direction (extended away from carbon)
        """
        random.seed(self.seed if self.seed is not None else int(time.time()))

        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()

        selected_carbons = self._select_carbons(symbols)
        if not selected_carbons:
            return atoms.copy()

        new_positions = positions.tolist()
        new_symbols = symbols[:]

        # neighbor list for H detection
        i_list, j_list = neighbor_list("ij", atoms, 1.5, self_interaction=False)

        # Typical O–H bond length (Å)
        OH_BOND = 0.96

        for carbon_index in selected_carbons:
            attached_hydrogens = [
                j for i, j in zip(i_list, j_list)
                if i == carbon_index and symbols[j] == "H"
            ]

            if len(attached_hydrogens) >= 2:
                # Delete one hydrogen (same as carbonyl)
                h_del = attached_hydrogens[0]
                new_positions[h_del] = None
                new_symbols[h_del] = None

                # Replace another hydrogen with oxygen (O sits where that H was)
                h_rep = attached_hydrogens[1]
                new_symbols[h_rep] = "O"

                c_pos = np.array(positions[carbon_index], dtype=float)
                o_pos = np.array(positions[h_rep], dtype=float)

                # Place the new hydroxyl H:
                # direction = away from carbon along the C->O vector
                v = o_pos - c_pos
                vnorm = np.linalg.norm(v)
                if vnorm < 1e-8:
                    # degenerate; skip adding OH H
                    continue
                u = v / vnorm
                h_oh_pos = (o_pos + u * OH_BOND).tolist()

                new_symbols.append("H")
                new_positions.append(h_oh_pos)

        # Remove deleted atoms
        new_positions = [p for p in new_positions if p is not None]
        new_symbols = [s for s in new_symbols if s is not None]

        return Atoms(
            symbols=new_symbols,
            positions=new_positions,
            cell=atoms.get_cell(),
            pbc=atoms.get_pbc(),
        )