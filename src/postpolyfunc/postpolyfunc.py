"""
postpolyfunc.postpolyfunc
=========================
High-level workflow functions that combine:
  1. Functionalization (via func.PolymerFunctionalizer)
  2. (Planned) LigParGen topology generation
  3. (Planned) Solvation and topology merge
"""

from __future__ import annotations
from pathlib import Path
from ase.io import read, write
from .func import PolymerFunctionalizer


def run_functionalization(
    solute_path: Path,
    outdir: Path,
    ratio: float = 0.1,
    seed: int = 42,
    mode: str = "carbonyl",
) -> Path:
    """
    Apply a chosen functionalization to a polymer structure and write output PDB.

    Parameters
    ----------
    solute_path : Path
        Input polymer structure file (.pdb or .xyz).
    outdir : Path
        Directory to write output.
    ratio : float
        Fraction of carbon atoms to functionalize.
    seed : int
        Random seed for reproducibility.
    mode : str
        Functionalization mode (e.g., 'carbonyl', 'hydroxyl', 'epoxide').

    Returns
    -------
    Path : Path
        Path to the functionalized polymer file.
    """
    atoms = read(solute_path)
    func = PolymerFunctionalizer(
        functionalization_ratio=ratio,
        seed=seed,
        mode=mode,
    )
    new_atoms = func.functionalize_carbons(atoms)
    outdir.mkdir(parents=True, exist_ok=True)

    out_path = outdir / f"{solute_path.stem}_func.pdb"
    write(out_path, new_atoms)
    return out_path


# --- placeholders for future workflow stages ---

def run_ligpargen(structure_path: Path) -> None:
    """Stub for future LigParGen integration."""
    raise NotImplementedError("LigParGen interface not yet implemented.")


def run_topology_merge(solute_top: Path, solvent_top: Path, outdir: Path) -> None:
    """Stub for future topology merge logic."""
    raise NotImplementedError("Topology merge not yet implemented.")
