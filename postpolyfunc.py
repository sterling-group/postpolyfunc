#!/usr/bin/env python3
"""
poly_topogen.py
Code for generating polymer + solvent topologies, building a solvated box,
and merging final topologies.

Roadmap:
1) Run LigParGen for polymer + solvent
2) Build/pack a box (Packmol or GROMACS)
3) Combine/clean topologies and coordinate files
4) Emit simulation-ready outputs

Author: Ignacio Migliaro
"""

from __future__ import annotations
import subprocess
import shutil
import logging
from dataclasses import dataclass, field
from pathlib import Path
import ase
from ase.visualize import view
from ase.io import read,write
from ase.neighborlist import neighbor_list
from ase import Atoms
import random
from pathlib import Path
from dataclasses import dataclass, field
import logging
import shutil
import subprocess

# --------------------------
# Logging
# --------------------------
def setup_logging(verbosity: int = 1) -> None:
    level = logging.WARNING if verbosity == 0 else logging.INFO if verbosity == 1 else logging.DEBUG
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

#!/usr/bin/env python3
"""
postpolyfunc.py
CLI entry point for polymer functionalization + topology generation workflow.

Current stage: performs functionalization only (C→C=O sites) using func.PolymerFunctionalizer.
Next stages (LigParGen, topology merge) will plug in below.
"""

import argparse
from pathlib import Path
from ase.io import read, write
from func import PolymerFunctionalizer


def main():
    parser = argparse.ArgumentParser(
        description="Functionalize a polymer and prepare for topology generation."
    )
    parser.add_argument(
        "--solute",
        required=True,
        type=Path,
        help="Path to solute (polymer) structure file (e.g., polymer.pdb).",
    )
    parser.add_argument(
        "--solvent",
        required=True,
        type=Path,
        help="Path to solvent structure file (e.g., ethanol.pdb).",
    )
    parser.add_argument(
        "--functionalization",
        "-f",
        type=str,
        default="carbonyl",
        choices=["carbonyl"],  # extend later: hydroxyl, epoxide, carboxyl
        help="Type of functionalization to apply to the polymer.",
    )
    parser.add_argument(
        "--ratio",
        "-r",
        type=float,
        default=0.1,
        help="Fraction of carbon atoms to functionalize (default: 0.1).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("outputs"),
        help="Output directory (default: ./outputs).",
    )

    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # --- 1. Functionalize polymer (solute) ---
    print(f"[INFO] Reading polymer structure: {args.solute}")
    atoms = read(args.solute)

    print(f"[INFO] Applying {args.functionalization} functionalization "
          f"to {args.solute.name} (ratio={args.ratio})")
    f = PolymerFunctionalizer(
        functionalization_ratio=args.ratio,
        seed=args.seed,
        mode=args.functionalization,
    )
    func_atoms = f.functionalize_carbons(atoms)

    func_path = args.outdir / f"{args.solute.stem}_func.pdb"
    write(func_path, func_atoms)
    print(f"[INFO] Functionalized polymer saved to {func_path}")

    # --- 2. Solvent stub (placeholder) ---
    print(f"[INFO] Solvent provided: {args.solvent}")
    # Later: run LigParGen for solvent and solute, pack box, merge topologies

    print("[INFO] Functionalization complete. "
          "Topology generation steps will follow in later integration.")


if __name__ == "__main__":
    main()

    
    

