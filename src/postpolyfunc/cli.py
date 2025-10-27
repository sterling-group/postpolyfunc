#!/usr/bin/env python3
"""
CLI entry point for postpolyfunc workflow.
Handles CLI parsing, then calls postpolyfunc.run_workflow().
"""
from __future__ import annotations
import argparse
from pathlib import Path
from .utils import setup_logging
from .postpolyfunc import run_workflow


import argparse
from pathlib import Path
import csv

def _ratio(value: str) -> float:
    x = float(value)
    if not (0.0 < x <= 1.0):
        raise argparse.ArgumentTypeError("ratio must be in (0, 1].")
    return x

def _positive_int(value: str) -> int:
    x = int(value)
    if x <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer.")
    return x

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="postpolyfunc",
        description="Functionalize a polymer, parameterize solute/solvent via LigParGen, then prep GMX.",
    )
    p.add_argument("--version", action="version", version="postpolyfunc 0.1.0")

    # --- I/O paths ---
    p.add_argument("--solute", required=True, type=Path,
                   help="Polymer (solute) structure file (e.g., .pdb or .mol2).")
    p.add_argument("--outdir", type=Path, default=Path("outputs"),
                   help="Output directory (default: outputs).")
    p.add_argument("--csv", type=Path, default=None, help="CSV file to override CLI arguments.")

    # --- Solvent definition (choose one) ---
    g_solvent = p.add_mutually_exclusive_group(required=True)
    g_solvent.add_argument("--solvent", type=Path,
                           help="Solvent structure file (e.g., LigParGen-ready .pdb/.mol2).")
    g_solvent.add_argument("--solvent-smiles", type=str,
                           help="Solvent SMILES string (LigParGen will build it).")

    # --- Functionalization options ---
    p.add_argument("-r", "--ratio", type=_ratio, default=0.1,
                   help="Fraction of carbon sites to functionalize, (0,1]. Default: 0.1")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    p.add_argument("--mode", type=str, default="carbonyl",
                   choices=["carbonyl"], help="Functionalization mode (default: carbonyl).")

    # --- LigParGen options ---
    p.add_argument("--lp-cgen", type=str, default="CM1A-LBCC",
                   choices=["CM1A", "CM1A-LBCC"], help="Charge model (default: CM1A-LBCC).")
    p.add_argument("--lp-opt", type=int, default=0, choices=[0, 1, 2, 3],
                   help="Geometry optimization level in LigParGen (default: 0).")
    p.add_argument("--lp-exe", type=str, default=None,
                   help="Override LigParGen executable name (optional).")
    p.add_argument("--solute-charge", type=int, default=0, help="Net charge of the solute (default: 0).")
    p.add_argument("--solvent-charge", type=int, default=0, help="Net charge of the solvent (default: 0).")

    # --- GMX prep (box + solvent packing + short equil) ---
    p.add_argument("--box", type=float, nargs=3, metavar=("X", "Y", "Z"),
                   help="Box dimensions (nm) for the solute (e.g., --box 10 10 12). Default: 12 12 12 if omitted.")
    p.add_argument(
    "--gmx",
    default="gmx_mpi",
    choices=["gmx_mpi", "gmx"],
    help="GROMACS frontend to use (default: gmx_mpi)."
)
    p.add_argument("--nsolv", type=_positive_int, required=True,
                   help="Target number of solvent molecules for the pure solvent box.")
    p.add_argument("--scale", type=float, default=0.33,
                   help="vdW radii scale for packing (gmx solvate -scale). Default: 0.57.")
    p.add_argument("--em-mdp",dest="em_mdp", type=Path, required=True,
                   help="Path to energy-minimization .mdp file.")
    p.add_argument("--nvt-mdp", type=Path, required=True,
                   help="Path to short NVT equilibration .mdp file.")
    p.add_argument("--npt-mdp", dest="npt_mdp", type=Path, required=False,
               help="Path to NPT .mdp file (optional).")
    p.add_argument("--prod-mdp", dest="prod_mdp", type=Path, required=False,
               help="Path to production MD (.mdp). If set, runs the production phase.")

    # --- Logging / flow control ---
    p.add_argument("--skip-ligpargen", action="store_true",
                   help="Skip parameterization (debug/dev only).")
    p.add_argument("-v", "--verbose", action="count", default=1,
                   help="Increase verbosity (-v, -vv).")

    return p

def validate_args(args: argparse.Namespace) -> None:
    # Ensure .mdp files exist early
    if not args.skip_ligpargen:
        # solvent source must be given (parser enforces), nothing to do here
        pass

    if not args.em_mdp.exists():
        raise FileNotFoundError(f"EM .mdp not found: {args.em_mdp}")
    if not args.nvt_mdp.exists():
        raise FileNotFoundError(f"NVT .mdp not found: {args.nvt_mdp}")

    if args.box is not None and len(args.box) != 3:
        raise ValueError("--box must provide exactly 3 numbers (X Y Z in nm).")

def override_args_with_csv(args: argparse.Namespace) -> argparse.Namespace:
    if args.csv is not None:
        with open(args.csv, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            row = next(reader)  # Only first row for now
            # Override CLI args with CSV values if present
            if "solute" in row and row["solute"]:
                args.solute = Path(row["solute"])
            if "solvent" in row and row["solvent"]:
                args.solvent = Path(row["solvent"])
                args.solvent_smiles = None
            if "solvent_smiles" in row and row["solvent_smiles"]:
                args.solvent_smiles = row["solvent_smiles"]
                args.solvent = None
            if "ratio" in row and row["ratio"]:
                args.ratio = float(row["ratio"])
            if "nsolv" in row and row["nsolv"]:
                args.nsolv = int(row["nsolv"])
            if "mode" in row and row["mode"]:
                args.mode = row["mode"]
    return args

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args = override_args_with_csv(args)
    setup_logging(args.verbose)
    return run_workflow(args)  # 👈 delegate all logic


if __name__ == "__main__":
    raise SystemExit(main())
