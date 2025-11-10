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
import logging

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
    p.add_argument("--solute", required=False, type=Path,
                   help="Polymer (solute) structure file (e.g., .pdb or .mol2).")
    p.add_argument("--outdir", type=Path, default=Path("outputs"),
                   help="Output directory (default: outputs).")
    p.add_argument("--csv", type=Path, default=None, help="CSV file to override CLI arguments.")

    # --- Solvent definition (choose one) ---
    g_solvent = p.add_mutually_exclusive_group(required=False)
    g_solvent.add_argument("--solvent", type=Path,
                           help="Solvent structure file (e.g., LigParGen-ready .pdb/.mol2).")
    g_solvent.add_argument("--solvent-smiles", type=str,
                           help="Solvent SMILES string (LigParGen will build it).")

    # --- Functionalization options ---
    p.add_argument("-r", "--ratio", type=_ratio, default=0.1,
                   help="Fraction of carbon sites to functionalize, (0,1]. Default: 0.1")
    p.add_argument("--seed", type=int, help="Random seed (default: random).")
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
    p.add_argument("--nsolv", type=_positive_int, required=False,
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

def override_args_with_csv(args: argparse.Namespace) -> list[argparse.Namespace] | argparse.Namespace:
    """
    If --csv is given, return a list of Namespace objects (one per row).
    Otherwise, return the original args Namespace.
    """
    if args.csv is not None:
        arglist = []
        with open(args.csv, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                row_args = argparse.Namespace(**vars(args))  # copy

                # Handle solute: path or SMILES
                if "solute" in row and row["solute"]:
                    solute_val = row["solute"]
                    solute_path = Path(solute_val)
                    if solute_path.exists():
                        row_args.solute = solute_path
                        row_args.solute_smiles = None
                    else:
                        row_args.solute = None
                        row_args.solute_smiles = solute_val
                elif "solute_smiles" in row and row["solute_smiles"]:
                    row_args.solute = None
                    row_args.solute_smiles = row["solute_smiles"]
                # Handle box parameters
                if "box" in row and row["box"]:
                    # Allow flexible formats: "12 12 12" or "12,12,12"
                    box_vals = row["box"].replace(",", " ").split()
                    if len(box_vals) != 3:
                        raise ValueError(f"CSV 'box' must have 3 dimensions (got {row['box']})")
                    row_args.box = [float(x) for x in box_vals]

                # Handle solvent: path or SMILES
                if "solvent" in row and row["solvent"]:
                    solvent_val = row["solvent"]
                    solvent_path = Path(solvent_val)
                    if solvent_path.exists():
                        row_args.solvent = solvent_path
                        row_args.solvent_smiles = None
                    else:
                        row_args.solvent = None
                        row_args.solvent_smiles = solvent_val
                elif "solvent_smiles" in row and row["solvent_smiles"]:
                    row_args.solvent = None
                    row_args.solvent_smiles = row["solvent_smiles"]

                if "ratio" in row and row["ratio"]:
                    row_args.ratio = float(row["ratio"])
                if "nsolv" in row and row["nsolv"]:
                    row_args.nsolv = int(row["nsolv"])
                if "mode" in row and row["mode"]:
                    row_args.mode = row["mode"]
                arglist.append(row_args)
        return arglist
    return args

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args_or_list = override_args_with_csv(args)
    setup_logging(args.verbose)
    logger = logging.getLogger("postpolyfunc.batch")

    # Enforce required arguments only if not using CSV
    if args.csv is None:
        if args.solute is None:
            parser.error("--solute is required unless --csv is given.")
        if args.solvent is None and args.solvent_smiles is None:
            parser.error("Either --solvent or --solvent-smiles is required unless --csv is given.")

    if isinstance(args_or_list, list):
        # Batch mode
        exit_codes = []
        for i, row_args in enumerate(args_or_list, 1):
            logger.info(f"=== Batch {i}/{len(args_or_list)} ===")
            # Create a unique output directory for each row
            if getattr(row_args, "solute", None):
                solute_name = Path(row_args.solute).stem
                solute_path = Path(row_args.solute)
                if not solute_path.exists():
                    logger.error(f"Solute file does not exist: {solute_path}")
                    exit_codes.append(1)
                    continue
                logger.info(f"Solute file found: {solute_path}")
            elif getattr(row_args, "solute_smiles", None):
                solute_name = row_args.solute_smiles.replace("/", "_").replace("\\", "_")
                logger.info(f"Solute provided as SMILES: {row_args.solute_smiles}")
            else:
                solute_name = "unknown_solute"
                logger.error("No solute or solute_smiles provided.")
                exit_codes.append(1)
                continue

            if getattr(row_args, "solvent", None):
                solvent_name = Path(row_args.solvent).stem
                solvent_path = Path(row_args.solvent)
                if not solvent_path.exists():
                    logger.error(f"Solvent file does not exist: {solvent_path}")
                    exit_codes.append(1)
                    continue
                logger.info(f"Solvent file found: {solvent_path}")
            elif getattr(row_args, "solvent_smiles", None):
                solvent_name = row_args.solvent_smiles.replace("/", "_").replace("\\", "_")
                logger.info(f"Solvent provided as SMILES: {row_args.solvent_smiles}")
            else:
                solvent_name = "unknown_solvent"
                logger.error("No solvent or solvent_smiles provided.")
                exit_codes.append(1)
                continue

            batch_outdir = Path(args.outdir) / f"batch_{i:02d}_{solute_name}_{solvent_name}"
            row_args.outdir = batch_outdir
            logger.info(f"Output directory: {row_args.outdir}")

            try:
                code = run_workflow(row_args)
                logger.info(f"Batch {i} completed with exit code {code}")
            except Exception as e:
                logger.exception(f"Batch {i} failed with exception: {e}")
                code = 1
            exit_codes.append(code)
        logger.info(f"Batch run complete. Exit codes: {exit_codes}")
        return max(exit_codes)
    else:
        return run_workflow(args_or_list)


if __name__ == "__main__":
    raise SystemExit(main())
