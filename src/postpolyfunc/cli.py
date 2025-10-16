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


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="postpolyfunc",
        description="Functionalize a polymer, then parameterize solute/solvent via LigParGen.",
    )
    p.add_argument("--version", action="version", version="postpolyfunc 0.1.0")

    # --- Input paths ---
    p.add_argument("--solute", required=True, type=Path, help="Polymer structure file (PDB/MOL2).")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--solvent", type=Path, help="Solvent structure file.")
    group.add_argument("--solvent-smiles", type=str, help="Solvent SMILES string.")
    p.add_argument("--outdir", type=Path, default=Path("outputs"), help="Output directory.")

    # --- Functionalization options ---
    p.add_argument("-r", "--ratio", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--mode", type=str, default="carbonyl", choices=["carbonyl"])

    # --- LigParGen options ---
    p.add_argument("--lp-cgen", type=str, default="CM1A-LBCC", choices=["CM1A", "CM1A-LBCC"])
    p.add_argument("--lp-opt", type=int, default=0, choices=[0,1,2,3])
    p.add_argument("--lp-exe", type=str, default=None, help="Override LigParGen executable name")
    p.add_argument("--solute-charge", type=int, default=0)
    p.add_argument("--solvent-charge", type=int, default=0)

    # --- Logging and flow ---
    p.add_argument("--skip-ligpargen", action="store_true")
    p.add_argument("-v", "--verbose", action="count", default=1)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)
    return run_workflow(args)  # 👈 delegate all logic


if __name__ == "__main__":
    raise SystemExit(main())
