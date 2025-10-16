#!/usr/bin/env python3
"""
postpolyfunc CLI
Sequential workflow:
1) Functionalize polymer
2) (future) Run LigParGen
3) (future) Pack solvent box
4) (future) Merge topologies
"""

from __future__ import annotations
import argparse
from pathlib import Path
from .utils import setup_logging
from .postpolyfunc import run_functionalization  # sequential step 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="postpolyfunc",
        description="Run the full post-polymer functionalization workflow.",
    )
    p.add_argument("--version", action="version", version="postpolyfunc 0.1.0")

    # --- Inputs and general options ---
    p.add_argument("--solute", required=True, type=Path,
                   help="Path to polymer structure file (e.g. polymer.pdb)")
    p.add_argument("--solvent", required=False, type=Path,
                   help="Path to solvent structure file (e.g. ethanol.pdb)")
    p.add_argument("--outdir", type=Path, default=Path("outputs"),
                   help="Output directory (default: ./outputs)")
    p.add_argument("-r", "--ratio", type=float, default=0.1,
                   help="Fraction of carbon atoms to functionalize (default: 0.1)")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--mode", type=str, default="carbonyl",
                   choices=["carbonyl"],  # extend later: hydroxyl, epoxide, etc.
                   help="Functionalization type (default: carbonyl)")
    p.add_argument("-v", "--verbose", action="count", default=1,
                   help="Increase output verbosity")

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(args.verbose)

    # --- Step 1: functionalization ---
    print(f"[INFO] Running polymer functionalization: {args.mode}")
    func_path = run_functionalization(
        solute_path=args.solute,
        outdir=args.outdir,
        ratio=args.ratio,
        seed=args.seed,
        mode=args.mode,
    )
    print(f"[INFO] Functionalized polymer saved to {func_path}")

    # --- Step 2: LigParGen stub (future) ---
    if args.solvent:
        print(f"[INFO] Solvent provided: {args.solvent} (LigParGen stage pending)")
    else:
        print("[INFO] No solvent provided. Skipping LigParGen stage.")

    print("[INFO] Workflow completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
