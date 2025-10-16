#!/usr/bin/env python3
"""
postpolyfunc CLI
- functionalize: apply polymer functionalization and write PDB
(roadmap: ligpargen, packing, merge)
"""

from __future__ import annotations
import argparse
from pathlib import Path
from ase.io import read, write

from .func import functionalize_carbons, PolymerFunctionalizer  # whichever you use
from .utils import setup_logging  # optional
from .postpolyfunc import run_functionalization

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="postpolyfunc",
        description="Functionalize a polymer and (roadmap) generate topologies.",
    )
    p.add_argument("--version", action="version", version="postpolyfunc 0.1.0")

    sub = p.add_subparsers(dest="cmd", required=True)

    # --- functionalize subcommand ---
    f = sub.add_parser("functionalize", help="Apply functionalization to a polymer PDB")
    f.add_argument("--solute", required=True, type=Path, help="Polymer PDB (input)")
    f.add_argument("--out", required=True, type=Path, help="Output PDB path")
    f.add_argument("-r", "--ratio", type=float, default=0.1, help="Fraction of eligible C to functionalize")
    f.add_argument("-n", "--n-sites", type=int, default=None, help="Exactly N sites (overrides --ratio)")
    f.add_argument("--seed", type=int, default=42)
    f.add_argument("-v", "--verbose", action="count", default=1)

    # --- roadmap: pack/ligpargen/merge subcommands later ---
    # sub.add_parser("topogen", help="Run LigParGen, pack box, and merge topologies")

    return p

def cmd_functionalize(args: argparse.Namespace) -> int:
    setup_logging(args.verbose)
    if not args.solute.exists():
        raise SystemExit(f"[ERROR] Input not found: {args.solute}")

    atoms = read(str(args.solute))
    f = PolymerFunctionalizer(
        functionalization_ratio=args.ratio,
        n_sites=args.n_sites,
        seed=args.seed,
        co_bond=args.co_bond,
        oh_bond=args.oh_bond,
        angle_COH_deg=args.angle_coh,
        cutoff_scale=args.cutoff_scale,
    )
    new_atoms = f.functionalize_carbons(atoms)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write(str(args.out), new_atoms)
    return 0

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "functionalize":
        return cmd_functionalize(args)

    parser.error("No subcommand specified")
    return 2

if __name__ == "__main__":
    raise SystemExit(main())
