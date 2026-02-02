#!/usr/bin/env python3
import argparse
import subprocess
import random
from typing import List, Set


def build_functionalized_smiles(
    n: int,
    func_ratio: float,
    seed: int | None,
    allow_terminal: bool
) -> tuple[str, List[int]]:
    if not (0.0 <= func_ratio <= 1.0):
        raise ValueError("--func-ratio must be between 0.0 and 1.0")

    rng = random.Random(seed)

    if allow_terminal:
        eligible = list(range(1, n + 1))          # 1-indexed
    else:
        eligible = list(range(2, n)) if n >= 3 else []

    k = int(round(func_ratio * n))
    k = min(k, len(eligible))

    chosen: Set[int] = set(rng.sample(eligible, k)) if k > 0 else set()
    chosen_sorted = sorted(chosen)

    tokens: List[str] = []
    for i in range(1, n + 1):
        tokens.append("C(=O)" if i in chosen else "C")

    smiles = "".join(tokens)
    return smiles, chosen_sorted


def main():
    parser = argparse.ArgumentParser(
        description="Generate a linear chain with optional random C=O functionalization using Open Babel."
    )
    parser.add_argument("n", type=int, help="Number of carbons in the backbone (e.g., 20).")

    parser.add_argument(
        "--format",
        choices=["pdb", "cml", "mol2"],
        default="pdb",
        help="Output format (default: pdb)."
    )
    parser.add_argument(
        "-o", "--outfile", type=str, default=None,
        help="Output filename. If not set, auto-generated with the right extension."
    )

    parser.add_argument(
        "--func-ratio", type=float, default=0.0,
        help="Fraction of backbone carbons to convert into carbonyls as C(=O). Example: 0.10."
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducibility (optional)."
    )
    parser.add_argument(
        "--allow-terminal", action="store_true",
        help="Allow carbonyls at terminal carbons (creates aldehyde ends). Default: internal only."
    )

    args = parser.parse_args()

    n = args.n
    if n < 1:
        raise ValueError("n must be ≥ 1")

    smiles, chosen_positions = build_functionalized_smiles(
        n=n,
        func_ratio=args.func_ratio,
        seed=args.seed,
        allow_terminal=args.allow_terminal,
    )

    ratio_tag = int(round(args.func_ratio * 100))

    # Auto outfile name if not provided
    if args.outfile is None:
        outfile = f"C{n}_{ratio_tag}.{args.format}"
    else:
        outfile = args.outfile

    # If user gave a name without an extension, add one
    if "." not in outfile.split("/")[-1]:
        outfile = f"{outfile}.{args.format}"

    # Explicitly set output format (-opdb/-ocml/-omol2) so it matches args.format
    cmd = ["obabel", f"-:{smiles}", f"-o{args.format}", "-O", outfile, "--gen3D", "-h"]

    print(f"[INFO] Backbone n={n}")
    print(f"[INFO] Functionalization ratio={args.func_ratio:.3f} (seed={args.seed})")
    print(f"[INFO] Carbonyl positions (1-indexed): {chosen_positions if chosen_positions else 'None'}")
    print(f"[INFO] SMILES: {smiles}")
    print(f"[INFO] Running: {' '.join(cmd)}")

    subprocess.run(cmd, check=True)

    print(f"[OK] Wrote {outfile}")


if __name__ == "__main__":
    main()

