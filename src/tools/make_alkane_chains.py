#!/usr/bin/env python3
import argparse
import subprocess
import random
from typing import List, Set


def build_functionalized_smiles(
    n: int,
    func_ratio: float,
    seed: int | None,
    allow_terminal: bool,
    func_type: str,
    oh_frac: float,
) -> tuple[str, List[int], List[int], List[int]]:
    """
    Build SMILES for a linear carbon backbone where selected positions are functionalized.

    Functional group options:
      carbonyl   -> backbone carbon becomes C(=O)
      hydroxyl   -> backbone carbon becomes C(O)
      benzamide  -> backbone carbon gets a benzamide-substituted phenyl side group
                    C(c1ccc(C(=O)N)cc1)
      mixed      -> selected sites are split between hydroxyl and carbonyl

    Returns:
      smiles,
      carbonyl_positions,
      hydroxyl_positions,
      benzamide_positions
    """
    if not (0.0 <= func_ratio <= 1.0):
        raise ValueError("--func-ratio must be between 0.0 and 1.0")

    if func_type not in {"carbonyl", "hydroxyl", "benzamide", "mixed"}:
        raise ValueError("--func-type must be one of: carbonyl, hydroxyl, benzamide, mixed")

    if not (0.0 <= oh_frac <= 1.0):
        raise ValueError("--oh-frac must be between 0.0 and 1.0")

    rng = random.Random(seed)

    # 1-indexed positions along the backbone
    if allow_terminal:
        eligible = list(range(1, n + 1))
    else:
        eligible = list(range(2, n)) if n >= 3 else []

    k = int(round(func_ratio * n))
    k = min(k, len(eligible))

    chosen: Set[int] = set(rng.sample(eligible, k)) if k > 0 else set()
    chosen_sorted = sorted(chosen)

    carbonyl_positions: List[int] = []
    hydroxyl_positions: List[int] = []
    benzamide_positions: List[int] = []

    if func_type == "carbonyl":
        carbonyl_positions = chosen_sorted

    elif func_type == "hydroxyl":
        hydroxyl_positions = chosen_sorted

    elif func_type == "benzamide":
        benzamide_positions = chosen_sorted

    else:
        # mixed: split chosen into OH vs carbonyl
        # NOTE: mixed currently only mixes hydroxyl and carbonyl,
        # matching your original behavior.
        n_oh = int(round(oh_frac * len(chosen_sorted)))
        n_oh = min(n_oh, len(chosen_sorted))

        oh_set = set(rng.sample(chosen_sorted, n_oh)) if n_oh > 0 else set()

        hydroxyl_positions = sorted(oh_set)
        carbonyl_positions = sorted(set(chosen_sorted) - oh_set)

    tokens: List[str] = []

    for i in range(1, n + 1):
        if i in carbonyl_positions:
            tokens.append("C(=O)")

        elif i in hydroxyl_positions:
            tokens.append("C(O)")

        elif i in benzamide_positions:
            # Backbone carbon attached through the amide nitrogen:
            #
            #     backbone-C-NH-C(=O)-phenyl
            #
            # This is N-substituted benzamide.
            #
            tokens.append("C(NC(=O)c1ccccc1)")

        else:
            tokens.append("C")

    smiles = "".join(tokens)
    return smiles, carbonyl_positions, hydroxyl_positions, benzamide_positions


def main():
    parser = argparse.ArgumentParser(
        description="Generate a linear chain with optional random functionalization using Open Babel."
    )

    parser.add_argument(
        "n",
        type=int,
        help="Number of carbons in the backbone, e.g. 20."
    )

    parser.add_argument(
        "--format",
        choices=["pdb", "cml", "mol2"],
        default="pdb",
        help="Output format. Default: pdb."
    )

    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default=None,
        help="Output filename. If not set, auto-generated with the right extension."
    )

    parser.add_argument(
        "--func-ratio",
        type=float,
        default=0.0,
        help="Fraction of backbone positions to functionalize, from 0.0 to 1.0. Example: 0.10."
    )

    parser.add_argument(
        "--func-type",
        choices=["carbonyl", "hydroxyl", "benzamide", "mixed"],
        default="carbonyl",
        help="Functional group type to apply at selected positions."
    )

    parser.add_argument(
        "--oh-frac",
        type=float,
        default=0.5,
        help="Only used if --func-type mixed. Fraction of functionalized sites that become hydroxyl."
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility."
    )

    parser.add_argument(
        "--allow-terminal",
        action="store_true",
        help="Allow functional groups at terminal carbons. Default: internal only."
    )

    args = parser.parse_args()

    n = args.n

    if n < 1:
        raise ValueError("n must be >= 1")

    smiles, carbonyl_pos, hydroxyl_pos, benzamide_pos = build_functionalized_smiles(
        n=n,
        func_ratio=args.func_ratio,
        seed=args.seed,
        allow_terminal=args.allow_terminal,
        func_type=args.func_type,
        oh_frac=args.oh_frac,
    )

    ratio_tag = int(round(args.func_ratio * 100))

    # Auto outfile name if not provided
    if args.outfile is None:
        tag = args.func_type

        if args.func_type == "mixed":
            tag = f"mixed_oh{int(round(args.oh_frac * 100))}"

        outfile = f"C{n}_{ratio_tag}_{tag}.{args.format}"

    else:
        outfile = args.outfile

    # If user gave a name without an extension, add one
    if "." not in outfile.split("/")[-1]:
        outfile = f"{outfile}.{args.format}"

    cmd = [
        "obabel",
        f"-:{smiles}",
        f"-o{args.format}",
        "-O",
        outfile,
        "--gen3D",
        "-h",
    ]

    print(f"[INFO] Backbone n={n}")
    print(f"[INFO] Functionalization ratio={args.func_ratio:.3f} seed={args.seed}")
    print(f"[INFO] func-type={args.func_type}")
    print(f"[INFO] oh-frac={args.oh_frac if args.func_type == 'mixed' else 'n/a'}")
    print(f"[INFO] Carbonyl positions  1-indexed: {carbonyl_pos if carbonyl_pos else 'None'}")
    print(f"[INFO] Hydroxyl positions   1-indexed: {hydroxyl_pos if hydroxyl_pos else 'None'}")
    print(f"[INFO] Benzamide positions  1-indexed: {benzamide_pos if benzamide_pos else 'None'}")
    print(f"[INFO] SMILES: {smiles}")
    print(f"[INFO] Running: {' '.join(cmd)}")

    subprocess.run(cmd, check=True)

    print(f"[OK] Wrote {outfile}")


if __name__ == "__main__":
    main()