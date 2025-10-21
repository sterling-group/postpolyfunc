"""
Main workflow orchestration for postpolyfunc.
"""

from __future__ import annotations
from pathlib import Path
from .func import PolymerFunctionalizer
from .ligpargen import generate_parameters
from .utils import keep_only_root_gmx, setup_logging

def run_functionalization(solute_path: Path, outdir: Path, ratio: float, seed: int, mode: str) -> Path:
    from ase.io import read, write
    atoms = read(solute_path)
    f = PolymerFunctionalizer(functionalization_ratio=ratio, seed=seed, mode=mode)
    new_atoms = f.functionalize_carbons(atoms)

    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"{solute_path.stem}_func.pdb"
    write(out_path, new_atoms)
    return out_path


def run_workflow(args) -> int:
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Functionalization
    print(f"[INFO] Functionalizing solute: {args.solute}")
    func_path = run_functionalization(
        solute_path=args.solute,
        outdir=outdir,
        ratio=args.ratio,
        seed=args.seed,
        mode=args.mode,
    )
    print(f"[INFO] Functionalized polymer saved to: {func_path}")

    if args.skip_ligpargen:
        print("[INFO] LigParGen skipped.")
        return 0

    # 2️⃣ LigParGen for solute
    print("[INFO] Running LigParGen for functionalized solute...")
    solute_artifacts = generate_parameters(
        workdir=outdir,
        resname=args.solute.stem[:3].upper(),
        molname="solute",
        ifile=func_path,
        charge=args.solute_charge,
        cgen=args.lp_cgen,
        opt=args.lp_opt,
        executable=args.lp_exe,
    )
  
    # 3️⃣ LigParGen for solvent
    print("[INFO] Running LigParGen for solvent...")
    solvent_artifacts = generate_parameters(
        workdir=outdir,
        resname="SOL",
        molname="solvent",
        ifile=args.solvent if args.solvent else None,
        smile=args.solvent_smiles if args.solvent_smiles else None,
        charge=args.solvent_charge,
        cgen=args.lp_cgen,
        opt=args.lp_opt,
        executable=args.lp_exe,
    )
   
    keep_only_root_gmx(outdir)

    # ✅ Summaries
    print("\n[SUMMARY] LigParGen outputs:")
    for group, files in solute_artifacts.items():
        print(f"  [SOLUTE/{group}] {len(files)} files")
    for group, files in solvent_artifacts.items():
        print(f"  [SOLVENT/{group}] {len(files)} files")

    print("\n[INFO] Workflow complete: functionalization + LigParGen parameterization.")
    return 0
