from __future__ import annotations
from pathlib import Path
import logging
from .func import PolymerFunctionalizer
from .ligpargen import generate_parameters
from .utils import keep_only_root_gmx
from .gmx import GmxAPI,StepOutput


def _first_with_suffix(artifacts: dict, *suffixes: str) -> Path | None:
    """
    Return the first file whose suffix matches any of `suffixes`.
    Accepts values that are either str or Path in the artifacts dict.
    """
    # normalize suffixes to lower-case with leading dot (e.g., ".gro")
    wanted = tuple(s.lower() if s.startswith(".") else f".{s.lower()}" for s in suffixes)

    for _grp, files in artifacts.items():
        for p in files:
            pp = Path(p)  # works for both str and Path
            if pp.suffix.lower() in wanted or any(str(pp).lower().endswith(s) for s in wanted):
                return pp
    return None


def _infer_moleculetype_from_itp(itp_path: Path) -> str:
    """Best-effort: read the first [ moleculetype ] name from an .itp (fallback to stem)."""
    try:
        lines = itp_path.read_text().splitlines()
        for i, line in enumerate(lines):
            if line.strip().lower().startswith("[ moleculetype ]"):
                # next non-empty, non-comment line: name  nrexcl
                for j in range(i + 1, min(i + 10, len(lines))):
                    row = lines[j].strip()
                    if row and not row.startswith(("#", ";")):
                        return row.split()[0]
                break
    except Exception:
        pass
    return itp_path.stem


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
    """
    Functionalization → LigParGen (solute + solvent) → box creation → solvation →
    topology preparation.  Stops before minimization or MD.
    """
    
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # 1️⃣ Functionalize polymer (keep exactly as you requested)
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
        resname=args.solvent.stem[:3].upper(),
        molname="solvent",
        ifile=args.solvent,
        charge=args.solvent_charge,
        cgen=args.lp_cgen,
        opt=args.lp_opt,
        executable=args.lp_exe,
    )

    outdir = Path(args.outdir).expanduser().resolve()
    keep_only_root_gmx(outdir)  
    gmx = GmxAPI(workdir=str(outdir), executable=args.gmx)
    gmx.args = args  # so create_boxed_structure_step can read --box if provided

    inputs = {
        "input_pdb": func_path,
        # no nmol/scale here; they’re for the solvent step only
    }

    topol = gmx.ensure_topol_top(path=outdir / "topol.top")   # ← create skeleton if missing

    overrides = {
        "create_box":   {"outname": "solute_boxed.gro"},
        "solvent_box":  {"outname": "solvent_box.gro", "nmol": args.nsolv, "scale": getattr(args, "scale", 0.33),
                        "box": tuple(args.box) if getattr(args, "box", None) else None},
        "solvate":      {"outname": "combined_test.gro", "topol_top": topol},  # now it exists
        "prepare_topology": {
            "solvent_itp": outdir / "solvent.gmx.itp",
            "solute_itp":  outdir / "solute.gmx.itp",
            "topol_top":   topol,
            "outdir":      "toppar",
            "solvent_outname": "dcb.itp",
            "solute_outname":  "c6.itp",
        },
    }


    results = gmx.orchestrate(
        steps=["create_box", "solvent_box", "solvate", "prepare_topology"],
        inputs=inputs,
        overrides=overrides,
    )

    print("[SUMMARY]")
    print(f"  Solute box:      {results['create_box'].files['gro']}")
    print(f"  Solvent box:     {results['solvent_box'].files['gro']}")
    print(f"  Combined system: {results['solvate'].files['gro']}")
    print(f"  Forcefield:      {results['prepare_topology'].files['forcefield_itp']}")
    print(f"  Solvent itp:     {results['prepare_topology'].files['solvent_itp']}")
    print(f"  Solute itp:      {results['prepare_topology'].files['solute_itp']}")
    print(f"  Updated topol:   {results['prepare_topology'].files['topol_top']}")
    print("[INFO] Workflow complete (stopped before minimization).")
    return 0
        
