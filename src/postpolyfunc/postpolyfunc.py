from __future__ import annotations
from pathlib import Path
import logging
from .func import PolymerFunctionalizer
from .ligpargen import generate_parameters
from .utils import keep_only_root_gmx
from .gmx import GmxAPI,StepOutput
from ase.io import read,write
import re
import shutil


def _safe_move(src: Path, dst_dir: Path) -> Path:
    """Move src into dst_dir without clobbering. Returns final destination path."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / src.name
    if not dst.exists():
        return Path(shutil.move(str(src), str(dst)))
    stem, suffix = src.stem, "".join(src.suffixes)
    i = 1
    while True:
        cand = dst_dir / f"{stem}.{i}{suffix}"
        if not cand.exists():
            return Path(shutil.move(str(src), str(cand)))
        i += 1

def tidy_outputs(workdir: str | Path) -> dict:
    """
    Organize MD outputs into phase directories and remove gmxapi/backup junk.

    Creates:
      workdir/min/, nvt/, npt/, prod/

    Moves (per phase):
      - Core: <phase>.* (tpr, gro, edr, log, cpt/xtc/trr, mdout.mdp)
      - Analysis: <phase>_*.xvg, <phase>_*.png, <phase>_*.dat
      - Restarts: <phase>_prev.cpt

    Deletes (recursively):
      - Emacs/Vim backups/locks: '#*#', '.*~', '.*.swp', '.#*'
      - Empty gmxapi.commandline.cli*_i0 directories

    Returns summary: {'moved': [...], 'deleted': [...], 'removed_dirs': [...]}
    """
    wd = Path(workdir).resolve()
    if not wd.exists():
        raise FileNotFoundError(f"workdir not found: {wd}")

    summary = {"moved": [], "deleted": [], "removed_dirs": []}

    phase_globs = {
        "min":  ["em.tpr", "em.mdout.mdp", "min.*", "min_*.xvg", "min_*.png", "min_*.dat", "min_prev.cpt"],
        "nvt":  ["nvt.*",  "nvt_*.xvg",    "nvt_*.png",          "nvt_*.dat",  "nvt_prev.cpt"],
        "npt":  ["npt.*",  "npt_*.xvg",    "npt_*.png",          "npt_*.dat",  "npt_prev.cpt"],
        "prod": ["prod.*", "prod_*.xvg",   "prod_*.png",         "prod_*.dat"],
    }

    # 1) Move files into phase directories (only from top-level wd)
    for phase, patterns in phase_globs.items():
        phase_dir = wd / phase
        for pat in patterns:
            for src in wd.glob(pat):
                if src.is_file() and src.parent != phase_dir:
                    try:
                        dst = _safe_move(src, phase_dir)
                        summary["moved"].append({"from": str(src), "to": str(dst)})
                    except Exception as e:
                        summary["moved"].append({"from": str(src), "to": None, "error": str(e)})

    # 2) Delete Emacs/Vim backup/lock files anywhere under workdir
    #    Matches examples like '#em.tpr.1#', '#npt_Total-Energy.xvg.1#', 'file~', '.#lock', '*.swp'
    backup_patterns = [
        re.compile(r"^#.*#$"),        # Emacs auto-save
        re.compile(r".*~$"),          # tilde backups
        re.compile(r"^\.\#.*$"),      # Emacs lock files
        re.compile(r".*\.swp$"),      # Vim swap
        re.compile(r".*\.swo$"),
        re.compile(r".*\.swx$"),
    ]
    for path in wd.rglob("*"):
        if path.is_file():
            name = path.name
            if any(rx.match(name) for rx in backup_patterns):
                try:
                    path.unlink()
                    summary["deleted"].append(str(path))
                except Exception as e:
                    summary["deleted"].append(f"{path} [error: {e}]")

    # 3) Remove empty gmxapi.commandline.cli*_i0 directories (top-level only)
    cli_dir_rx = re.compile(r"^gmxapi\.commandline\.cli\d+_i0$")
    for sub in wd.iterdir():
        if sub.is_dir() and cli_dir_rx.match(sub.name):
            try:
                # remove only if empty
                if not any(sub.iterdir()):
                    sub.rmdir()
                    summary["removed_dirs"].append(str(sub))
            except Exception as e:
                summary["removed_dirs"].append(f"{sub} [error: {e}]")

    return summary

def run_functionalization(solute_path: Path, outdir: Path, ratio: float, seed: int, mode: str) -> Path:
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
    topology prep → EM (grompp + mdrun).
    """
    from pathlib import Path

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Functionalize polymer
    func_path = run_functionalization(
        solute_path=args.solute,
        outdir=outdir,
        ratio=args.ratio,
        seed=args.seed,
        mode=args.mode,
    )
    print(f"[INFO] Functionalized polymer saved to: {func_path}")
    
    # ============================================================
    # Functionalization-only mode: stop here
    # ============================================================
    if getattr(args, "functionalize_only", False):
        print("[INFO] --functionalize-only set. Skipping LigParGen, solvent, and GMX.")
        return 0

    if args.skip_ligpargen:
        print("[INFO] LigParGen skipped.")
        return 0

    # 2) LigParGen: solute
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

    # 3) LigParGen: solvent
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

    # After LigParGen: solvent
    solvent_stem = args.solvent.stem
    src_gro = outdir / f"{solvent_stem}.gmx.gro"
    dst_gro = outdir / "solvent.gmx.gro"
    # if src_gro.exists():
    #     shutil.copy(src_gro, dst_gro)
    # else:
    #     print(f"[ERROR] Expected {src_gro} not found after LigParGen for solvent.")

    src_itp = outdir / f"{solvent_stem}.gmx.itp"
    dst_itp = outdir / "solvent.gmx.itp"
    if src_itp.exists():
        shutil.copy(src_itp, dst_itp)
    else:
        print(f"[ERROR] Expected {src_itp} not found after LigParGen for solvent.")

    keep_only_root_gmx(outdir)

    gmx = GmxAPI(workdir=str(outdir), executable=args.gmx)
    gmx.args = args  # for create_box to read --box if provided

    # Include em_mdp (and nvt_mdp for later) in inputs; grompp will take mdp via overrides
    inputs = {
        "input_pdb": func_path,
        "em_mdp": Path(args.em_mdp),
        "nvt_mdp": Path(args.nvt_mdp),
    }

    topol = gmx.ensure_topol_top(path=outdir / "topol.top")

    overrides = {
    "create_box": {
        "outname": "solute_boxed.gro",
    },
    "solvent_box": {
        "outname": "solvent_box.gro",
        "nmol": args.nsolv,
        "scale": getattr(args, "scale", 0.33),
        "box": tuple(args.box) if getattr(args, "box", None) else None,
    },
    "solvate": {
        "outname": "solvated.gro",
        "topol_top": topol,
    },
    "prepare_topology": {
        "solvent_itp": outdir / "solvent.gmx.itp",
        "solute_itp":  outdir / "solute.gmx.itp",
        "topol_top":   topol,
        "outdir":      "toppar",
        "solvent_outname": "solvent.itp",
        "solute_outname":  "solute.itp",
    },

    # EM
    "grompp:em": {
        "mdp": Path(args.em_mdp),
        "out_tpr": outdir / "em.tpr",
        "mdout_mdp": outdir / "em.mdout.mdp",
        "maxwarn": 1,
    },
    "mdrun:em": {
        "tpr": outdir / "em.tpr",   # defaults to em.tpr anyway; ok to keep explicit
        "deffnm": "min",
        "gpu_id": args.gpu_id,
        "extra_args": ["-pin", "on"],
    },

    # NVT
    "grompp:nvt": {
        "mdp": Path(args.nvt_mdp),
        "out_tpr": outdir / "nvt.tpr",
        "gpu_id": args.gpu_id,
        "mdout_mdp": outdir / "nvt.mdout.mdp",
        "maxwarn": 1,
    },
    "mdrun:nvt": {
        "tpr": outdir / "nvt.tpr",
        "deffnm": "nvt",
        "gpu_id": args.gpu_id,
        "np": 1,
        "extra_args": ["-pin", "on"],
    },

    # Optional NPT
    "grompp:npt": ({
        "mdp": Path(args.npt_mdp),
        "out_tpr": outdir / "npt.tpr",
        "mdout_mdp": outdir / "npt.mdout.mdp",
        "maxwarn": 1,
    } if getattr(args, "npt_mdp", None) else {}),
    "mdrun:npt": ({
        "tpr": outdir / "npt.tpr",
        "deffnm": "npt",
        "np": 1,
        "gpu_id": args.gpu_id,
        "extra_args": ["-pin", "on"],
    } if getattr(args, "npt_mdp", None) else {}),

    # Optional PROD
    "grompp:prod": ({
        "mdp": Path(args.prod_mdp),
        "out_tpr": outdir / "prod.tpr",
        "mdout_mdp": outdir / "prod.mdout.mdp",  
        "maxwarn": 1,
    } if getattr(args, "prod_mdp", None) else {}),
}


    results = gmx.orchestrate(
        steps=[
            "create_box",
            "solvent_box",
            "solvate",
            "prepare_topology",
            "normalize_resnames",
            "normalize_atomnames",
            "normalize_topol",
            "grompp:em",
            "mdrun:em",
            "grompp:nvt",
            "mdrun:nvt",
            "grompp:npt",
            "mdrun:npt",
            "grompp:prod",
        ],
        inputs=inputs,
        overrides=overrides,
    )
    summ = tidy_outputs(outdir)
    print("[SUMMARY]")
    print(f"  Solute box:      {results['create_box'].files['gro']}")
    print(f"  Solvent box:     {results['solvent_box'].files['gro']}")
    print(f"  Combined system: {results['solvate'].files['gro']}")
    print(f"  Forcefield:      {results['prepare_topology'].files['forcefield_itp']}")
    print(f"  Solvent itp:     {results['prepare_topology'].files['solvent_itp']}")
    print(f"  Solute  itp:     {results['prepare_topology'].files['solute_itp']}")
    print(f"  Updated topol:   {results['prepare_topology'].files['topol_top']}")
    if "grompp:em" in results:
        print(f"  EM TPR:          {results['grompp:em'].files['tpr']}")
    if "mdrun:em" in results:
        print(f"  Minimized GRO:   {results['mdrun:em'].files.get('gro', 'N/A')}")
        print(f"  EM log:          {results['mdrun:em'].files.get('log', 'N/A')}")
    
    print("[INFO] Workflow complete (minimization finished).")
    return 0



