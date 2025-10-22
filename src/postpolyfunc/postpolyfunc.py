from __future__ import annotations
from pathlib import Path
import logging
from .func import PolymerFunctionalizer
from .ligpargen import generate_parameters
from .utils import keep_only_root_gmx
from .gmx import GmxAPI,StepOutput
from ase.io import read,write
import re

def _read_moleculetype_name(itp_path: Path) -> str:
    txt = Path(itp_path).read_text()
    m = re.search(r'^\s*\[\s*moleculetype\s*\]\s*(?:;.*)?$([\s\S]*?)(?=^\s*\[|\Z)', txt, re.MULTILINE)
    if not m:
        raise ValueError(f"[moleculetype] not found in {itp_path}")
    for line in m.group(1).splitlines():
        s = line.strip()
        if s and not s.startswith((';', '#')):
            return s.split()[0]
    raise ValueError(f"Empty [moleculetype] in {itp_path}")

def _rewrite_gro_resname(gro_path: Path, target_resname: str, only_if_name_in: set[str] | None = None) -> int:
    p = Path(gro_path)
    lines = p.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid .gro: {gro_path}")
    title, natoms = lines[0], int(lines[1].strip())
    atom_lines = lines[2:2+natoms]
    box_line   = lines[2+natoms] if len(lines) >= 3+natoms else ""
    tname = (target_resname[:5]).ljust(5)

    changed, fixed = 0, []
    for L in atom_lines:
        if len(L) < 20:
            fixed.append(L); continue
        resid, resname, atom, atomnr, rest = L[0:5], L[5:10], L[10:15], L[15:20], L[20:]
        cur = resname.strip()
        if (only_if_name_in is None) or (cur in only_if_name_in):
            resid_s  = f"{int(resid):5d}" if resid.strip().isdigit() else resid
            atomnr_s = f"{int(atomnr):5d}" if atomnr.strip().isdigit() else atomnr
            fixed.append(f"{resid_s}{tname}{atom}{atomnr_s}{rest}")
            changed += 1
        else:
            fixed.append(L)

    p.write_text("\n".join([title, f"{natoms}", *fixed, box_line]) + "\n")
    return changed

def _normalize_combined_gro_with_itps(combined_gro: Path, solute_itp: Path, solvent_itp: Path) -> tuple[int, int]:
    solute_name  = _read_moleculetype_name(solute_itp)
    solvent_name = _read_moleculetype_name(solvent_itp)
    # Change MOL → solute_name; also coerce generic waters (if any) → solvent_name
    solvent_like = {"SOL", "WAT", "HOH"}
    n_mol  = _rewrite_gro_resname(combined_gro, solute_name, only_if_name_in={"MOL"})
    n_solv = _rewrite_gro_resname(combined_gro, solvent_name, only_if_name_in=solvent_like)
    return n_mol, n_solv


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
        "solvate":      {"outname": "solvated.gro", "topol_top": topol},  # now it exists
        "prepare_topology": {
            "solvent_itp": outdir / "solvent.gmx.itp",
            "solute_itp":  outdir / "solute.gmx.itp",
            "topol_top":   topol,
            "outdir":      "toppar",
            "solvent_outname": "solvent.itp",
            "solute_outname":  "solute.itp",
        },
    }


    results = gmx.orchestrate(
        steps=["create_box", 
               "solvent_box", 
               "solvate", 
               "prepare_topology",
               "normalize_resnames",
               "normalize_atomnames",
                "normalize_topol"
               ],
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
        
