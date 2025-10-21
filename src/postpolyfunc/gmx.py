# gmx.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Union, List, Tuple
import shutil
import gmxapi as gmx
import logging
import time
import shlex
import subprocess
# ──────────────────────────
# Internal GROMACS helpers
# ──────────────────────────

def _read_box_xyz_from_gro(gro_path: str) -> tuple[float, float, float]:
    """
    Parse box vectors (nm) from last line of a .gro file.
    For orthorhombic boxes expect 3 floats; for triclinic, we take first 3.
    """
    last = Path(gro_path).read_text().strip().splitlines()[-1].split()
    if len(last) < 3:
        raise ValueError(f"Could not parse box from {gro_path}")
    return float(last[0]), float(last[1]), float(last[2])


def _count_residues_in_gro(gro_path: str, resname: str) -> int:
    """
    Count residues by residue number for a given residue name in a .gro.
    Assumes classic .gro formatting; uses columns [0:5]=resnr, [5:10]=resname.
    """
    with open(gro_path, "r") as f:
        lines = f.readlines()
    if len(lines) < 3:
        return 0
    natoms = int(lines[1].strip())
    seen: set[int] = set()
    for line in lines[2:2 + natoms]:
        rn = line[5:10].strip()
        if rn == resname:
            try:
                resnr = int(line[:5])
                seen.add(resnr)
            except ValueError:
                # tolerate malformed lines
                continue
    return len(seen)


def write_merged_topology(
    out_top: Path,
    *,
    forcefield_includes: list[str],
    solute_itp: Path,
    solvent_itp: Path,
    solute_name: str,
    solvent_name: str,
    solvated_gro: Path,
    solvent_resname: str = "SOL",
    n_solute: int = 1,
) -> Path:
    """
    Minimal, reproducible merged topology that includes both ITPs and writes a [ molecules ] table.
    Counts solvent molecules directly from the solvated .gro.
    """
    ns = _count_residues_in_gro(str(solvated_gro), solvent_resname)

    lines: list[str] = []
    lines += [*forcefield_includes, ""]
    lines += [f'#include "{solute_itp.name}"', f'#include "{solvent_itp.name}"', ""]
    lines += ["[ system ]", "Mixed solute + solvent", ""]
    lines += ["[ molecules ]", f"{solute_name}    {n_solute}", f"{solvent_name}   {ns}", ""]
    out_top.write_text("\n".join(lines))
    return out_top


# ──────────────────────────
# Orchestration types
# ──────────────────────────

@dataclass
class StepOutput:
    name: str
    files: Dict[str, str] = field(default_factory=dict)  # logical key -> absolute path
    meta: Dict[str, Any] = field(default_factory=dict)


# ──────────────────────────
# Main API
# ──────────────────────────

class GmxAPI:
    """
    Thin, explicit wrapper around gmxapi with a small step registry and an orchestrator.

    Design principles:
      • Filenames are internal/stable; only physical params (box, nmol, scale) come from CLI (via args).
      • Each step returns a StepOutput with realized file paths and optional metadata.
      • Orchestrate by step names; downstream steps receive upstream outputs via context promotion.
    """

    def __init__(
        self,
        executable: str = "gmx_mpi",
        workdir: Union[str, Path, None] = None,
        args: Any = None,  # argparse.Namespace-like; may hold .box, .scale, etc.
    ):
        self.executable = executable
        self.workdir = Path(workdir) if workdir else Path.cwd()
        self.args = args
        self.workdir.mkdir(parents=True, exist_ok=True)

        # Step registry: map name -> bound method
        self._registry: Dict[str, Any] = {
            "create_box": self.create_boxed_structure_step,  # solute-only box (editconf)
            "solvent_box": self.solvent_box_step,            # build pure solvent box sized to solute box
            "solvate_step": self.solvate_step,              # combine solute + solvent boxes (solvate)
            # "grompp_em": self.grompp_em_step,                # EM preproc
            # "mdrun_em": self.mdrun_em_step,                  # EM run
            # "grompp_nvt": self.grompp_nvt_step,              # NVT preproc
            # "mdrun_nvt": self.mdrun_nvt_step,                # NVT run
            # "combine_solvate": self.combine_solvate_step,    # merge solute + solvent via gmx solvate
        }

    # ── steps ─────────────────────────────────────────────────────────────────

    def create_boxed_structure_step(
    self,
    *,
    input_pdb: str | Path,
    box: list[float] | tuple[float, float, float] | None = None,
    boxtype: str = "cubic",
    center: bool = False,                 # upstream centers; leave False unless needed
    outname: str = "solute_boxed.gro",
) -> StepOutput:
       

        in_pdb = str(Path(input_pdb).expanduser().resolve())
        out_path = str((self.workdir / outname).resolve())

        # Resolve box
        if box is None and self.args is not None:
            box = getattr(self.args, "box", None)
        if box is None:
            box = (12.0, 12.0, 12.0)
        bx, by, bz = (str(float(box[0])), str(float(box[1])), str(float(box[2])))

        # Build arguments (exactly like your working snippet)
        args = ["editconf"]
        if center:
            args += ["-c"]
        # "-d 0.0" is optional when using -box; include only if you want it:
        # args += ["-d", "0.0"]
        args += ["-bt", boxtype, "-box", bx, by, bz]

        editconf_op = gmx.commandline_operation(
            self.executable,                  # e.g., "gmx_mpi"
            args,
            input_files={"-f": in_pdb},
            output_files={"-o": out_path},
        )
        editconf_op.run()

        # Where did gmxapi actually write it?
        try:
            produced = Path(editconf_op.output.file["-o"].result()).resolve()
        except Exception:
            produced = Path(out_path)

        target = Path(out_path)
        if produced != target and produced.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(produced), str(target))

        return StepOutput(
            name="create_box",
            files={"gro": str(target)},
            meta={"box": (float(bx), float(by), float(bz)), "boxtype": boxtype, "center": center},
        )

    def solvent_box_step(
    self,
    *,
    # NEW preferred inputs
    box: Optional[Tuple[float, float, float]] = None,
    nmol: int,
    scale: float = 0.33,
    outname: str = "solvent_box.gro",
    solvent_resname: str = "SOL",
    # LEGACY inputs (kept for backward-compat with orchestrate())
    solvent_gro: Optional[Union[str, Path]] = None,       # ignored, kept for compatibility
    solute_box_gro: Optional[Union[str, Path]] = None,    # if provided, read box from last line
) -> StepOutput:
        """
        Build a pure solvent box with `gmx_mpi insert-molecules`.

        Accepts either:
        - NEW: box=(lx,ly,lz)
        - LEGACY: solute_box_gro=<.gro with box at last line>  (box parsed from file)

        Always uses the hardcoded template: <workdir>/solvent.gmx.gro

        CLI equivalent:
        gmx_mpi insert-molecules -nmol {nmol} -ci solvent.gmx.gro -scale {scale} \
            -o {outname} -box {lx} {ly} {lz}
        """
        workdir = Path(self.workdir).expanduser().resolve()
        workdir.mkdir(parents=True, exist_ok=True)

        # Hardcoded solvent template produced by the previous step
        solvent_ci = (workdir / "solvent.gmx.gro").resolve()
        if not solvent_ci.exists():
            raise FileNotFoundError(f"Required solvent template not found: {solvent_ci}")

        # Resolve box from either new arg or legacy file
        if box is not None:
            lx, ly, lz = (float(box[0]), float(box[1]), float(box[2]))
        else:
            if solute_box_gro is None:
                raise TypeError("Either `box` must be provided or `solute_box_gro` must point to a .gro with box.")
            solute_box_gro = Path(solute_box_gro).expanduser().resolve()
            last = solute_box_gro.read_text().strip().splitlines()[-1].split()
            lx, ly, lz = (float(last[0]), float(last[1]), float(last[2]))

        out_path = str((workdir / outname).resolve())

        # gmx insert-molecules args
        args = [
            "insert-molecules",
            "-nmol", str(int(nmol)),
            "-ci", str(solvent_ci),
            "-scale", str(float(scale)),
            "-box", str(lx), str(ly), str(lz),
        ]

        logging.info(
            "Building pure solvent box (insert-molecules): %s %s",
            self.executable, " ".join(shlex.quote(a) for a in args)
        )

        # Use keyword args to avoid signature mismatch
        op = gmx.commandline_operation(
            executable=self.executable,      # e.g., "gmx_mpi"
            arguments=args,
            input_files={},                   # none
            output_files={"-o": out_path},
        )
        op.run()
        produced = Path(op.output.file["-o"].result()).resolve()

        # Count actual inserted molecules by unique residue IDs with given resname
        nmol_actual = 0
        try:
            with produced.open() as fh:
                lines = fh.readlines()
            if len(lines) >= 3:
                natoms = int(lines[1].strip())
                seen_resids = set()
                # GRO columns: resid (0:5), resname (5:10)
                for line in lines[2:2 + natoms]:
                    if line[5:10].strip() == solvent_resname:
                        resid = line[:5].strip()
                        if resid:
                            seen_resids.add(int(resid))
                nmol_actual = len(seen_resids)
        except Exception as e:
            logging.warning(f"Could not count molecules from {produced.name}: {e}")

        return StepOutput(
            name="solvent_box",
            files={"gro": str(produced)},
            meta={
                "box": (lx, ly, lz),
                "scale": float(scale),
                "nmol_target": int(nmol),
                "nmol_actual": int(nmol_actual),
            },
        )

    def solvate_step(
    self,
    *,
    solute_boxed_gro: str | Path = "solute_boxed.gro",
    solvent_box_gro: str | Path = "solvent_box.gro",
    outname: str = "combined_test.gro",
    solvent_resname: str = "SOL",   # used only for counting
) -> StepOutput:
        """
        Combine a boxed solute with a solvent box using `gmx_mpi solvate`.

        Equivalent CLI:
        gmx_mpi solvate -cp solute_boxed.gro -cs solvent_box.gro -o combined_test.gro
        """
        workdir = Path(self.workdir).expanduser().resolve()
        workdir.mkdir(parents=True, exist_ok=True)

        # Resolve inputs under workdir unless absolute paths are given
        solute_gro = Path(solute_boxed_gro)
        if not solute_gro.is_absolute():
            solute_gro = workdir / solute_gro
        solvent_gro = Path(solvent_box_gro)
        if not solvent_gro.is_absolute():
            solvent_gro = workdir / solvent_gro
        out_path = str((workdir / outname).resolve())

        if not solute_gro.exists():
            raise FileNotFoundError(f"Solute box not found: {solute_gro}")
        if not solvent_gro.exists():
            raise FileNotFoundError(f"Solvent box not found: {solvent_gro}")

        # Build args (note: DO NOT include '-o' here; gmxapi sets it via output_files)
        args = [
            "solvate",
            "-cp", str(solute_gro),
            "-cs", str(solvent_gro),
        ]

        logging.info(
            "Solvating (combine polymer + solvent): %s %s -o %s",
            self.executable,
            " ".join(shlex.quote(a) for a in args),
            shlex.quote(out_path),
        )

        op = gmx.commandline_operation(
            executable=self.executable,      # e.g., "gmx_mpi"
            arguments=args,
            input_files={
                "-cp": str(solute_gro),
                "-cs": str(solvent_gro),
            },
            output_files={"-o": out_path},
        )

        rc = op.run()

        combined = Path(out_path)
        if not combined.exists():
            try:
                stdout = op.output.stdout.result()
            except Exception:
                stdout = ""
            try:
                stderr = op.output.stderr.result()
            except Exception:
                stderr = ""
            raise RuntimeError(
                "gmx solvate failed (no output file created).\n"
                f"Return code: {rc}\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}"
            )

        # Optional: count solvent molecules (unique residue ids with given resname)
        nmol_solvent = 0
        try:
            with combined.open() as fh:
                lines = fh.readlines()
            natoms = int(lines[1].strip()) if len(lines) >= 3 else 0
            seen = set()
            for line in lines[2:2 + natoms]:
                if line[5:10].strip() == solvent_resname:
                    resid = line[:5].strip()
                    if resid:
                        seen.add(int(resid))
            nmol_solvent = len(seen)
        except Exception as e:
            logging.warning(f"Could not count solvent molecules in {combined.name}: {e}")

        return StepOutput(
            name="solvate",
            files={"gro": str(combined)},
            meta={"solvent_resname": solvent_resname, "nmol_solvent": int(nmol_solvent)},
        )

    def grompp_em_step(
    self,
    *,
    gro: str | Path,
    top: str | Path,
    mdp: str | Path,
    out_tpr: str = "em.tpr",
) -> StepOutput:
     
        gro = str(Path(gro).expanduser().resolve())
        top = str(Path(top).expanduser().resolve())
        mdp = str(Path(mdp).expanduser().resolve())
        out_tpr = str((self.workdir / out_tpr).resolve())

        op = gmx.commandline_operation(
            self.executable,
            ["grompp"],
            input_files={"-f": mdp, "-c": gro, "-p": top},
            output_files={"-o": out_tpr},
        )
        logging.info(f"Running grompp for EM: {self.executable} grompp -f {mdp} -c {gro} -p {top} -o {out_tpr}")
        op.run()
        produced = Path(op.output.file["-o"].result()).resolve()

        return StepOutput(name="grompp_em", files={"tpr": str(produced)}, meta={})

    def grompp_em_step(
    self,
    *,
    gro: str | Path,
    top: str | Path,
    mdp: str | Path,
    out_tpr: str = "em.tpr",
) -> StepOutput:
    
        gro = str(Path(gro).expanduser().resolve())
        top = str(Path(top).expanduser().resolve())
        mdp = str(Path(mdp).expanduser().resolve())
        out_tpr = str((self.workdir / out_tpr).resolve())

        op = gmx.commandline_operation(
            self.executable,
            ["grompp"],
            input_files={"-f": mdp, "-c": gro, "-p": top},
            output_files={"-o": out_tpr},
        )
        logging.info(f"Running grompp for EM: {self.executable} grompp -f {mdp} -c {gro} -p {top} -o {out_tpr}")
        op.run()
        produced = Path(op.output.file["-o"].result()).resolve()

        return StepOutput(name="grompp_em", files={"tpr": str(produced)}, meta={})

    def mdrun_em_step(
    self,
    *,
    tpr: str | Path,
    deffnm: str = "em",
) -> StepOutput:
        from pathlib import Path
        import gmxapi as gmx

        tpr = str(Path(tpr).expanduser().resolve())
        # Tell gmxapi where we expect the primary outputs to land
        out_gro = str((self.workdir / f"{deffnm}.gro").resolve())
        out_edr = str((self.workdir / f"{deffnm}.edr").resolve())
        out_log = str((self.workdir / f"{deffnm}.log").resolve())
        out_trr = str((self.workdir / f"{deffnm}.trr").resolve())

        args = ["mdrun", "-deffnm", deffnm]
        op = gmx.commandline_operation(
            self.executable,
            args,
            input_files={"-s": tpr},
            output_files={"-c": out_gro, "-e": out_edr, "-g": out_log, "-o": out_trr},
        )
        logging.info(f"Running mdrun for EM: {self.executable} {' '.join(shlex.quote(a) for a in args)}")
        op.run()

        # Collect outputs (may already be exactly at those paths)
        out = op.output.file
        return StepOutput(
            name="mdrun_em",
            files={
                "gro": str(Path(out["-c"].result()).resolve()),
                "edr": str(Path(out["-e"].result()).resolve()),
                "log": str(Path(out["-g"].result()).resolve()),
                "trr": str(Path(out["-o"].result()).resolve()),
            },
            meta={"deffnm": deffnm},
        )

    def grompp_nvt_step(
    self,
    *,
    gro: str | Path,
    top: str | Path,
    mdp: str | Path,
    out_tpr: str = "nvt.tpr",
) -> StepOutput:
        from pathlib import Path
        import gmxapi as gmx

        gro = str(Path(gro).expanduser().resolve())
        top = str(Path(top).expanduser().resolve())
        mdp = str(Path(mdp).expanduser().resolve())
        out_tpr = str((self.workdir / out_tpr).resolve())

        op = gmx.commandline_operation(
            self.executable,
            ["grompp"],
            input_files={"-f": mdp, "-c": gro, "-p": top},
            output_files={"-o": out_tpr},
        )
        op.run()
        produced = Path(op.output.file["-o"].result()).resolve()
        return StepOutput(name="grompp_nvt", files={"tpr": str(produced)}, meta={})

    def mdrun_nvt_step(
    self,
    *,
    tpr: str | Path,
    deffnm: str = "nvt",
) -> StepOutput:
      

        tpr = str(Path(tpr).expanduser().resolve())
        out_gro = str((self.workdir / f"{deffnm}.gro").resolve())
        out_edr = str((self.workdir / f"{deffnm}.edr").resolve())
        out_log = str((self.workdir / f"{deffnm}.log").resolve())
        out_trr = str((self.workdir / f"{deffnm}.trr").resolve())

        args = ["mdrun", "-deffnm", deffnm]
        op = gmx.commandline_operation(
            self.executable,
            args,
            input_files={"-s": tpr},
            output_files={"-c": out_gro, "-e": out_edr, "-g": out_log, "-o": out_trr},
        )
        op.run()

        out = op.output.file
        return StepOutput(
            name="mdrun_nvt",
            files={
                "gro": str(Path(out["-c"].result()).resolve()),
                "edr": str(Path(out["-e"].result()).resolve()),
                "log": str(Path(out["-g"].result()).resolve()),
                "trr": str(Path(out["-o"].result()).resolve()),
            },
            meta={"deffnm": deffnm},
        )

    def combine_solvate_step(
    self,
    *,
    solute_gro: str | Path,
    solvent_gro: str | Path,
    outname: str = "solvated.gro",
    scale: float = 0.57,
) -> StepOutput:
        from pathlib import Path
        import gmxapi as gmx

        solute_gro = str(Path(solute_gro).expanduser().resolve())
        solvent_gro = str(Path(solvent_gro).expanduser().resolve())
        out_path = str((self.workdir / outname).resolve())

        args = ["solvate", "-cp", solute_gro, "-cs", solvent_gro, "-scale", str(scale)]
        op = gmx.commandline_operation(
            self.executable,
            args,
            input_files={},                      # none
            output_files={"-o": out_path},
        )
        op.run()
        produced = Path(op.output.file["-o"].result()).resolve()

        return StepOutput(
            name="combine_solvate",
            files={"gro": str(produced)},
            meta={"scale": float(scale)},
        )

    # ── orchestrator ───────────────────────────────────────────────────────────

    def orchestrate(
        self,
        steps: Iterable[str],
        *,
        inputs: Dict[str, Any],
        overrides: Optional[Dict[str, Dict[str, Any]]] = None,
        strict: bool = True,
    ) -> Dict[str, StepOutput]:
        """
        Execute named steps in order. Promotes common outputs into context for downstream steps.
        """
        ctx: Dict[str, Any] = dict(inputs)
        results: Dict[str, StepOutput] = {}
        overrides = overrides or {}

        for name in steps:
            fn = self._registry.get(name)
            if fn is None:
                if strict:
                    raise KeyError(f"Unknown step: {name}")
                print(f"[WARN] Skipping unknown step: {name}")
                continue

            # Merge context + per-step overrides and call the step
            kwargs = {**ctx, **overrides.get(name, {})}
            out: StepOutput = fn(**kwargs)  # type: ignore[misc]
            results[name] = out

            # Promote canonical outputs into ctx for typical chains
            if name == "create_box":
                ctx["solute_box_gro"] = out.files["gro"]
                ctx["gro"] = out.files["gro"]
            elif name == "solvent_box":
                ctx["gro"] = out.files["gro"]
            elif name in ("mdrun_em", "mdrun_nvt"):
                if out.files.get("gro"):
                    ctx["gro"] = out.files["gro"]
            elif name == "combine_solvate":
                ctx["gro"] = out.files["gro"]

            # Stash metadata (optional)
            ctx[f"{name}__meta"] = out.meta

        return results
