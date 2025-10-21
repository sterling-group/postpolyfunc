# gmx.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Union, List, Tuple
import shutil
import gmxapi as gmx


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
            "grompp_em": self.grompp_em_step,                # EM preproc
            "mdrun_em": self.mdrun_em_step,                  # EM run
            "grompp_nvt": self.grompp_nvt_step,              # NVT preproc
            "mdrun_nvt": self.mdrun_nvt_step,                # NVT run
            "combine_solvate": self.combine_solvate_step,    # merge solute + solvent via gmx solvate
        }

    # ── low-level runner ───────────────────────────────────────────────────────

    def _run_cmd(
        self,
        *,
        arguments: List[str],
        input_files: Optional[Dict[str, str]] = None,
        output_files: Optional[Dict[str, str]] = None,
    ) -> Dict[str, str]:
        """
        Wrap gmx.commandline_operation with keyword args (correct signature).
        Returns a dict mapping output flags (e.g., '-o') to realized absolute paths.
        """
        op = gmx.commandline_operation(
            command=[self.executable],
            arguments=arguments,
            input_files=input_files or {},
            output_files=output_files or {},
        )
        op.run()

        realized: Dict[str, str] = {}
        for flag in (output_files or {}):
            realized[flag] = str(Path(op.output.file[flag].result()).resolve())
        return realized

    # ── steps ─────────────────────────────────────────────────────────────────

    def create_boxed_structure_step(
        self,
        *,
        input_pdb: str,
        box: Optional[Iterable[float]] = None,
        boxtype: str = "cubic",
        center: bool = False,                  # upstream already centers
        outname: str = "solute_boxed.gro",     # internal
    ) -> StepOutput:
        """Create a box around the (already-centered) solute using editconf."""
        # Resolve box
        if box is None and self.args is not None:
            box = getattr(self.args, "box", None)
        if box is None:
            box = (12.0, 12.0, 12.0)
        box = list(box)
        if len(box) != 3:
            raise ValueError(f"box must be 3 numbers, got {box}")

        out_path = (self.workdir / outname).resolve()

        args = ["editconf"]
        if center:
            args.append("-c")
        args += ["-bt", boxtype, "-box", str(box[0]), str(box[1]), str(box[2])]

        realized = self._run_cmd(
            arguments=args,
            input_files={"-f": str(input_pdb)},
            output_files={"-o": str(out_path)},
        )

        produced = Path(realized["-o"]).resolve()
        if produced != out_path:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(produced), str(out_path))

        return StepOutput(
            name="create_box",
            files={"gro": str(out_path)},
            meta={"box": tuple(float(x) for x in box), "boxtype": boxtype, "center": center},
        )

    def solvent_box_step(
        self,
        *,
        solvent_gro: str,            # LigParGen solvent GRO (single or small box)
        solute_box_gro: str,         # the solute-box .gro to copy box size from
        nmol: int,                   # target # of solvent molecules
        scale: Optional[float] = None,
        solvent_resname: str = "SOL",
        outname: str = "solvent_box.gro",
    ) -> StepOutput:
        """
        Build a pure solvent box with the same dimensions as the solute box.
        Uses `gmx solvate -box ... -cs ... -maxsol N -scale S`.
        """
        lx, ly, lz = _read_box_xyz_from_gro(solute_box_gro)
        if scale is None and self.args is not None:
            scale = getattr(self.args, "scale", None)
        if scale is None:
            scale = 0.57

        out_path = (self.workdir / outname).resolve()
        realized = self._run_cmd(
            arguments=[
                "solvate",
                "-box", str(lx), str(ly), str(lz),
                "-cs", str(solvent_gro),
                "-scale", str(scale),
                "-maxsol", str(nmol),
            ],
            output_files={"-o": str(out_path)},
        )

        produced = realized["-o"]
        nmol_actual = _count_residues_in_gro(produced, solvent_resname)
        return StepOutput(
            name="solvent_box",
            files={"gro": produced},
            meta={"box": (lx, ly, lz), "scale": float(scale), "nmol_target": int(nmol), "nmol_actual": nmol_actual},
        )

    def grompp_em_step(
        self,
        *,
        gro: str,
        top: str,
        mdp: str,
        out_tpr: str = "em.tpr",
    ) -> StepOutput:
        """Preprocess for energy minimization."""
        out_tpr_path = (self.workdir / out_tpr).resolve()
        realized = self._run_cmd(
            arguments=["grompp"],
            input_files={"-f": str(mdp), "-c": str(gro), "-p": str(top)},
            output_files={"-o": str(out_tpr_path)},
        )
        return StepOutput(name="grompp_em", files={"tpr": realized["-o"]})

    def mdrun_em_step(
        self,
        *,
        tpr: str,
        deffnm: str = "em",
    ) -> StepOutput:
        """Run energy minimization with mdrun (-deffnm em)."""
        realized = self._run_cmd(
            arguments=["mdrun", "-deffnm", deffnm],
            input_files={"-s": str(tpr)},
            output_files={
                "-c": str((self.workdir / f"{deffnm}.gro").resolve()),
                "-e": str((self.workdir / f"{deffnm}.edr").resolve()),
                "-g": str((self.workdir / f"{deffnm}.log").resolve()),
                "-o": str((self.workdir / f"{deffnm}.trr").resolve()),
            },
        )
        return StepOutput(
            name="mdrun_em",
            files={"gro": realized.get("-c"), "edr": realized.get("-e"), "log": realized.get("-g"), "trr": realized.get("-o")},
            meta={"deffnm": deffnm},
        )

    def grompp_nvt_step(
        self,
        *,
        gro: str,
        top: str,
        mdp: str,
        out_tpr: str = "nvt.tpr",
    ) -> StepOutput:
        """Preprocess for short NVT equilibration."""
        out_tpr_path = (self.workdir / out_tpr).resolve()
        realized = self._run_cmd(
            arguments=["grompp"],
            input_files={"-f": str(mdp), "-c": str(gro), "-p": str(top)},
            output_files={"-o": str(out_tpr_path)},
        )
        return StepOutput(name="grompp_nvt", files={"tpr": realized["-o"]})

    def mdrun_nvt_step(
        self,
        *,
        tpr: str,
        deffnm: str = "nvt",
    ) -> StepOutput:
        """Run short NVT equilibration with mdrun (-deffnm nvt)."""
        realized = self._run_cmd(
            arguments=["mdrun", "-deffnm", deffnm],
            input_files={"-s": str(tpr)},
            output_files={
                "-c": str((self.workdir / f"{deffnm}.gro").resolve()),
                "-e": str((self.workdir / f"{deffnm}.edr").resolve()),
                "-g": str((self.workdir / f"{deffnm}.log").resolve()),
                "-o": str((self.workdir / f"{deffnm}.trr").resolve()),
            },
        )
        return StepOutput(
            name="mdrun_nvt",
            files={"gro": realized.get("-c"), "edr": realized.get("-e"), "log": realized.get("-g"), "trr": realized.get("-o")},
            meta={"deffnm": deffnm},
        )

    def combine_solvate_step(
        self,
        *,
        solute_gro: str,             # boxed solute
        solvent_gro: str,            # (equilibrated) pure solvent box
        outname: str = "solvated.gro",
        scale: Optional[float] = None,
    ) -> StepOutput:
        """
        Combine solute + solvent using `gmx solvate -cp solute -cs solvent [-scale S]`.
        """
        if scale is None and self.args is not None:
            scale = getattr(self.args, "scale", None)
        if scale is None:
            scale = 0.57

        out_path = (self.workdir / outname).resolve()
        realized = self._run_cmd(
            arguments=["solvate", "-cp", str(solute_gro), "-cs", str(solvent_gro), "-scale", str(scale)],
            output_files={"-o": str(out_path)},
        )
        return StepOutput(
            name="combine_solvate",
            files={"gro": realized["-o"]},
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
