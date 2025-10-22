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
import inspect


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
            "solvate": self.solvate_step,              # combine solute + solvent boxes (solvate)
            "combine_solvate": self.combine_solvate_step,    # combine solute + solvent boxes (solvate) 
            "prepare_topology": self.prepare_topology_step,  # combine topology files from LigParGen outputs
        }
#TODO "grompp_em": self.grompp_em_step,                # EM preproc
            # "mdrun_em": self.mdrun_em_step,                  # EM run
            # "grompp_nvt": self.grompp_nvt_step,              # NVT preproc
            # "mdrun_nvt": self.mdrun_nvt_step,                # NVT run
            # "combine_solvate": self.combine_solvate_step,    # merge solute + solvent via gmx solvate
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

    def ensure_topol_top(self, *, path: str | Path, system_name: str = "Polymer in solvent") -> Path:
        p = Path(path).expanduser().resolve()
        if p.exists():
            return p
        p.write_text(
            '; auto-generated skeleton\n'
            '#include "./toppar/forcefield.itp"\n'
            '#include "./toppar/dcb.itp"\n'
            '#include "./toppar/c6.itp"\n\n'
            '[ system ]\n'
            f'{system_name}\n\n'
            '[ molecules ]\n'
            '; name  number\n'
        )
        return p


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
    solvent_resname: str = "SOL",
    topol_top: str | Path | None = None,  # optional: let GROMACS update [molecules]
) -> StepOutput:
        workdir = Path(self.workdir).expanduser().resolve()
        workdir.mkdir(parents=True, exist_ok=True)

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

        # Build args (ONLY args; DO NOT duplicate with input_files)
        args = [
            "solvate",
            "-cp", str(solute_gro),
            "-cs", str(solvent_gro),
            
        ]
        # inside solvate_step, around where you handle topol_top
        if topol_top is not None:
            topol_top = Path(topol_top) if Path(topol_top).is_absolute() else workdir / topol_top
            if topol_top.exists():
                args += ["-p", str(topol_top)]
            else:
                logging.warning(f"topol.top not found at {topol_top}; proceeding without -p so [molecules] won't be auto-updated.")


        logging.info("Solvating (combine polymer + solvent): %s %s",
                    self.executable, " ".join(shlex.quote(a) for a in args))

        op = gmx.commandline_operation(
            executable=self.executable,
            arguments=args,
            input_files={},                 # <— IMPORTANT: nothing here to avoid duplicate flags
            output_files={"-o": out_path},  # ok to keep; alternatively remove "-o" from args above
        )
        rc = op.run()

        combined = Path(out_path)
        if not combined.exists():
            try: stdout = op.output.stdout.result()
            except Exception: stdout = ""
            try: stderr = op.output.stderr.result()
            except Exception: stderr = ""
            raise RuntimeError(
                "gmx solvate failed (no output file created).\n"
                f"Return code: {rc}\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}"
            )

        # Optional: count solvent molecules
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

    # Combine topology files from LigParGen outputs
    
    def prepare_topology_step(
        self,
        *,
        solvent_itp: str | Path = "solvent.gmx.itp",
        solute_itp: str | Path = "solute.gmx.itp",
        topol_top: str | Path = "topol.top",
        outdir: str | Path = "toppar",
        solvent_outname: str = "dcb.itp",
        solute_outname: str = "c6.itp",
    ) -> StepOutput:
        """
        Consolidate atomtypes into toppar/forcefield.itp, suffix solvent types with 's'
        and solute (polymer) types with 'p', rewrite [ atoms ] in both .itp files to
        use the suffixed types, and ensure #includes exist in topol.top.

        Outputs:
        toppar/forcefield.itp
        toppar/dcb.itp
        toppar/c6.itp
        """
        import re
        from pathlib import Path

        workdir = Path(self.workdir).expanduser().resolve()
        outdir = (workdir / outdir).resolve()
        outdir.mkdir(parents=True, exist_ok=True)

        solvent_itp = Path(solvent_itp) if Path(solvent_itp).is_absolute() else workdir / solvent_itp
        solute_itp  = Path(solute_itp)  if Path(solute_itp).is_absolute()  else workdir / solute_itp
        topol_top   = Path(topol_top)   if Path(topol_top).is_absolute()   else workdir / topol_top

        if not solvent_itp.exists():
            raise FileNotFoundError(f"Solvent ITP not found: {solvent_itp}")
        if not solute_itp.exists():
            raise FileNotFoundError(f"Solute ITP not found: {solute_itp}")
        if not topol_top.exists():
            raise FileNotFoundError(f"topol.top not found: {topol_top}")

        SECTION_RE = re.compile(r'^\s*\[([^\]]+)\]\s*(?:;.*)?$')
        COMMENT_RE = re.compile(r';.*$')

        def split_sections(text: str) -> dict[str, list[str]]:
            sections: dict[str, list[str]] = {}
            current = None
            for raw in text.splitlines():
                line = raw.rstrip('\n')
                m = SECTION_RE.match(line)
                if m:
                    current = m.group(1).strip().lower()
                    sections.setdefault(current, [])
                else:
                    if current is not None:
                        sections[current].append(line)
            return sections

        def is_data_line(line: str) -> bool:
            s = line.strip()
            return bool(s) and not s.startswith(';') and not s.startswith('#')

        def rename_opls_token(token: str, suffix: str) -> str:
            # only rename opls_### tokens
            return re.sub(r'^(opls_\d+)$', rf'\1{suffix}', token)

        def reatomtype_line(line: str, suffix: str) -> str:
            stripped = COMMENT_RE.sub('', line).strip()
            if not stripped:
                return line
            parts = stripped.split()
            if not parts:
                return line
            parts[0] = rename_opls_token(parts[0], suffix)
            cm = COMMENT_RE.search(line)
            comment = (' ' + cm.group(0)) if cm else ''
            return (' '.join(parts) + comment).rstrip()

        def parse_atomtypes_block(lines: list[str]) -> list[str]:
            return [ln for ln in lines if is_data_line(ln)]

        def rewrite_atoms_section(lines: list[str], suffix: str) -> list[str]:
            out: list[str] = []
            for ln in lines:
                if not is_data_line(ln):
                    out.append(ln)
                    continue
                pre, comment = ln, ""
                m = COMMENT_RE.search(ln)
                if m:
                    comment = m.group(0)
                    pre = ln[:m.start()]
                toks = pre.split()
                # [ atoms ]: nr  type  resnr  resid  atom  cgnr  charge  mass ...
                if len(toks) >= 2:
                    toks[1] = rename_opls_token(toks[1], suffix)
                    new_line = "{:<6} {:<16} {}".format(toks[0], toks[1], " ".join(toks[2:])).rstrip()
                    if comment:
                        new_line += " " + comment
                    out.append(new_line)
                else:
                    out.append(ln)
            return out

        def write_itp(dest: Path, sections: dict[str, list[str]]):
            # include common sections (order isn’t critical, but keep it tidy)
            order = [
                'moleculetype','atoms','bonds','pairs','angles','dihedrals','constraints',
                'exclusions','virtual_sites2','virtual_sites3','virtual_sites4','settles',
                'system','molecules','atomtypes','nonbond_params','bondtypes','angletypes',
                'dihedraltypes','constrainttypes'
            ]
            with dest.open('w') as fh:
                for name in order:
                    if name in sections and sections[name]:
                        fh.write(f"[ {name} ]\n")
                        for ln in sections[name]:
                            fh.write(ln.rstrip() + "\n")
                        fh.write("\n")

        # Load & split
        solv_sec = split_sections(solvent_itp.read_text())
        sol_sec  = split_sections(solute_itp.read_text())
        if 'atomtypes' not in solv_sec or 'atomtypes' not in sol_sec:
            raise RuntimeError("Both solvent and solute ITPs must contain an [ atomtypes ] section.")

        # Build forcefield.itp (defaults + suffixed atomtypes)
        solv_atomtypes = [reatomtype_line(ln, 's') for ln in parse_atomtypes_block(solv_sec['atomtypes'])]
        sol_atomtypes  = [reatomtype_line(ln, 'p') for ln in parse_atomtypes_block(sol_sec['atomtypes'])]

        forcefield_txt = (
            "[ defaults ]\n"
            "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n"
            "1         3          yes        0.5      0.5\n\n"
            "[ atomtypes ]\n" +
            "\n".join(solv_atomtypes) + "\n\n" +
            "\n".join(sol_atomtypes)  + "\n"
        )
        (outdir / "forcefield.itp").write_text(forcefield_txt)

        # Rewrite [ atoms ] in both itps to match suffixed types
        solv_mod = dict(solv_sec)
        if 'atoms' in solv_mod and solv_mod['atoms']:
            solv_mod['atoms'] = rewrite_atoms_section(solv_mod['atoms'], 's')
        sol_mod = dict(sol_sec)
        if 'atoms' in sol_mod and sol_mod['atoms']:
            sol_mod['atoms'] = rewrite_atoms_section(sol_mod['atoms'], 'p')

        # Write updated ITPs with requested names
        solvent_out = outdir / solvent_outname
        solute_out  = outdir / solute_outname
        write_itp(solvent_out, solv_mod)
        write_itp(solute_out,  sol_mod)

        # Ensure #includes in topol.top (idempotent)
        include_lines = [
            '#include "./toppar/forcefield.itp"',
            f'#include "./toppar/{solvent_outname}"',
            f'#include "./toppar/{solute_outname}"',
        ]
        top_txt = topol_top.read_text()
        missing = [L for L in include_lines if L not in top_txt]
        if missing:
            lines = top_txt.splitlines()
            insertion_idx = None
            for i, ln in enumerate(lines):
                s = ln.strip().lower()
                if s.startswith("[ system ]") or s.startswith("[ molecules ]"):
                    insertion_idx = i
                    break
            if insertion_idx is None:
                new_txt = top_txt.rstrip() + "\n\n" + "\n".join(include_lines) + "\n"
            else:
                new_txt = "\n".join(lines[:insertion_idx] + include_lines + [""] + lines[insertion_idx:])
            topol_top.write_text(new_txt)

        return StepOutput(
            name="prepare_topology",
            files={
                "forcefield_itp": str((outdir / "forcefield.itp").resolve()),
                "solvent_itp": str(solvent_out.resolve()),
                "solute_itp":  str(solute_out.resolve()),
                "topol_top":   str(topol_top.resolve()),
            },
            meta={
                "solvent_suffix": "s",
                "solute_suffix": "p",
                "includes_added": bool(missing),
            },
        )

    # ── orchestrator ───────────────────────────────────────────────────────────

    def orchestrate(self, steps, *, inputs, overrides=None, strict=True) -> Dict[str, StepOutput]:
        ctx = dict(inputs)
        results = {}
        overrides = overrides or {}

        for name in steps:
            fn = self._registry.get(name)
            if fn is None:
                if strict:
                    raise KeyError(f"Unknown step: {name}")
                print(f"[WARN] Skipping unknown step: {name}")
                continue

            # merge then filter by function signature
            merged = {**ctx, **overrides.get(name, {})}
            sig = inspect.signature(fn)
            allowed = {k: v for k, v in merged.items() if k in sig.parameters}

            out: StepOutput = fn(**allowed)  # type: ignore[misc]
            results[name] = out

            # promote outputs (your existing logic)
            # promote outputs
            if name == "create_box":
                ctx["solute_box_gro"] = out.files["gro"]
                ctx["gro"] = out.files["gro"]
            elif name == "solvent_box":                 # <-- match the registry key
                ctx["solvent_box_gro"] = out.files["gro"]
                ctx["gro"] = out.files["gro"]
            elif name in ("combine_solvate", "solvate"):
                ctx["gro"] = out.files["gro"]
            elif name in ("mdrun_em", "mdrun_nvt"):
                if out.files.get("gro"):
                    ctx["gro"] = out.files["gro"]

        return results
