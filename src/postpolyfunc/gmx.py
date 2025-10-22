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
import re
from types import SimpleNamespace


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
            "normalize_resnames": self.normalize_resnames_step,  # fix .gro resnames to match .itp names  
            "normalize_atomnames": self.normalize_atomnames_step,   #fix .atom names to match .itp names
            "normalize_topol": self.normalize_topol_step,      # fix topol.top to match .itp and .gro names

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
            '#include "./toppar/solvent.itp"\n'
            '#include "./toppar/solute.itp"\n\n'
            '[ system ]\n'
            f'{system_name}\n\n'
            '[ molecules ]\n'
            '; name  number\n'
        )
        return p
   
        """
        Rewrite the resname (columns 6–10) of every atom line in a .gro file to target_resname.
        Keeps fixed-width formatting per GROMACS convention.
        """
        lines = gro_path.read_text().splitlines()
        if len(lines) < 3:
            raise ValueError(f"Not a valid .gro: {gro_path}")

        title = lines[0]
        natoms = int(lines[1].strip())
        atom_lines = lines[2:2+natoms]
        box_line = lines[2+natoms] if len(lines) >= 3+natoms else ""

        tname = (target_resname[:5]).ljust(5)  # GRO resname is width 5

        fixed = []
        for L in atom_lines:
            # GRO atom line canonical layout:
            # %5d %-5s %5s %5d %8.3f %8.3f %8.3f (plus optional v fields)
            # We’ll parse minimally and rebuild the left part with fixed widths,
            # then append any trailing coords/velocities as-is.
            if len(L) < 20:
                fixed.append(L)  # leave weird lines untouched
                continue
            resid   = L[0:5]
            resname = L[5:10]
            atom    = L[10:15]
            atomnr  = L[15:20]
            rest    = L[20:]  # includes coords (and velocities if present)

            # normalize resid and atomnr spacing
            resid_s  = f"{int(resid):5d}" if resid.strip().isdigit() else resid
            atomnr_s = f"{int(atomnr):5d}" if atomnr.strip().isdigit() else atomnr

            newL = f"{resid_s}{tname}{atom}{atomnr_s}{rest}"
            fixed.append(newL)

        out = [title, f"{natoms}"] + fixed + [box_line]
        gro_path.write_text("\n".join(out) + "\n")

    def normalize_resnames_to_itp(solute_gro: str | Path,
                                solvent_gro: str | Path,
                                solute_itp_out: str | Path,
                                solvent_itp_out: str | Path) -> tuple[str, str]:
        """
        Make .gro residue names match the [ moleculetype ] names in the given .itp files.
        Returns (solute_name, solvent_name).
        """
        solute_gro   = Path(solute_gro)
        solvent_gro  = Path(solvent_gro)
        solute_itp   = Path(solute_itp_out)
        solvent_itp  = Path(solvent_itp_out)

        solute_name  = _read_moleculetype_name(solute_itp)
        solvent_name = _read_moleculetype_name(solvent_itp)

        _rewrite_gro_resname(solute_gro, solute_name)
        _rewrite_gro_resname(solvent_gro, solvent_name)

        return solute_name, solvent_name

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
    outname: str = "solvated_polymer.gro",
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
    solvent_outname: str = "solvent.itp",
    solute_outname: str = "solute.itp",
):
        """
        Build toppar/forcefield.itp with [ defaults ] (once) and two [ atomtypes ] blocks
        (solute 'p' then solvent 's'). Rewrite per-molecule ITPs to remove both [ atomtypes ]
        and [ defaults ], and suffix opls_#### types in their [ atoms ] sections.
        Ensures #include lines in topol.top.

        Returns an object with .files and .log.
        """
        import re
        from pathlib import Path
        from types import SimpleNamespace

        # --- setup paths ---
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

        # --- regex helpers ---
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
            # only rename tokens that look exactly like opls_###
            return re.sub(r'^(opls_\d+)$', rf'\1{suffix}', token)

        def rewrite_atoms_section(lines: list[str], suffix: str) -> list[str]:
            """Suffix TYPE (2nd field) if it matches opls_###."""
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
                if len(toks) >= 2:
                    toks[1] = rename_opls_token(toks[1], suffix)
                    new_line = " ".join(toks).rstrip()
                    if comment:
                        new_line += " " + comment
                    out.append(new_line)
                else:
                    out.append(ln)
            return out

        def collect_atomtypes(itp_text: str, suffix: str) -> list[str]:
            """Get [ atomtypes ] lines, suffix name if opls_###, keep comments, sort by opls number."""
            sections = split_sections(itp_text)
            lines = sections.get('atomtypes', [])
            out = []
            for ln in lines:
                if not is_data_line(ln):
                    continue
                cm = COMMENT_RE.search(ln)
                comment = cm.group(0) if cm else ""
                core = ln if not cm else ln[:cm.start()]
                toks = core.split()
                if len(toks) < 7:
                    continue
                toks[0] = rename_opls_token(toks[0], suffix)
                newline = "  " + " ".join(toks)
                if comment:
                    newline += " " + comment
                out.append(newline.strip())

            def opls_num(line: str) -> int:
                m = re.search(r'opls_(\d+)', line)
                return int(m.group(1)) if m else 0

            out.sort(key=opls_num)
            return out

        def collect_defaults(itp_text: str) -> list[str]:
            """Return the [ defaults ] data+comment lines (no blank-only lines)."""
            sections = split_sections(itp_text)
            lines = sections.get('defaults', [])
            out = []
            for ln in lines:
                if ln.strip() == "":
                    continue
                out.append(ln.rstrip())
            return out

        def rewrite_itp_without_ff_sections(itp_text: str, suffix: str) -> str:
            """
            Remove [ atomtypes ] and [ defaults ] sections entirely.
            Rewrite [ atoms ] types with suffix. Leave others unchanged.
            """
            sections = split_sections(itp_text)
            out_lines: list[str] = []
            for sec_name, sec_lines in sections.items():
                if sec_name in ('atomtypes', 'defaults'):
                    continue  # drop; consolidated into forcefield.itp
                out_lines.append(f"[ {sec_name} ]")
                if sec_name == 'atoms':
                    out_lines.extend(rewrite_atoms_section(sec_lines, suffix))
                else:
                    out_lines.extend(sec_lines)
                out_lines.append("")  # blank after section
            return "\n".join(out_lines).rstrip() + "\n"

        def ensure_includes(top_path: Path, inc_paths: list[Path]) -> None:
            """Ensure #include lines (paths relative to topol.top folder) exist."""
            txt = top_path.read_text()
            rels = [str(p.relative_to(top_path.parent)) for p in inc_paths]
            missing = [r for r in rels if r not in txt]
            if not missing:
                return
            lines = txt.splitlines()
            insert_idx = 0
            for i, L in enumerate(lines):
                if L.strip().startswith("#include"):
                    insert_idx = i + 1
            new_lines = lines[:insert_idx] + [f'#include "{r}"' for r in missing] + lines[insert_idx:]
            top_path.write_text("\n".join(new_lines) + ("\n" if not txt.endswith("\n") else ""))

        # --- read inputs ---
        solvent_txt = solvent_itp.read_text()
        solute_txt  = solute_itp.read_text()

        # --- gather [ defaults ] once for forcefield.itp ---
        defaults_solute  = collect_defaults(solute_txt)
        defaults_solvent = collect_defaults(solvent_txt)
        if defaults_solute and defaults_solvent and defaults_solute != defaults_solvent:
            defaults_block = defaults_solute  # policy: prefer solute if they differ
        elif defaults_solute:
            defaults_block = defaults_solute
        elif defaults_solvent:
            defaults_block = defaults_solvent
        else:
            defaults_block = [
                "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ",
                "  1       3          yes        0.5       0.5",
            ]

        # --- gather atomtypes for both, with suffixes ---
        atomtypes_s = collect_atomtypes(solvent_txt, 's')
        atomtypes_p = collect_atomtypes(solute_txt,  'p')

        # --- write toppar/forcefield.itp ---
        ff_path = outdir / "forcefield.itp"
        with ff_path.open("w") as fh:
            fh.write("[ defaults ]\n")
            for l in defaults_block:
                fh.write(l.rstrip() + "\n")
            fh.write("\n[ atomtypes ]\n")
            for l in atomtypes_p:
                fh.write(l + "\n")
            fh.write("\n[ atomtypes ]\n")
            for l in atomtypes_s:
                fh.write(l + "\n")

        # --- write cleaned per-molecule ITPs (no defaults, no atomtypes) ---
        solvent_out = outdir / solvent_outname
        solute_out  = outdir / solute_outname
        solvent_out.write_text(rewrite_itp_without_ff_sections(solvent_txt, 's'))
        solute_out.write_text(rewrite_itp_without_ff_sections(solute_txt,  'p'))

        # --- ensure #includes in topol.top ---
        ensure_includes(topol_top, [ff_path, solvent_out, solute_out])

        # --- return small object with .files + .log ---
        files_map = {
            "forcefield_itp": str(ff_path),
            "solvent_itp": str(solvent_out),
            "solute_itp": str(solute_out),
            "topol_top": str(topol_top),
        }
        msg = (
            "Topology prepared:\n"
            f" - {ff_path}\n"
            f" - {solvent_out}\n"
            f" - {solute_out}\n"
            f"Including lines ensured in {topol_top}."
        )
        return SimpleNamespace(files=files_map, log=msg)

    def _read_moleculetype_name(itp_path: Path) -> str:
        txt = Path(itp_path).read_text()
        m = re.search(r'^\s*\[\s*moleculetype\s*\]\s*(?:;.*)?$([\s\S]*?)(?=^\s*\[|\Z)', txt, re.MULTILINE)
        for line in m.group(1).splitlines():
            s = line.strip()
            if s and not s.startswith((';', '#')):
                return s.split()[0]
        raise ValueError(f"Bad [moleculetype] in {itp_path}")

    def _rewrite_gro_resname(gro_path: Path, target_resname: str, only_if_name_in: set[str] | None = None) -> int:
        p = Path(gro_path)
        lines = p.read_text().splitlines()
        title, natoms = lines[0], int(lines[1].strip())
        atom_lines = lines[2:2+natoms]
        box_line   = lines[2+natoms]
        tname = (target_resname[:5]).ljust(5)
        changed, fixed = 0, []
        for L in atom_lines:
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
   
    def normalize_resnames_step(
        self,
        *,
        gro: str | Path,
        solute_itp: str | Path,
        solvent_itp: str | Path,
    ):
        """
        Normalize the combined GRO’s residue names *in place* to the moleculetype names
        from the ITPs. Changes only resnames that are generic:
        MOL  -> <solute moleculetype>
        SOL/WAT/HOH -> <solvent moleculetype>
        Returns an object with .files['gro'] so orchestrate can promote it.
        """
        import re
        from pathlib import Path
        from types import SimpleNamespace

        gro = Path(gro).resolve()
        solute_itp = Path(solute_itp).resolve()
        solvent_itp = Path(solvent_itp).resolve()

        def read_mtype(p: Path) -> str:
            txt = p.read_text()
            m = re.search(r'^\s*\[\s*moleculetype\s*\]\s*(?:;.*)?$([\s\S]*?)(?=^\s*\[|\Z)', txt, re.MULTILINE)
            if not m:
                raise ValueError(f"[moleculetype] not found in {p}")
            for line in m.group(1).splitlines():
                s = line.strip()
                if s and not s.startswith((';', '#')):
                    return s.split()[0]
            raise ValueError(f"Empty [moleculetype] in {p}")

        def rewrite_gro_resname(gro_path: Path, target_resname: str, only_if: set[str] | None = None) -> int:
            lines = gro_path.read_text().splitlines()
            if len(lines) < 3:
                raise ValueError(f"Invalid .gro: {gro_path}")
            title, natoms = lines[0], int(lines[1].strip())
            atom_lines = lines[2:2+natoms]
            box_line   = lines[2+natoms] if len(lines) >= 3+natoms else ""
            tname = (target_resname[:5]).ljust(5)

            changed = 0
            fixed = []
            for L in atom_lines:
                if len(L) < 20:
                    fixed.append(L); continue
                resid, resname, atom, atomnr, rest = L[0:5], L[5:10], L[10:15], L[15:20], L[20:]
                cur = resname.strip()
                if (only_if is None) or (cur in only_if):
                    resid_s  = f"{int(resid):5d}" if resid.strip().isdigit() else resid
                    atomnr_s = f"{int(atomnr):5d}" if atomnr.strip().isdigit() else atomnr
                    fixed.append(f"{resid_s}{tname}{atom}{atomnr_s}{rest}")
                    changed += 1
                else:
                    fixed.append(L)

            gro_path.write_text("\n".join([title, f"{natoms}", *fixed, box_line]) + "\n")
            return changed

        solute_name  = read_mtype(solute_itp)
        solvent_name = read_mtype(solvent_itp)
        n1 = rewrite_gro_resname(gro, solute_name, only_if={"MOL"})
        n2 = rewrite_gro_resname(gro, solvent_name, only_if={"SOL", "WAT", "HOH"})
        msg = f"Normalized {gro.name}: MOL→{solute_name} ({n1}), SOL/WAT→{solvent_name} ({n2})"

        return SimpleNamespace(files={"gro": str(gro)}, log=msg)

    def normalize_atomnames_step(
        self,
        *,
        gro: str | Path,
        solute_itp: str | Path,
    ):
        """
        Rename the atom-name column in the GRO for the *solute* residue so it matches the
        atom names from the solute ITP [ atoms ] section (5th column). Operates in place.
        Assumes a single solute residue (count=1).
        """
        import re
        from pathlib import Path
        from types import SimpleNamespace

        gro_path = Path(gro).resolve()
        itp_path = Path(solute_itp).resolve()

        if not gro_path.exists() or not itp_path.exists():
            msg = f"normalize_atomnames: missing file(s) gro={gro_path.exists()} itp={itp_path.exists()}"
            print(f"[WARN] {msg}")
            return SimpleNamespace(files={"gro": str(gro_path)}, log=msg)

        txt = itp_path.read_text()

        # moleculetype name (also your solute resname)
        mtype_m = re.search(r'^\s*\[\s*moleculetype\s*\]\s*(?:;.*)?$([\s\S]*?)(?=^\s*\[|\Z)', txt, re.MULTILINE)
        if not mtype_m:
            msg = f"normalize_atomnames: [moleculetype] not found in {itp_path.name}"
            print(f"[WARN] {msg}")
            return SimpleNamespace(files={"gro": str(gro_path)}, log=msg)
        solute_resname = None
        for line in mtype_m.group(1).splitlines():
            s = line.strip()
            if s and not s.startswith((';', '#')):
                solute_resname = s.split()[0]
                break
        if not solute_resname:
            msg = f"normalize_atomnames: empty [moleculetype] in {itp_path.name}"
            print(f"[WARN] {msg}")
            return SimpleNamespace(files={"gro": str(gro_path)}, log=msg)

        # atom names from [ atoms ] (5th column)
        atoms_m = re.search(r'^\s*\[\s*atoms\s*\]\s*(?:;.*)?$([\s\S]*?)(?=^\s*\[|\Z)', txt, re.MULTILINE)
        if not atoms_m:
            msg = f"normalize_atomnames: [ atoms ] not found in {itp_path.name}"
            print(f"[WARN] {msg}")
            return SimpleNamespace(files={"gro": str(gro_path)}, log=msg)

        atom_names = []
        for line in atoms_m.group(1).splitlines():
            s = line.strip()
            if not s or s.startswith((';', '#')):
                continue
            toks = s.split()
            # [ atoms ]: nr type resnr resid atom cgnr charge mass ...
            if len(toks) >= 5:
                atom_names.append(toks[4])

        if not atom_names:
            msg = f"normalize_atomnames: no atom names parsed from [ atoms ] in {itp_path.name}"
            print(f"[WARN] {msg}")
            return SimpleNamespace(files={"gro": str(gro_path)}, log=msg)

        # rewrite the GRO atom-name field for the first residue with resname==solute_resname
        lines = gro_path.read_text().splitlines()
        if len(lines) < 3:
            msg = f"normalize_atomnames: invalid GRO {gro_path.name}"
            print(f"[WARN] {msg}")
            return SimpleNamespace(files={"gro": str(gro_path)}, log=msg)

        title = lines[0]
        natoms = int(lines[1].strip())
        atom_lines = lines[2:2+natoms]
        box_line = lines[2+natoms] if len(lines) >= 3+natoms else ""

        # Find indices of atoms in the first solute residue block
        first_idx = None
        target_resid = None
        idxs = []
        for i, L in enumerate(atom_lines):
            if len(L) < 20:
                continue
            resid, resname = L[0:5], L[5:10]
            rname = resname.strip()
            if rname == solute_resname:
                rid = resid.strip()
                if first_idx is None:
                    first_idx = i
                    target_resid = rid  # stick to the first residue number we see
                if resid.strip() == target_resid:
                    idxs.append(i)
                else:
                    # we reached another residue, stop collecting
                    break

        if not idxs:
            msg = f"normalize_atomnames: no atoms with resname '{solute_resname}' found in {gro_path.name}"
            print(f"[WARN] {msg}")
            return SimpleNamespace(files={"gro": str(gro_path)}, log=msg)

        if len(idxs) != len(atom_names):
            msg = (f"normalize_atomnames: solute atom count mismatch "
                f"(GRO {len(idxs)} vs ITP {len(atom_names)}); skipping rename.")
            print(f"[WARN] {msg}")
            return SimpleNamespace(files={"gro": str(gro_path)}, log=msg)

        # Apply renaming (column 10–15 in GRO)
        fixed = atom_lines[:]  # copy
        for k, i in enumerate(idxs):
            L = atom_lines[i]
            resid, resname, atom, atomnr, rest = L[0:5], L[5:10], L[10:15], L[15:20], L[20:]
            new_atom = (atom_names[k][:5]).rjust(5)  # atom field is width 5, right-justified
            fixed[i] = f"{resid}{resname}{new_atom}{atomnr}{rest}"

        out_lines = [title, f"{natoms}", *fixed, box_line]
        gro_path.write_text("\n".join(out_lines) + "\n")

        msg = f"normalize_atomnames: renamed {len(idxs)} atoms in residue {solute_resname} from GRO to match ITP"
        print(f"[INFO] {msg}")
        return SimpleNamespace(files={"gro": str(gro_path)}, log=msg)

    def _read_moleculetype_name(itp_path: Path) -> str:
        """
        Return the first data token in [ moleculetype ] from an .itp file.
        """
        text = itp_path.read_text()
        m = re.search(r'^\s*\[\s*moleculetype\s*\]\s*(?:;.*)?$([\s\S]*?)(?=^\s*\[|\Z)',
                    text, re.MULTILINE)
        if not m:
            raise ValueError(f"[moleculetype] not found in {itp_path}")
        for line in m.group(1).splitlines():
            s = line.strip()
            if s and not s.startswith((';', '#')):
                return s.split()[0]
        raise ValueError(f"Empty [moleculetype] block in {itp_path}")

    def _count_residues_in_gro(gro_path: Path) -> Dict[str, int]:
        """
        Count residues in a .gro by resname. We count a residue whenever the (resid,resname)
        tuple changes (standard GRO layout).
        Returns dict: {resname: count}
        """
        lines = gro_path.read_text().splitlines()
        if len(lines) < 3:
            raise ValueError(f"Invalid .gro: {gro_path}")
        natoms = int(lines[1].strip())
        atom_lines = lines[2:2+natoms]

        counts: Dict[str, int] = {}
        prev_key: Optional[Tuple[str, str]] = None
        for L in atom_lines:
            if len(L) < 20:
                continue
            resid = L[0:5].strip()
            resnm = L[5:10].strip()
            key = (resid, resnm)
            if key != prev_key:
                counts[resnm] = counts.get(resnm, 0) + 1
                prev_key = key
        return counts

    def _rel_include(from_file: Path, target: Path) -> str:
        """
        Return a quoted include path relative to 'from_file' parent.
        """
        rel = target.resolve().relative_to(from_file.parent.resolve())
        return f"#include \"{rel.as_posix()}\""

    def normalize_topol_step(
    self,
    *,
    gro: str | Path,                        # combined, normalized GRO (e.g., solvated.gro)
    forcefield_itp: str | Path,             # toppar/forcefield.itp
    solvent_itp: str | Path,                # toppar/dcb.itp (or your solvent)
    solute_itp: str | Path,                 # toppar/c6.itp (or your solute)
    topol_top: str | Path = "topol.top",    # output to write (rebuilt)
    system_title: str = "Polymer in solvent in water",
    include_others: bool = False,           # set True to append any extra residue types seen in GRO
):
        """
        Rebuild topol.top with:
        - 3 includes (forcefield, solvent, solute) using paths relative to topol.top
        - [ system ] title
        - [ molecules ] with counts taken from GRO resnames (must already be normalized)
        """
        import re
        from pathlib import Path
        from types import SimpleNamespace

        gro = Path(gro).resolve()
        forcefield_itp = Path(forcefield_itp).resolve()
        solvent_itp    = Path(solvent_itp).resolve()
        solute_itp     = Path(solute_itp).resolve()
        topol_top      = Path(topol_top).resolve()

        if not gro.exists():
            raise FileNotFoundError(f"GRO not found: {gro}")
        for p in (forcefield_itp, solvent_itp, solute_itp):
            if not p.exists():
                raise FileNotFoundError(f"ITP not found: {p}")

        def read_moleculetype_name(itp_path: Path) -> str:
            txt = itp_path.read_text()
            m = re.search(r'^\s*\[\s*moleculetype\s*\]\s*(?:;.*)?$([\s\S]*?)(?=^\s*\[|\Z)', txt, re.MULTILINE)
            if not m:
                raise ValueError(f"[moleculetype] not found in {itp_path}")
            for line in m.group(1).splitlines():
                s = line.strip()
                if s and not s.startswith((';', '#')):
                    return s.split()[0]
            raise ValueError(f"Empty [moleculetype] block in {itp_path}")

        def count_residues_in_gro(gro_path: Path) -> dict[str, int]:
            lines = gro_path.read_text().splitlines()
            if len(lines) < 3:
                raise ValueError(f"Invalid .gro: {gro_path}")
            natoms = int(lines[1].strip())
            atom_lines = lines[2:2+natoms]
            counts: dict[str, int] = {}
            prev_key = None
            for L in atom_lines:
                if len(L) < 20:
                    continue
                resid = L[0:5].strip()
                resnm = L[5:10].strip()
                key = (resid, resnm)
                if key != prev_key:
                    counts[resnm] = counts.get(resnm, 0) + 1
                    prev_key = key
            return counts

        def rel_include(from_file: Path, target: Path) -> str:
            rel = target.resolve().relative_to(from_file.parent.resolve())
            return f'#include "{rel.as_posix()}"'

        # names & counts
        solute_name  = read_moleculetype_name(solute_itp)
        solvent_name = read_moleculetype_name(solvent_itp)
        counts       = count_residues_in_gro(gro)
        solute_count  = counts.get(solute_name, 0)
        solvent_count = counts.get(solvent_name, 0)

        # includes relative to topol.top
        inc_ff  = rel_include(topol_top, forcefield_itp)
        inc_sol = rel_include(topol_top, solvent_itp)
        inc_solute = rel_include(topol_top, solute_itp)

        mol_lines = [
            "; name  number",
            f"{solute_name:<16} {solute_count}",
            f"{solvent_name:<16} {solvent_count}",
        ]
        if include_others:
            for rn, n in sorted(counts.items()):
                if rn not in (solute_name, solvent_name):
                    mol_lines.append(f"{rn:<16} {n}")

        out_lines = [
            "; auto-generated skeleton",
            inc_ff,
            inc_sol,
            inc_solute,
            "",
            "[ system ]",
            system_title,
            "",
            "[ molecules ]",
            *mol_lines,
            "",
        ]
        topol_top.parent.mkdir(parents=True, exist_ok=True)
        topol_top.write_text("\n".join(out_lines))

        others = {rn: n for rn, n in counts.items() if rn not in (solute_name, solvent_name)}
        msg = (f"normalize_topol: wrote {topol_top.name} with "
        f"{solute_name}={solute_count}, {solvent_name}={solvent_count}; others={others}")

        print("[INFO] " + msg)
        return SimpleNamespace(files={"topol_top": str(topol_top), "gro": str(gro)}, log=msg)

    # ── orchestrator ───────────────────────────────────────────────────────────

    def orchestrate(self, steps, *, inputs, overrides=None, strict=True) -> Dict[str, StepOutput]:
        """
        Run the registered steps in order, passing only the parameters each step accepts,
        and promote key outputs into the shared ctx for downstream steps.
        No GRO/ITP normalization logic lives here.
        """
        import inspect
        from typing import Dict

        ctx = dict(inputs)
        results: Dict[str, StepOutput] = {}
        overrides = overrides or {}

        for name in steps:
            fn = self._registry.get(name)
            if fn is None:
                if strict:
                    raise KeyError(f"Unknown step: {name}")
                print(f"[WARN] Skipping unknown step: {name}")
                continue

            # Merge inputs with per-step overrides, filter by fn signature
            merged = {**ctx, **overrides.get(name, {})}
            sig = inspect.signature(fn)
            allowed = {k: v for k, v in merged.items() if k in sig.parameters}

            # Run the step
            out: StepOutput = fn(**allowed)  # type: ignore[misc]
            results[name] = out

            # Promote common outputs for downstream steps
            if name == "create_box":
                ctx["solute_box_gro"] = out.files["gro"]
                ctx["gro"] = out.files["gro"]

            elif name == "solvent_box":
                ctx["solvent_box_gro"] = out.files["gro"]
                ctx["gro"] = out.files["gro"]

            elif name == "prepare_topology":
                # expose paths produced by the topology-prep step
                ctx["forcefield_itp"] = out.files.get("forcefield_itp")
                ctx["solvent_itp"]    = out.files.get("solvent_itp")
                ctx["solute_itp"]     = out.files.get("solute_itp")
                ctx["topol_top"]      = out.files.get("topol_top")
                
            elif name == "normalize_resnames":
                ctx["gro"] = out.files["gro"]    # ensure EM/NVT read the normalized GRO

            elif name == "normalize_atomnames":
                ctx["gro"] = out.files["gro"]        
            elif name in ("combine_solvate", "solvate"):
                ctx["gro"] = out.files["gro"]

            elif name in ("mdrun_em", "mdrun_nvt"):
                if out.files.get("gro"):
                    ctx["gro"] = out.files["gro"]
            elif name == "normalize_topol":
                ctx["topol_top"] = out.files["topol_top"]
                # ctx["gro"] already set earlier; keep as-is
        return results


