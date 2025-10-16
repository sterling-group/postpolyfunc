from __future__ import annotations
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# -------------------------------
# Helpers
# -------------------------------

def _find_executable(candidates: Tuple[str, ...] = ("ligpargen", "LigParGen", "LigParGen2.1")) -> str:
    """Return first found executable name or raise FileNotFoundError."""
    for name in candidates:
        exe = shutil.which(name)
        if exe:
            return exe
    raise FileNotFoundError(
        "Could not find LigParGen executable in PATH. "
        "Tried: ligpargen, LigParGen, LigParGen2.1"
    )


def _bool_flag(flag: str, enabled: bool) -> List[str]:
    return [flag] if enabled else []


# -------------------------------
# Dataclass config (mirrors CLI)
# -------------------------------

@dataclass
class LigParGenConfig:
    # Required-ish
    path: Path                              # -p / --path (working dir; outputs go here)
    resname: Optional[str] = None           # -r / --resname (3-letter)
    molname: Optional[str] = None           # -n / --molname

    # Choose ONE: smiles or input file
    smile: Optional[str] = None             # -s / --smile
    ifile: Optional[Path] = None            # -i / --ifile

    # Charge/gen options
    cgen: str = "CM1A-LBCC"                 # -cgen {CM1A,CM1A-LBCC}
    charge: int = 0                         # -c  -10..10
    opt: int = 0                            # -o  {0,1,2,3}

    # Flags
    debug: bool = False                     # --debug
    verbose: bool = False                   # --verbose

    # Optional Alchemical B counterpart
    smileB: Optional[str] = None            # -sb / --smileB
    ifileB: Optional[Path] = None           # -ib / --ifileB
    cgenB: Optional[str] = None             # -cgenB {CM1A,CM1A-LBCC}
    chargeB: Optional[int] = None           # -cb
    optB: Optional[int] = None              # -ob

    # Consistency checker
    checker: bool = False                   # -check / --checker

    # Execution
    executable: Optional[str] = None        # override binary name; otherwise auto-discovered
    extra_args: List[str] = field(default_factory=list)
    env: Optional[Dict[str, str]] = None
    timeout_s: Optional[int] = None         # subprocess timeout

    def validate(self) -> None:
        if not (self.smile or self.ifile):
            raise ValueError("Provide either 'smile' or 'ifile' to LigParGenConfig.")
        if self.smile and self.ifile:
            raise ValueError("Provide only one of 'smile' or 'ifile', not both.")
        if self.smileB and self.ifileB:
            raise ValueError("Provide only one of 'smileB' or 'ifileB', not both.")
        if self.cgen not in ("CM1A", "CM1A-LBCC"):
            raise ValueError("cgen must be 'CM1A' or 'CM1A-LBCC'.")
        if self.cgenB and self.cgenB not in ("CM1A", "CM1A-LBCC"):
            raise ValueError("cgenB must be 'CM1A' or 'CM1A-LBCC'.")
        if not self.path:
            raise ValueError("'path' (working directory) is required.")
        if self.opt not in (0, 1, 2, 3):
            raise ValueError("opt must be one of {0,1,2,3}.")
        if self.optB is not None and self.optB not in (0, 1, 2, 3):
            raise ValueError("optB must be one of {0,1,2,3}.")


# -------------------------------
# Command builder + runner
# -------------------------------

def build_ligpargen_cmd(cfg: LigParGenConfig) -> List[str]:
    cfg.validate()
    exe = cfg.executable or _find_executable()
    cmd: List[str] = [exe]

    # Working directory
    cmd += ["-p", str(cfg.path)]

    # Names
    if cfg.resname:
        cmd += ["-r", cfg.resname]
    if cfg.molname:
        cmd += ["-n", cfg.molname]

    # Primary input
    if cfg.smile:
        cmd += ["-s", cfg.smile]
    elif cfg.ifile:
        cmd += ["-i", str(cfg.ifile)]

    # Charge/gen
    cmd += ["-cgen", cfg.cgen, "-c", str(cfg.charge), "-o", str(cfg.opt)]

    # Flags
    cmd += _bool_flag("--debug", cfg.debug)
    cmd += _bool_flag("--verbose", cfg.verbose)

    # Alchemical inputs (optional)
    if cfg.smileB:
        cmd += ["-sb", cfg.smileB]
    if cfg.ifileB:
        cmd += ["-ib", str(cfg.ifileB)]
    if cfg.cgenB:
        cmd += ["-cgenB", cfg.cgenB]
    if cfg.chargeB is not None:
        cmd += ["-cb", str(cfg.chargeB)]
    if cfg.optB is not None:
        cmd += ["-ob", str(cfg.optB)]

    # Checker
    if cfg.checker:
        cmd += ["-check"]

    # Extra passthrough
    cmd += cfg.extra_args

    return cmd


def run_ligpargen(cfg: LigParGenConfig) -> subprocess.CompletedProcess:
    """
    Execute LigParGen with the given config.
    Returns the CompletedProcess. Raises on non-zero exit.
    """
    cfg.path.mkdir(parents=True, exist_ok=True)
    cmd = build_ligpargen_cmd(cfg)

    # Merge custom env over current env if provided
    env = os.environ.copy()
    if cfg.env:
        env.update(cfg.env)

    proc = subprocess.run(
        cmd,
        cwd=str(cfg.path),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=cfg.timeout_s,
        check=False,
    )

    if proc.returncode != 0:
        raise RuntimeError(
            "LigParGen failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return code: {proc.returncode}\n"
            f"STDOUT:\n{proc.stdout}\n"
            f"STDERR:\n{proc.stderr}"
        )
    return proc


# -------------------------------
# Output collector
# -------------------------------

def collect_outputs(workdir: Path) -> Dict[str, List[Path]]:
    """
    Inspect workdir and collect common LigParGen artifacts by family.
    Returns a dict e.g.:
      {
        'gromacs': [*.top, *.itp, *.gro, *.itp.gz, *.ndx, ...],
        'openmm':  [*.xml],
        'charmm':  [*.prm, *.rtf, ...],
        'boss':    [*.mol, *.z, Q/*],
        'pdb2pqr': [*.pqr],
        'general': [*.pdb, *.mol2, *.sdf, *.log, ...]
      }
    """
    workdir = Path(workdir)
    out: Dict[str, List[Path]] = {k: [] for k in ["gromacs", "openmm", "charmm", "boss", "pdb2pqr", "general"]}

    # GROMACS
    for pat in ("*.top", "*.itp", "*.gro", "*.ndx"):
        out["gromacs"] += list(workdir.glob(pat))

    # OpenMM
    out["openmm"] += list(workdir.glob("*.xml"))

    # CHARMM/NAMD
    for pat in ("*.prm", "*.rtf", "*.psf"):
        out["charmm"] += list(workdir.glob(pat))

    # BOSS/MCPRO/Q
    out["boss"] += list(workdir.glob("*.mol"))
    out["boss"] += list(workdir.glob("*.z"))
    qdir = workdir / "Q"
    if qdir.is_dir():
        out["boss"] += list(qdir.glob("*"))

    # PDB2PQR
    out["pdb2pqr"] += list(workdir.glob("*.pqr"))

    # General / coordinates / logs
    for pat in ("*.pdb", "*.mol2", "*.sdf", "*.log", "*.out", "*.in", "*.dat", "*.csv"):
        out["general"] += list(workdir.glob(pat))

    # Filter empties
    return {k: v for k, v in out.items() if v}


# -------------------------------
# High-level convenience
# -------------------------------

def generate_parameters(
    *,
    workdir: Path,
    resname: str,
    molname: str,
    smile: str | None = None,
    ifile: Path | None = None,
    charge: int = 0,
    cgen: str = "CM1A-LBCC",
    opt: int = 0,
    debug: bool = False,
    verbose: bool = False,
    smilesB: str | None = None,
    ifileB: Path | None = None,
    cgenB: str | None = None,
    chargeB: int | None = None,
    optB: int | None = None,
    checker: bool = False,
    executable: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
    timeout_s: Optional[int] = None,
) -> Dict[str, List[Path]]:
    """
    One-shot: run LigParGen and return a dict of discovered output files by tool family.
    """
    cfg = LigParGenConfig(
        path=workdir,
        resname=resname,
        molname=molname,
        smile=smile,
        ifile=ifile,
        cgen=cgen,
        charge=charge,
        opt=opt,
        debug=debug,
        verbose=verbose,
        smileB=smilesB,
        ifileB=ifileB,
        cgenB=cgenB,
        chargeB=chargeB,
        optB=optB,
        checker=checker,
        executable=executable,
        extra_args=extra_args or [],
        timeout_s=timeout_s,
    )
    run_ligpargen(cfg)
    return collect_outputs(workdir)
