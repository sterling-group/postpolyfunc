import logging
import shutil
from pathlib import Path

def setup_logging(verbosity: int = 1) -> None:
    level = logging.WARNING if verbosity == 0 else logging.INFO if verbosity == 1 else logging.DEBUG
    logging.basicConfig(level=level, format="%(asctime)s [%(levelname)s] %(message)s")

def keep_only_root_gmx(outdir: Path) -> None:
    """
    In `outdir`, keep only *.gmx.gro and *.gmx.itp at the top level.
    Move all other solute.* / solvent.* files and the solute/ and solvent/ dirs
    into outdir/other_tops/.
    """
    outdir = Path(outdir)
    other = outdir / "other_tops"
    other.mkdir(parents=True, exist_ok=True)

    # Keepers at root
    keep = {"solute.gmx.gro", "solute.gmx.itp", "solvent.gmx.gro", "solvent.gmx.itp"}

    # Move files
    for item in outdir.iterdir():
        if item.is_file():
            name = item.name
            # Only consider solute.* and solvent.* at root
            if (name.startswith("solute.") or name.startswith("solvent.")) and name not in keep:
                shutil.move(str(item), other / name)

    # Move the LigParGen subfolders themselves if they exist (they contain non-GMX formats)
    for d in ("solute", "solvent"):
        p = outdir / d
        if p.exists() and p.is_dir():
            shutil.move(str(p), other / d)

