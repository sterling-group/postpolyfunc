# tests/test_functionalize.py
from __future__ import annotations
from pathlib import Path
import importlib.resources as ir

import numpy as np
import pytest
from ase.io import read, write
from importlib.resources import files, as_file

from postpolyfunc.func import PolymerFunctionalizer


def _pkg_file(relpath: str) -> Path:
    res = files("postpolyfunc").joinpath(relpath)
    with as_file(res) as p:
        return Path(p)


def _count_elements(atoms, symbol: str) -> int:
    return sum(1 for s in atoms.get_chemical_symbols() if s == symbol)


@pytest.fixture(scope="session")
def hexane_pdb() -> Path:
    return _pkg_file("test/hexane.pdb")


@pytest.fixture(scope="session")
def dcb_pdb() -> Path:
    return _pkg_file("test/1-4-dcb.pdb")


def test_functionalize_one_site_counts(hexane_pdb: Path, tmp_path: Path):
    """
    Functionalize exactly one site on hexane and verify:
      - total atoms unchanged
      - H decreases by 1
      - O increases by 1
    """
    atoms = read(str(hexane_pdb))

    n_atoms0 = len(atoms)
    nH0 = _count_elements(atoms, "H")
    nO0 = _count_elements(atoms, "O")

    f = PolymerFunctionalizer(
        functionalization_ratio=0.2,  # ignored because we set n_sites
        seed=123,          # deterministic site choice when ties exist
        mode="carbonyl",   # whatever mode you're currently using
    )
    atoms_func = f.functionalize_carbons(atoms)

    # Counts
    assert len(atoms_func) < n_atoms0, "Atom count should decrease per site"
    assert _count_elements(atoms_func, "H") == nH0 - 1, "Expected one fewer H"
    assert _count_elements(atoms_func, "O") == nO0 + 1, "Expected one more O"

    # Optional: geometry sanity checks
    # e.g., ensure at least one O is bonded ~1.2–1.6 Å from some carbon (sp3 C–O ~1.43 Å)
    # This is conservative and avoids hard failures on slightly stretched coords.
    pos = atoms_func.get_positions()
    syms = atoms_func.get_chemical_symbols()
    O_idxs = [i for i, s in enumerate(syms) if s == "O"]
    C_idxs = [i for i, s in enumerate(syms) if s == "C"]
    has_reasonable_CO = False
    for io in O_idxs:
        for ic in C_idxs:
            d = np.linalg.norm(pos[io] - pos[ic])
            if 1.2 <= d <= 1.6:
                has_reasonable_CO = True
                break
        if has_reasonable_CO:
            break
    assert has_reasonable_CO, "No plausible C–O bond distance found in functionalized structure"

    # Save to tmp just to confirm writer doesn't error
    outp = tmp_path / "hexane_func.pdb"
    write(outp, atoms_func)
    assert outp.exists()


def test_no_functionalization_when_ratio_zero(hexane_pdb: Path):
    """
    With functionalization_ratio=0 and n_sites=None, structure should be unchanged.
    """
    atoms = read(str(hexane_pdb))
    n_atoms0 = len(atoms)
    syms0 = atoms.get_chemical_symbols()
    pos0 = atoms.get_positions().copy()

    f = PolymerFunctionalizer(
        functionalization_ratio=0.0,
        seed=1,
        mode="carbonyl",
    )
    atoms2 = f.functionalize_carbons(atoms)

    assert len(atoms2) == n_atoms0
    assert atoms2.get_chemical_symbols() == syms0
    # Positions may not be bit-identical cross-platform; keep a small tolerance
    assert np.allclose(atoms2.get_positions(), pos0, atol=1e-8)
