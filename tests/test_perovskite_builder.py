from pyrovskite.perovskite import Perovskite as pk
from pyrovskite.builder import make_2drp, make_2d_double, make_monolayer, make_dj, make_bulk, make_double
import ase.io, os, numpy as np

# make pytest work from either base directory, or tests directory.
if os.path.isdir("data"):
    ddir = os.path.join(os.getcwd(), "data")
else:
    ddir = os.path.join(os.path.join(os.getcwd(), "tests"), "data")

Ap = ase.io.read(os.path.join(ddir, 'benzylammonium.xyz'))
A = ase.io.read(os.path.join(ddir, 'methylammonium.xyz'))


def test_monolayer():
    new_atoms = make_monolayer(Ap, 'Cs', 'Pb', 'I', 2, 3.1, exp=True)
    comp_atoms = ase.io.read(os.path.join(ddir, "test_mono.json"))
    assert comp_atoms.get_chemical_symbols() == comp_atoms.get_chemical_symbols() , "Chemical symbols don't match test data"
    assert np.allclose(comp_atoms.positions, new_atoms.positions) , "Atomic positions don't match test data."

def test_2drp():
    new_atoms = make_2drp(Ap, 'Cs', 'Pb', 'I', 3, 3.1)
    comp_atoms = ase.io.read(os.path.join(ddir, "test_rp.json"))

    assert comp_atoms.get_chemical_symbols() == comp_atoms.get_chemical_symbols() , "Chemical symbols don't match test data"
    assert np.allclose(comp_atoms.positions, new_atoms.positions) , "Atomic positions don't match test data."

def test_dj():
    new_atoms = make_dj(Ap, 'Cs', 'Pb', 'I', 3, 3.1)
    comp_atoms = ase.io.read(os.path.join(ddir, "test_dj.json"))

    assert comp_atoms.get_chemical_symbols() == comp_atoms.get_chemical_symbols() , "Chemical symbols don't match test data"
    assert np.allclose(comp_atoms.positions, new_atoms.positions) , "Atomic positions don't match test data."

def test_2d_double_rp():
    new_atoms = make_2d_double(Ap, 'Cs', 'Pb', 'Sn', 'I', 3, 3.1)
    comp_atoms = ase.io.read(os.path.join(ddir, "test_2ddouble_rp.json"))

    assert comp_atoms.get_chemical_symbols() == comp_atoms.get_chemical_symbols() , "Chemical symbols don't match test data"
    assert np.allclose(comp_atoms.positions, new_atoms.positions) , "Atomic positions don't match test data."

def test_2d_double_dj():
    new_atoms = make_2d_double(Ap, 'Cs', 'Pb', 'Sn', 'I', 3, 3.1, phase = 'dj')
    comp_atoms = ase.io.read(os.path.join(ddir, "test_2ddouble_dj.json"))

    assert comp_atoms.get_chemical_symbols() == comp_atoms.get_chemical_symbols() , "Chemical symbols don't match test data"
    assert np.allclose(comp_atoms.positions, new_atoms.positions) , "Atomic positions don't match test data."

def test_double():
    new_atoms = make_double('Cs', 'Pb', 'Sn', 'I', 3.1)
    comp_atoms = ase.io.read(os.path.join(ddir, "test_double.json"))
    assert comp_atoms.get_chemical_symbols() == comp_atoms.get_chemical_symbols() , "Chemical symbols don't match test data"
    assert np.allclose(comp_atoms.positions, new_atoms.positions) , "Atomic positions don't match test data."

def test_bulk():
    new_atoms = make_bulk(A, 'Pb', 'I', 3.1)
    comp_atoms = ase.io.read(os.path.join(ddir, "test_bulk.json"))
    assert comp_atoms.get_chemical_symbols() == comp_atoms.get_chemical_symbols() , "Chemical symbols don't match test data"
    assert np.allclose(comp_atoms.positions, new_atoms.positions) , "Atomic positions don't match test data."




