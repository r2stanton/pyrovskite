from pyrovskite.perovskite import Perovskite as pk
import numpy as np
import os, ase.io

try:
    dirr = os.path.join(os.getcwd(), "data")
    file_name = os.path.join(dirr, "test_perovskite_1.cif")
    test_atoms = ase.io.read(file_name)
except:
    dirr = os.path.join(os.getcwd(), "tests/data")
    file_name = os.path.join(dirr, "test_perovskite_1.cif")
    test_atoms = ase.io.read(file_name)

def test_load_from_file_name():
    perovskite = pk(file_name)
    assert perovskite.atoms.get_chemical_symbols() == test_atoms.get_chemical_symbols() , "Chemical Symbols don't match."
    assert np.allclose(perovskite.atoms.positions, test_atoms.positions), "Positions don't match."

def test_default_constructor():
    perovskite = pk(test_atoms.copy())
    assert perovskite.atoms.get_chemical_symbols() == test_atoms.get_chemical_symbols() , "Chemical Symbols don't match."
    assert np.allclose(perovskite.atoms.positions, test_atoms.positions), "Positions don't match."


