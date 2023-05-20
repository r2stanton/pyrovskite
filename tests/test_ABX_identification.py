from pyrovskite.perovskite import Perovskite as pk
import numpy as np
import os, ase.io

try:
    dirr = os.path.join(os.getcwd(), "data")
    file_name = os.path.join(dirr, "test_perovskite_1.cif")
    test_atoms = ase.io.read(file_name)
    test_doublepk = ase.io.read(os.path.join(dirr, "3AMP_CsInSbI_dj.cif"))
except:
    dirr = os.path.join(os.getcwd(), "tests/data")
    file_name = os.path.join(dirr, "test_perovskite_1.cif")
    test_doublepk = ase.io.read(os.path.join(dirr, "3AMP_CsInSbI_dj.cif"))
    test_atoms = ase.io.read(file_name)



def test_X_detection():
    perov = pk(test_atoms)
    assert perov.X == "I" , "X-anion detection went wrong."

def test_B_detection():
    perov = pk(test_atoms)
    assert perov.B == "Pb" , "B-cation detection went wrong."

def test_double_perovskite_detection():
    perov = pk(test_doublepk)
    b1 = perov.B == "In" or perov.B == "Sb"
    b2 = perov.Bp == "In" or perov.Bp == "Sb"
    b3 = perov.Bp == perov.B

    assert b1 and b2 and not b3 , "B, Bp detection in double perovskite is not correct."
    assert perov.X == "I" , "X-anion detection in double perovskite is not correct."

