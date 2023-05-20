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

# This test ensures that the indexing of the PMG structure is the same of the ASE Atoms object.
# This convenience is used for BXB bond angle determination, and corner/edge/face sharing detection.
def test_pmg_ase_indexing():
    perov = pk(test_atoms)
    pmg_sites = np.zeros((len(perov.atoms), 3))
    for i in range(len(perov.atoms)):
        pmg_sites[i][0] = perov.structure[i].x
        pmg_sites[i][1] = perov.structure[i].y
        pmg_sites[i][2] = perov.structure[i].z
    assert np.allclose(perov.atoms.positions, pmg_sites)
    
