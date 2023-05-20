from pyrovskite.perovskite import Perovskite as pk
import numpy as np
import os, ase.io

try:
    dirr = os.path.join(os.getcwd(), "data")
    file_name = os.path.join(dirr, "3AMP_CsInSbI_dj.cif")
    test_atoms = ase.io.read(file_name)
except:
    dirr = os.path.join(os.getcwd(), "tests/data")
    file_name = os.path.join(dirr, "3AMP_CsInSbI_dj.cif")
    test_atoms = ase.io.read(file_name)

def test_delta():
    perov = pk(test_atoms.copy())

    # Test that all return types work.
    delta = perov.compute_delta(return_type = "delta")
    octahedra_delta = perov.compute_delta(return_type = "otahedra_delta")
    delta, octahedra_delta = perov.compute_delta(return_type = "both")

    # Compare to test data.
    assert np.abs(delta - 0.00034644998540341607) < 1e-5 , "Average delta not matching test data."
    assert np.allclose(octahedra_delta, np.array([4.63862260e-04, 2.44981753e-04, 6.07915749e-04, 6.90401800e-05])) , "octahedra_delta not matchign test data."



def test_sigma():
    perov = pk(test_atoms.copy())

    # Test that all return types work.
    sigma = perov.compute_sigma(return_type = "sigma")
    octahedra_sigma = perov.compute_sigma(return_type = "octahedra_sigma")
    sigma, octahedra_sigma = perov.compute_sigma(return_type = "both")

    # Compare to test data.
    assert np.abs(sigma - 3.1056874220863913) < 1e-5 , "Average sigma not matching test data."
    assert np.allclose(octahedra_sigma, [3.99967035, 1.68098048, 2.14586783, 4.59623103])
