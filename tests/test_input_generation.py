from pyrovskite.io import xTB_input, GPAW_input
import ase.io
import os

try:
    dirr = os.path.join(os.getcwd(), "data")
    file_name = os.path.join(dirr, "3AMP_CsInSbI_dj.cif")
    test_atoms = ase.io.read(file_name)
except:
    dirr = os.path.join(os.getcwd(), "tests/data")
    file_name = os.path.join(dirr, "3AMP_CsInSbI_dj.cif")
    test_atoms = ase.io.read(file_name)


# Holding off on this test while whether to require the paths for the params
# as environment variables is sorted.

# def test_vcopt_xtb():
    # atoms = test_atoms.copy()
    # xTB_input(calc_type = "vcopt", prefix = dirr+"/testvcopt", atoms = atoms,
            # xtb_path = "placeholder", d3_path = "placeholder2")


def test_fcopt_gpaw():
    atoms = test_atoms.copy()
    GPAW_input(calc_type = "fcopt", prefix = dirr+"/testvcopt", atoms = atoms)
