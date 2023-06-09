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

def test_octahedra():
    perov = pk(test_atoms.copy())
    print(perov.identify_octahedra())

    data_arr = np.array([[[ 0.,          8.333691,    8.21675   ],
                      [ 2.29018543,  6.25001254,  7.9603874 ],
                      [-2.29018543,  6.25001254,  8.4731126 ],
                      [-0.31899962,  8.32411365,  5.02240627],
                      [ 0.31899962,  8.32411365, 11.41109373],
                      [ 2.00321457, 10.89921254,  7.9603874 ],
                      [-2.00321457, 10.89921254,  8.4731126 ]],
                     [[ 4.2934    ,  5.613909  , 24.65025   ],
                      [ 2.00321457,  7.69758746, 24.3938874 ],
                      [ 6.58358543,  7.69758746, 24.9066126 ],
                      [ 3.97440038,  5.62348635, 27.84459373],
                      [ 4.61239962,  5.62348635, 21.45590627],
                      [ 2.29018543,  3.04838746, 24.3938874 ],
                      [ 6.29661457,  3.04838746, 24.9066126 ]],
                     [[ 0.        ,  0.964709  , 24.65025   ],
                      [-2.29018543,  3.04838746, 24.9066126 ],
                      [ 2.29018543,  3.04838746, 24.3938874 ],
                      [-0.31899962,  0.97428635, 21.45590627],
                      [ 0.31899962,  0.97428635, 27.84459373],
                      [-2.00321457, -1.60081254, 24.9066126 ],
                      [ 2.00321457, -1.60081254, 24.3938874 ]],
                     [[ 4.2934    ,  3.684491  ,  8.21675   ],
                      [ 2.00321457,  1.60081254,  7.9603874 ],
                      [ 6.58358543,  1.60081254,  8.4731126 ],
                      [ 4.61239962,  3.67491365,  5.02240627],
                      [ 3.97440038,  3.67491365, 11.41109373],
                      [ 2.29018543,  6.25001254,  7.9603874 ],
                      [ 6.29661457,  6.25001254,  8.4731126 ]]])

    assert np.allclose(perov.octahedra, data_arr) , "Octahedra data does not match test case"


