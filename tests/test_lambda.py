from pyrovskite.perovskite import Perovskite as pk
import numpy as np
import os, ase.io

try:
    dirr = os.path.join(os.getcwd(), "data")
    data3d = os.path.join(dirr, "mapbi_3d_diag.cif")
    data2d = os.path.join(dirr, "mapbi_2d_diag.cif")
    data3d = ase.io.read(data3d)
    data2d = ase.io.read(data2d)
except:
    dirr = os.path.join(os.getcwd(), "tests/data")
    data3d = os.path.join(dirr, "mapbi_3d_diag.cif")
    data2d = os.path.join(dirr, "mapbi_2d_diag.cif")
    data3d = ase.io.read(data3d)
    data2d = ase.io.read(data2d)
    
def test_lambda_2d():
    d2 = data2d.copy()
    d2 = pk(d2)
    ol3, ol2, ol1, l3, l2, l1 = d2.compute_lambda(visualize=False, return_type = 'both', scaled = True, ratio = True)
    assert abs(l3-0.030802420477229) < 1e-5 , "Scaled lambda3 computation doesn't match data."
    assert abs(l2-1.801065580608829) < 1e-5 , "Scaled lambda2 computation doesn't match data."
    assert abs(l1-0.017102331424721) < 1e-5 , "Scaled lambda1 computation doesn't match data."

    ol3, ol2, ol1, l3, l2, l1 = d2.compute_lambda(visualize=False, return_type = 'both', scaled = False, ratio = True)
    assert abs(l3-0.0165290806673103) < 1e-5 , "Unscaled lambda3 computation doesn't match data."
    assert abs(l2-0.9994425621873911) < 1e-5 , "Unscaled lambda2 computation doesn't match data."
    assert abs(l1-0.0165382997409421) < 1e-5 , "Unscaled lambda1 computation doesn't match data."

    ol3, ol2, l3, l2= d2.compute_lambda(visualize=False, return_type = 'both', scaled = False, ratio = False)
    assert abs(l3-0.0165290806673103) < 1e-5 , "Unscaled lambda3 computation doesn't match data."
    assert abs(l2-0.9994425621873911) < 1e-5 , "Unscaled lambda2 computation doesn't match data."

    ol3, ol2, l3, l2 = d2.compute_lambda(visualize=False, return_type = 'both', scaled = True, ratio = False)
    assert abs(l3-0.030802420477229) < 1e-5 , "Scaled lambda3 computation doesn't match data, for ratioless option."
    assert abs(l2-1.801065580608829) < 1e-5 , "Scaled lambda2 computation doesn't match data, for ratioless option."



def test_lambda_3d():
    d3 = data3d.copy()
    d3 = pk(d3)
    ol3, ol2, ol1, l3, l2, l1 = d3.compute_lambda(visualize=False, return_type = 'both', scaled = True, ratio = True)
    assert abs(l3-1.9207281847632536) < 1e-5 , "Scaled lambda3 computation doesn't match data."
    assert abs(l2-1.2643724892960577) < 1e-5 , "Scaled lambda2 computation doesn't match data."
    assert abs(l1-1.5191157677217602) < 1e-5 , "Scaled lambda1 computation doesn't match data."

    ol3, ol2, ol1, l3, l2, l1 = d3.compute_lambda(visualize=False, return_type = 'both', scaled = False, ratio = True)
    assert abs(l3-0.8525261226466953) < 1e-5 , "Unscaled lambda3 computation doesn't match data."
    assert abs(l2-0.9743159081299232) < 1e-5 , "Unscaled lambda2 computation doesn't match data."
    assert abs(l1-0.87499969520462) < 1e-5 ,   "Unscaled lambda1 computation doesn't match data."


    ol3, ol2, l3, l2 = d3.compute_lambda(visualize=False, return_type = 'both', scaled = True, ratio = False)
    assert abs(l3-1.9207281847632536) < 1e-5 , "Scaled lambda3 computation doesn't match data."
    assert abs(l2-1.2643724892960577) < 1e-5 , "Scaled lambda2 computation doesn't match data."

    ol3, ol2, l3, l2 = d3.compute_lambda(visualize=False, return_type = 'both', scaled = False, ratio = False)
    assert abs(l3-0.8525261226466953) < 1e-5 , "Unscaled lambda3 computation doesn't match data."
    assert abs(l2-0.9743159081299232) < 1e-5 , "Unscaled lambda2 computation doesn't match data."
