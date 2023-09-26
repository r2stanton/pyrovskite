from pyrovskite.featurizer import Featurizer
import pandas as pd
import os


try:
    dirr = os.path.join(os.getcwd(), "data")
    df = pd.read_pickle(f"{dirr}/test_df.pkl")
except:
    dirr = os.path.join(os.getcwd(), "tests/data")
    df = pd.read_pickle(f"{dirr}/test_df.pkl")

def test_organic_inorganic_ratio():
    ft = Featurizer(df = df, ase_col = 'ase_atoms')
    ft.org_inorg_ratio()
    assert abs(ft.df['org_inorg_ratio'].iloc[0] - 0.6) < 1e-6 , "Test data doesn't match reference data"
    assert abs(ft.df['org_inorg_ratio'].iloc[1] - 0.636364) < 1e-6 , "Test data doesn't match reference data"
    assert abs(ft.df['org_inorg_ratio'].iloc[2] - 0.666667) < 1e-6 , "Test data doesn't match reference data"

def test_organic_inorganic_weight_ratio():
    ft = Featurizer(df = df, ase_col = 'ase_atoms')
    ft.org_inorg_weight_ratio()
    assert abs(ft.df['org_inorg_weight_ratio'].iloc[0] - 0.207976) < 1e-6 , "Test data doesn't match reference data"
    assert abs(ft.df['org_inorg_weight_ratio'].iloc[1] - 0.203185) < 1e-6 , "Test data doesn't match reference data"
    assert abs(ft.df['org_inorg_weight_ratio'].iloc[2] - 0.198316) < 1e-6 , "Test data doesn't match reference data"
    assert abs(ft.df['org_inorg_weight_ratio'].iloc[3] - 0.316702) < 1e-6 , "Test data doesn't match reference data"

