import numpy as np
import pandas as pd
import os
import pytest
import src.smfs.estimation.smfs_fit_num as smfs_fit_num

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data")

def test_smfs_fit_num_regression():
    data_file = os.path.join(DATA_DIR, "input_data.csv")
    sites_file = os.path.join(DATA_DIR, "expected_poisson.csv")
    ref_file = os.path.join(DATA_DIR, "smfs_num_reference.csv")
    data_df = pd.read_csv(data_file)
    sites_df = pd.read_csv(sites_file)
    ac_limit = 10
    result = smfs_fit_num.compute(data_df, sites_df, ac_limit)
    ref_df = pd.read_csv(ref_file)
    expected = ref_df['pred prob'].values
    assert np.allclose(result, expected, rtol=1e-3, atol=1e-6), f"Result: {result}, Expected: {expected}"

