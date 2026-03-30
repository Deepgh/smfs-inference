import numpy as np
import pandas as pd
import os
import pytest
import src.smfs.estimation.smfs_fit_num as smfs_fit_num
import src.smfs.estimation.smfs_fit_loglik as smfs_fit_loglik

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data")

def test_smfs_fit_num_regression():
    data_file = os.path.join(DATA_DIR, "input_data.csv")
    sites_file = os.path.join(DATA_DIR, "expected_poisson.csv")
    ref_file = os.path.join(DATA_DIR, "smfs_num_reference.csv")
    data_df = pd.read_csv(data_file)
    sites_df = pd.read_csv(sites_file)
    ac_limit = 10
    allele_counts = data_df['AC'].values
    pmf_unq_mut_sites = sites_df['Unq muts sites'].values
    result = smfs_fit_num.compute(allele_counts, pmf_unq_mut_sites, ac_limit)
    ref_df = pd.read_csv(ref_file)
    expected = ref_df['pred prob'].values
    assert np.allclose(result, expected, rtol=1e-3, atol=1e-6), f"Result: {result}, Expected: {expected}"


def test_smfs_fit_num_weighted_rows_matches_expanded_rows():
    allele_counts = np.array([1, 2, 3])
    site_weights = np.array([5, 2, 1])
    expanded_counts = np.repeat(allele_counts, site_weights)
    expected_sites_per_mutation_count = np.array([0.0, 10.0, 3.0, 1.0, 0.5])
    ac_limit = 3

    weighted_result = smfs_fit_num.compute(
        allele_counts,
        expected_sites_per_mutation_count,
        ac_limit,
        site_weights=site_weights,
    )
    expanded_result = smfs_fit_num.compute(
        expanded_counts,
        expected_sites_per_mutation_count,
        ac_limit,
    )

    assert np.allclose(weighted_result, expanded_result)

    
def test_smfs_fit_loglik_regression():
    data_file = os.path.join(DATA_DIR, "input_data.csv")
    sites_file = os.path.join(DATA_DIR, "poisson_per_mut_type.csv")
    ref_file = os.path.join(DATA_DIR, "smfs_loglik_reference.csv")
    data_df = pd.read_csv(data_file)
    sites_df = pd.read_csv(sites_file)
    ac_limit = 10
    result = smfs_fit_loglik.compute_sfms(data_df, sites_df, ac_limit)
    ref_df = pd.read_csv(ref_file)
    expected = ref_df['pred prob'].values
    assert np.allclose(result, expected, rtol=1e-3, atol=1e-6), f"Result: {result}, Expected: {expected}"
