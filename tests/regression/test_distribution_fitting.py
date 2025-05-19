import numpy as np
import pandas as pd
import os
import pytest
from smfs.estimation.distribution_fitting import *
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data")

@pytest.fixture
def poisson_data():
    input_df = pd.read_csv(os.path.join(DATA_DIR, "input_data.csv"))
    expected_output = pd.read_csv(os.path.join(DATA_DIR, "expected_poisson.csv"))
    expected_lambda = pd.read_csv(os.path.join(DATA_DIR, "expected_lambda.csv"))["Lambda poisson"].values[0]
    return input_df, expected_output, expected_lambda

@pytest.fixture
def gamma_data():
    input_df = pd.read_csv(os.path.join(DATA_DIR, "input_data.csv"))
    expected_output = pd.read_csv(os.path.join(DATA_DIR, "expected_gamma.csv"))
    expected_lambda = pd.read_csv(os.path.join(DATA_DIR, "expected_lambda.csv"))["Lambda gamma"].values[0]
    return input_df, expected_output, expected_lambda

def test_poisson_lambda(poisson_data):
    df, _, expected_lambda = poisson_data
    poisson_pmf = poisson_pmf_model(df['Mutation rate'])

    lam = fit_model_scaling(
        p0_func=lambda lam: poisson_pmf(0, lam),
        is_seg= df["AC"] != 0,
        initial_lambda=1e8
    )
    assert round(lam, 8) == round(expected_lambda, 8)

def test_poisson_distribution(poisson_data):
    df, expected_output, expected_lambda = poisson_data
    poisson_pmf = poisson_pmf_model(df['Mutation rate'])

    expected_site_counts = expected_sites_per_mutation_count(
        site_mutation_probability= lambda k: poisson_pmf(k, expected_lambda),
        max_mutations=len(expected_output) - 1
    )
    result = pd.DataFrame({
        'Unique mutation': np.arange(len(expected_site_counts)),
        'Unq muts sites': expected_site_counts
    })

    pd.testing.assert_frame_equal(
        result.round(8),
        expected_output.round(8),
        rtol=1e-8,
        check_dtype=False
    )

def test_gamma_lambda(gamma_data):
    input_df, _, expected_lambda = gamma_data
    gamma_pmf = poisson_gamma_pmf_model(input_df['Mean theta'], input_df['SD theta'])

    lam = fit_model_scaling(
        p0_func= lambda lam: gamma_pmf(0, lam),
        is_seg= input_df["AC"] != 0,
        initial_lambda=1e3
    )
    assert round(lam, 8) == round(expected_lambda, 8)

def test_gamma_distribution(gamma_data):
    input_df, expected_output, expected_lambda = gamma_data
    gamma_pmf = poisson_gamma_pmf_model(input_df['Mean theta'], input_df['SD theta'])

    expected_site_counts = expected_sites_per_mutation_count(
        site_mutation_probability= lambda k: gamma_pmf(k, expected_lambda),
        max_mutations=len(expected_output) - 1
    )
    result = pd.DataFrame({
        'Unique mutation': np.arange(len(expected_site_counts)),
        'Unq muts sites': expected_site_counts
    })
    pd.testing.assert_frame_equal(
        result.round(8),
        expected_output.round(8),
        rtol=1e-8,
        check_dtype=False
    )
