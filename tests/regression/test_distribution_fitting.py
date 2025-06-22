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

def test_prob_sites_per_mutation_regression():
    reference_path = os.path.join(DATA_DIR, "poisson_per_mut_type.csv")
    input_data_path = os.path.join(DATA_DIR, "input_data.csv")
    lambda_data_path = os.path.join(DATA_DIR, "expected_lambda.csv")
    ref_df = pd.read_csv(reference_path)
    input_df = pd.read_csv(input_data_path)
    lambda_df = pd.read_csv(lambda_data_path)
    # Use the same unq_mut_limit as in the reference file
    unq_mut_limit = ref_df['Unique mutation'].max()
    lambda_val = lambda_df['Lambda poisson'].iloc[0]

    grouped_rates = (
        input_df.groupby('Mutation number')['Mutation rate']
        .apply(np.array)
        .to_dict()
    )

    # Compute probabilities
    pmf_model_by_mut_type = {
        mut_type: poisson_pmf_model(rates)
        for mut_type, rates in grouped_rates.items()
    }

    result_dict = {
        mut_type : unique_mutation_distribution(lambda k: pmf_model(k, lambda_val), unq_mut_limit)
        for mut_type, pmf_model in pmf_model_by_mut_type.items()
    }   

    # Convert result for comparison
    result_df = pd.DataFrame({'Unique mutation': np.arange(unq_mut_limit + 1)})
    for mut_num, probs in result_dict.items():
        result_df[f"Mutation number {mut_num}"] = probs

    # Validate
    pd.testing.assert_frame_equal(
        result_df.round(8),
        ref_df.round(8),
        rtol=1e-8,
        check_dtype=False
    )
