import pandas as pd
import os
import pytest
from smfs.estimation.estimate_lambda_mew import estimate_lambda, compute_poisson_distribution

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
    input_df, _, expected_lambda = poisson_data
    mu_seg = input_df[input_df["AC"] != 0]["Mutation rate"].to_numpy()
    mu_non_seg = input_df[input_df["AC"] == 0]["Mutation rate"].to_numpy()

    lam = estimate_lambda(mu_seg, mu_non_seg, method="poisson")
    assert round(lam, 8) == round(expected_lambda, 8)

def test_poisson_distribution(poisson_data):
    input_df, expected_output, expected_lambda = poisson_data
    unq_mut_limit = len(expected_output) - 1

    result = compute_poisson_distribution(input_df, expected_lambda, unq_mut_limit)

    pd.testing.assert_frame_equal(
        result.round(8),
        expected_output.round(8),
        rtol=1e-8,
        check_dtype=False
    )

def test_gamma_lambda(gamma_data):
    input_df, _, expected_lambda = gamma_data
    mu_seg = input_df[input_df["AC"] != 0]["Mutation rate"].to_numpy()
    mu_non_seg = input_df[input_df["AC"] == 0]["Mutation rate"].to_numpy()

    lam = estimate_lambda(mu_seg, mu_non_seg, method="gamma")
    assert round(lam, 8) == round(expected_lambda, 8)

def test_gamma_distribution(gamma_data):
    input_df, expected_output, expected_lambda = gamma_data
    unq_mut_limit = len(expected_output) - 1

    result = compute_poisson_distribution(input_df, expected_lambda, unq_mut_limit)

    pd.testing.assert_frame_equal(
        result.round(8),
        expected_output.round(8),
        rtol=1e-8,
        check_dtype=False
    )
