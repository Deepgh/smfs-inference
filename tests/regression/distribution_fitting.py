import unittest
import pandas as pd
import os
from smfs.estimation.estimate_lambda_mew import estimate_lambda, compute_poisson_distribution

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data")

class TestPoissonRegression(unittest.TestCase):

    def setUp(self):
        self.input_df = pd.read_csv(os.path.join(DATA_DIR, "input_data.csv"))
        self.expected_output = pd.read_csv(os.path.join(DATA_DIR, "expected_poisson.csv"))
        self.expected_lambda = pd.read_csv(os.path.join(DATA_DIR, "expected_lambda.csv"))
        self.expected_lambda = self.expected_lambda["Lambda poisson"].values[0]

    def test_lambda_regression(self):
        mu_seg = self.input_df[self.input_df["AC"] != 0]["Mutation rate"].to_numpy()
        mu_non_seg = self.input_df[self.input_df["AC"] == 0]["Mutation rate"].to_numpy()

        lam = estimate_lambda(mu_seg, mu_non_seg)
        self.assertAlmostEqual(lam, self.expected_lambda, places=8)

    def test_poisson_regression(self):
        unq_mut_limit = len(self.expected_output) - 1

        result = compute_poisson_distribution(self.input_df, self.expected_lambda, unq_mut_limit)

        pd.testing.assert_frame_equal(
            result.round(8),
            self.expected_output.round(8),
            rtol=1e-8,
            check_dtype=False
        )

class TestGammaRegression(unittest.TestCase):

    def setUp(self):
        self.input_df = pd.read_csv(os.path.join(DATA_DIR, "input_data.csv"))
        self.expected_output = pd.read_csv(os.path.join(DATA_DIR, "expected_gamma.csv"))
        self.expected_lambda = pd.read_csv(os.path.join(DATA_DIR, "expected_lambda.csv"))
        self.expected_lambda = self.expected_lambda["Lambda gamma"].values[0]

    def test_lambda_regression(self):
        mu_seg = self.input_df[self.input_df["AC"] != 0]["Mutation rate"].to_numpy()
        mu_non_seg = self.input_df[self.input_df["AC"] == 0]["Mutation rate"].to_numpy()

        lam = estimate_lambda(mu_seg, mu_non_seg)
        self.assertAlmostEqual(lam, self.expected_lambda, places=8)

    def test_poisson_regression(self):
        unq_mut_limit = len(self.expected_output) - 1

        result = compute_poisson_distribution(self.input_df, self.expected_lambda, unq_mut_limit)

        pd.testing.assert_frame_equal(
            result.round(8),
            self.expected_output.round(8),
            rtol=1e-8,
            check_dtype=False
        )


if __name__ == "__main__":
    unittest.main()
