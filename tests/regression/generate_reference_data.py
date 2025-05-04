import pandas as pd
import os
import numpy as np
from smfs.estimation.estimate_lambda_mew import estimate_lambda, compute_poisson_distribution
from smfs.estimation.estimate_lambda_mew_sd import estimate_lambda_sd, compute_gamma_distribution

def generate_poisson(DATA_DIR, input_data):
    input_data_path = os.path.join(DATA_DIR, input_data)
    expected_lambda_path = os.path.join(DATA_DIR, "expected_lambda.csv")
    expected_output_path = os.path.join(DATA_DIR, "expected_poisson.csv")

    data_df = pd.read_csv(input_data_path)

    mu_seg = data_df[data_df['AC']!=0]['Mutation rate'].to_numpy()
    mu_non_seg = data_df[data_df['AC']==0]['Mutation rate'].to_numpy()

    lambd = estimate_lambda(mu_seg, mu_non_seg, initial_lambda=1e8)
    lam_val = pd.DataFrame([{'Lambda poisson':lambd}])
    lam_val.to_csv(expected_lambda_path, index=False)

    unq_mut_limit = 10
    np.random.seed(42)
    df1 = compute_poisson_distribution(data_df, lambd, unq_mut_limit)
    df1.to_csv(expected_output_path, index=False)
    

def generate_gamma(DATA_DIR, input_data):
    input_data_path = os.path.join(DATA_DIR, input_data)
    expected_lambda_path = os.path.join(DATA_DIR, "expected_lambda.csv")
    expected_output_path = os.path.join(DATA_DIR, "expected_gamma.csv")

    data_df = pd.read_csv(input_data_path)

    mu_seg = data_df[data_df['AC']!=0]['Mean theta'].to_numpy()
    mu_non_seg = data_df[data_df['AC']==0]['Mean theta'].to_numpy()
    sig_seg = data_df[data_df['AC']!=0]['SD theta'].to_numpy()
    sig_non_seg = data_df[data_df['AC']==0]['SD theta'].to_numpy()

    lambd = estimate_lambda_sd(mu_seg, mu_non_seg, sig_seg, sig_non_seg, initial_lambda=1e3)
    print('Estimated lambda:', lambd)   

    if os.path.exists(expected_lambda_path):
        df = pd.read_csv(expected_lambda_path)
    else:
        df = pd.DataFrame()
    df['Lambda gamma'] = lambd
    df.to_csv(expected_lambda_path, index=False)

    unq_mut_limit = 10
    np.random.seed(42)
    df1 = compute_gamma_distribution(data_df, lambd, unq_mut_limit)
    df1.to_csv(expected_output_path, index=False)


if __name__ == "__main__":
    DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data")
    generate_poisson(DATA_DIR, "input_data.csv")
    generate_gamma(DATA_DIR, "input_data.csv")

