import pandas as pd
import os
import numpy as np
from smfs.estimation.distribution_fitting import *
from smfs.estimation.smfs_fit_num import compute as smfs_fit_num_compute
from smfs.estimation.unq_mut_sites_per_mut_mew import generate_prob_sites_per_mutation

def generate_poisson(DATA_DIR, input_data, max_mutations=10):
    input_data_path = os.path.join(DATA_DIR, input_data)
    expected_lambda_path = os.path.join(DATA_DIR, "expected_lambda.csv")
    expected_output_path = os.path.join(DATA_DIR, "expected_poisson.csv")

    data_df = pd.read_csv(input_data_path)
    
    np.random.seed(42)

    poisson_pmf = poisson_pmf_model(data_df["Mutation rate"])

    lambd = fit_model_scaling(
        p0_func=lambda lam: poisson_pmf(0, lam),
        is_seg= data_df["AC"] != 0,
        initial_lambda=1e8
    )
    expected_site_counts = expected_sites_per_mutation_count(
        site_mutation_probability= lambda k: poisson_pmf(k, lambd),
        max_mutations=max_mutations
    )
    df1 = pd.DataFrame({
        'Unique mutation': np.arange(len(expected_site_counts)),
        'Unq muts sites': expected_site_counts
    })

    if os.path.exists(expected_lambda_path):
        df = pd.read_csv(expected_lambda_path)
    else:
        df = pd.DataFrame()
        
    df['Lambda poisson'] = lambd
    df.to_csv(expected_lambda_path, index=False)
    df1.to_csv(expected_output_path, index=False)
    

def generate_gamma(DATA_DIR, input_data, max_mutations=10):
    input_data_path = os.path.join(DATA_DIR, input_data)
    expected_lambda_path = os.path.join(DATA_DIR, "expected_lambda.csv")
    expected_output_path = os.path.join(DATA_DIR, "expected_gamma.csv")

    data_df = pd.read_csv(input_data_path)
    np.random.seed(42)

    poisson_gamma_pmf = poisson_gamma_pmf_model(data_df["Mean theta"], data_df["SD theta"])   

    lambd = fit_model_scaling(
        p0_func=lambda lam: poisson_gamma_pmf(0, lam),
        is_seg= data_df["AC"] != 0,
        initial_lambda=1e3
    )
    expected_site_counts = expected_sites_per_mutation_count(
        site_mutation_probability= lambda k: poisson_gamma_pmf(k, lambd),
        max_mutations=max_mutations
    )
    df1 = pd.DataFrame({
        'Unique mutation': np.arange(len(expected_site_counts)),
        'Unq muts sites': expected_site_counts
    })

    if os.path.exists(expected_lambda_path):
        df = pd.read_csv(expected_lambda_path)
    else:
        df = pd.DataFrame()

    df['Lambda gamma'] = lambd
    df.to_csv(expected_lambda_path, index=False)
    df1.to_csv(expected_output_path, index=False)

def generate_smfs_num_reference(DATA_DIR, input_data, sites_file, ac_limit=10):
    data_file = os.path.join(DATA_DIR, input_data)
    sites_file_path = os.path.join(DATA_DIR, sites_file)
    output_path = os.path.join(DATA_DIR, "smfs_num_reference.csv")
    data_df = pd.read_csv(data_file)
    sites_df = pd.read_csv(sites_file_path)
    allele_counts = data_df['AC'].values
    expected_sites_per_mutation_count = sites_df['Unq muts sites'].values
    result = smfs_fit_num_compute(allele_counts, expected_sites_per_mutation_count, ac_limit)
    df = pd.DataFrame({
        'AC': np.arange(1, ac_limit+1),
        'pred prob': result
    })
    df.to_csv(output_path, index=False)

def generate_prob_sites_per_mutation_reference(DATA_DIR, input_data, lambda_data, unq_mut_limit):
    input_data_path = os.path.join(DATA_DIR, input_data)
    lambda_data_path = os.path.join(DATA_DIR, lambda_data)
    output_path = os.path.join(DATA_DIR, f"poisson_per_mut_type.csv")

    df = pd.read_csv(input_data_path)
    lam_df = pd.read_csv(lambda_data_path)
    df_result = generate_prob_sites_per_mutation(df, lam_df, unq_mut_limit=unq_mut_limit, number_of_cores=1)
    df_result.to_csv(output_path, index=False)

if __name__ == "__main__":
    DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data")
    # generate_poisson(DATA_DIR, "input_data.csv")
    # generate_gamma(DATA_DIR, "input_data.csv")
    # generate_smfs_num_reference(DATA_DIR, "input_data.csv", "expected_poisson.csv", ac_limit=10)
    generate_prob_sites_per_mutation_reference(DATA_DIR, "input_data.csv", "expected_lambda.csv", unq_mut_limit=10)

