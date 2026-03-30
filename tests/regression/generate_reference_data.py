import pandas as pd
import os
import numpy as np
from smfs.estimation.distribution_fitting import *
from smfs.estimation.smfs_fit_loglik import compute as smfs_fit_loglik_compute
from smfs.estimation.smfs_fit_num import compute as smfs_fit_num_compute
# from examples.unq_mut_sites_per_mut_mew import generate_prob_sites_per_mutation

def generate_poisson(DATA_DIR, input_data, max_mutations=10):
    input_data_path = os.path.join(DATA_DIR, input_data)
    expected_xi_path = os.path.join(DATA_DIR, "expected_lambda.csv")
    expected_output_path = os.path.join(DATA_DIR, "expected_poisson.csv")

    data_df = pd.read_csv(input_data_path)
    
    np.random.seed(42)

    poisson_pmf = poisson_pmf_model(data_df["Mutation rate"])

    xi = fit_model_scaling(
        p0_func=lambda xi: poisson_pmf(0, xi),
        is_seg=data_df["AC"] != 0,
        initial_xi=1e8
    )
    expected_site_counts = expected_sites_per_mutation_count(
        site_mutation_probability=lambda k: poisson_pmf(k, xi),
        max_mutations=max_mutations
    )
    df1 = pd.DataFrame({
        'Unique mutation': np.arange(len(expected_site_counts)),
        'Unq muts sites': expected_site_counts
    })

    if os.path.exists(expected_xi_path):
        df = pd.read_csv(expected_xi_path)
    else:
        df = pd.DataFrame()
        
    df['Xi poisson'] = xi
    df.to_csv(expected_xi_path, index=False)
    df1.to_csv(expected_output_path, index=False)
    

def generate_gamma(DATA_DIR, input_data, max_mutations=10):
    input_data_path = os.path.join(DATA_DIR, input_data)
    expected_xi_path = os.path.join(DATA_DIR, "expected_lambda.csv")
    expected_output_path = os.path.join(DATA_DIR, "expected_gamma.csv")

    data_df = pd.read_csv(input_data_path)
    np.random.seed(42)

    poisson_gamma_pmf = poisson_gamma_pmf_model(data_df["Mean theta"], data_df["SD theta"])   

    xi = fit_model_scaling(
        p0_func=lambda xi: poisson_gamma_pmf(0, xi),
        is_seg=data_df["AC"] != 0,
        initial_xi=1e3
    )
    expected_site_counts = expected_sites_per_mutation_count(
        site_mutation_probability=lambda k: poisson_gamma_pmf(k, xi),
        max_mutations=max_mutations
    )
    df1 = pd.DataFrame({
        'Unique mutation': np.arange(len(expected_site_counts)),
        'Unq muts sites': expected_site_counts
    })

    if os.path.exists(expected_xi_path):
        df = pd.read_csv(expected_xi_path)
    else:
        df = pd.DataFrame()

    df['Xi gamma'] = xi
    df.to_csv(expected_xi_path, index=False)
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

def generate_prob_sites_per_mutation_reference(DATA_DIR, input_data, xi_data, unq_mut_limit):
    input_data_path = os.path.join(DATA_DIR, input_data)
    xi_data_path = os.path.join(DATA_DIR, xi_data)

    df = pd.read_csv(input_data_path)
    xi_df = pd.read_csv(xi_data_path)
    xi = xi_df['Xi poisson'].iloc[0]

    grouped_rates = (
        df.groupby('Mutation number')['Mutation rate']
        .apply(np.array) 
        .to_dict()
    )

    pmf_model_by_mut_type = {
        mut_type: poisson_pmf_model(rates)
        for mut_type, rates in grouped_rates.items()
    }

    df_result_dict = {
        mut_type: unique_mutation_distribution(lambda k: pmf_model_by_mut_type[mut_type](k, xi), unq_mut_limit)
        for mut_type in pmf_model_by_mut_type
    }   

    df_result = pd.DataFrame({'Unique mutation': range(unq_mut_limit + 1)})
    for mut_num, probs in df_result_dict.items():
        df_result[f"Mutation number {mut_num}"] = probs

    df_result.to_csv(os.path.join(DATA_DIR, f"poisson_per_mut_type.csv"), index=False)

def generate_smfs_loglik_reference(DATA_DIR, input_data, sites_file, ac_limit=10):
    data_file = os.path.join(DATA_DIR, input_data)
    sites_file_path = os.path.join(DATA_DIR, sites_file)
    output_path = os.path.join(DATA_DIR, "smfs_loglik_reference.csv")

    data_df = pd.read_csv(data_file)
    sites_df = pd.read_csv(sites_file_path)

    # Fit the model to get the expected probabilities for each unique mutation
    prob_ac_one_unq_mut = smfs_fit_loglik_compute(data_df, sites_df, ac_limit)

    # Create a DataFrame with the results
    df = pd.DataFrame({
        'AC': np.arange(1, ac_limit + 1),
        'pred prob': prob_ac_one_unq_mut
    })

    print(prob_ac_one_unq_mut)

    # Save the results to a CSV file
    df.to_csv(output_path, index=False)

if __name__ == "__main__":
    DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data")
    # generate_poisson(DATA_DIR, "input_data.csv")
    # generate_gamma(DATA_DIR, "input_data.csv")
    # generate_smfs_num_reference(DATA_DIR, "input_data.csv", "expected_poisson.csv", ac_limit=10)
    # generate_prob_sites_per_mutation_reference(DATA_DIR, "input_data.csv", "expected_lambda.csv", unq_mut_limit=10)
    # generate_smfs_loglik_reference(DATA_DIR, "input_data.csv", "poisson_per_mut_type.csv", ac_limit=10)

