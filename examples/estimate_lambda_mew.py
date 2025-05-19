
import os
import pandas as pd
from pathlib import Path
import numpy as np

from smfs.estimation.distribution_fitting import poisson_pmf_model

def load_data(data_path, filename):
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data path {data_path} does not exist.")
    if not os.path.exists(os.path.join(data_path, filename)):   
        raise FileNotFoundError(f"File {filename} does not exist in {data_path}.")
    return pd.read_csv(os.path.join(data_path, filename))


def save_csv(df, output_path, filename, overwrite=False):
    if not overwrite and os.path.exists(os.path.join(output_path, filename)):
        raise FileExistsError(f"{filename} already exists in {output_path}. Use overwrite=True to replace it.")
    Path(output_path).mkdir(parents=True, exist_ok=True)
    df.to_csv(os.path.join(output_path, filename), index=False)


if __name__ == "__main__":
    from smfs.estimation.distribution_fitting import fit_model_scaling, expected_sites_per_mutation_count, poisson_pmf_model

    path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
    data_file = 'edited_sim_data_smp_mu_diff_site.csv'

    output_path = path
    output_lam_filename = 'lam_mu_sim_data_diff_site.csv'

    unq_mut_limit = 1000
    sites_unq_mut_filename = f'unq_mut_sites_mu_{unq_mut_limit}_sim_data_diff_site.csv'


    data_df = load_data(path, data_file)

    mu_seg = data_df[data_df['AC']!=0]['Sampled mu'].to_numpy()
    mu_non_seg = data_df[data_df['AC']==0]['Sampled mu'].to_numpy()
    print(len(mu_seg), len(mu_non_seg), (len(mu_seg)+len(mu_non_seg))/10**6)

    poisson_pmf = poisson_pmf_model(data_df["Sampled mu"])

    lambd = fit_model_scaling(
        p0_func=lambda lam: poisson_pmf(0, lam),
        is_seg= data_df["AC"] != 0,
        initial_lambda=1e8
    )
    expected_site_counts = expected_sites_per_mutation_count(
        site_mutation_probability= lambda k: poisson_pmf(k, lambd),
        max_mutations=unq_mut_limit
    )
    result = pd.DataFrame({
        'Unique mutation': np.arange(len(expected_site_counts)),
        'Unq muts sites': expected_site_counts
    })

    print(lambd)
    lam_val = pd.DataFrame([{'Lambda':lambd}])
    save_csv(lam_val, output_path, output_lam_filename)
    save_csv(result, output_path, sites_unq_mut_filename)

