#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 16:43:15 2024

@author: ghosh1
"""

import pandas as pd
import numpy as np
from pathlib import Path
import os

from smfs.estimation.distribution_fitting import fit_model_scaling, expected_sites_per_mutation_count, poisson_gamma_pmf_model


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

    path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
    data_df = load_data(path, 'edited_sim_data_smp_mu_diff_site.csv')

    output_path = path

    unq_mut_limit = 1000

    poisson_gamma_pmf = poisson_gamma_pmf_model(data_df["Mean theta"], data_df["SD theta"])

    lambda0 = fit_model_scaling(
        p0_func=lambda lam: poisson_gamma_pmf(0, lam),
        is_seg= data_df["AC"] != 0,
        initial_lambda=1e8
    )
    expected_site_counts = expected_sites_per_mutation_count(
        site_mutation_probability= lambda k: poisson_gamma_pmf(k, lambda0),
        max_mutations=unq_mut_limit
    )
    df1 = pd.DataFrame({
        'Unique mutation': np.arange(len(expected_site_counts)),
        'Unq muts sites': expected_site_counts
    })

    print('Estimated lambda:', lambda0)   
    lam_val = pd.DataFrame([{'Lambda':lambda0}])

    save_csv(lam_val, output_path, 'lam_all_sim_data_smp_mu_diff_site.csv')
    save_csv(df1, output_path, f'unq_mut_sites_{unq_mut_limit}_ukbb_data_mew_sd.csv')




















