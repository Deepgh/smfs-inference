#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 14:58:16 2025

@author: ghosh1
"""

import pandas as pd
import numpy as np
from smfs.estimation.distribution_fitting import poisson_pmf_model, unique_mutation_distribution
import os
import multiprocessing

number_of_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()))
#number_of_cores = multiprocessing.cpu_count()

data_path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
data_file = 'edited_sim_data_smp_mu_diff_site.csv'
df = pd.read_csv(os.path.join(data_path, data_file))

xi_path = data_path
xi_file = 'xi_mu_sim_data_diff_site.csv'

unq_mut_limit = 1000
output_filename = f'prob_sites_per_mut_{unq_mut_limit}_sim_data_diff_mu_site.csv'


xi_df = pd.read_csv(os.path.join(xi_path, xi_file))
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

df_result.to_csv(os.path.join(xi_path, output_filename), index=False)
    
    



