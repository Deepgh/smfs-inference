#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 17:47:52 2025

@author: ghosh1
"""

import pandas as pd
import os

der_allele_path = '/project/yuvalsim/Deep/project2/josh_full_data/error_data/simulation_results/sampled_data/'
der_allel_file = 'der_allele_sum_per_site_smp_mu_diff_site.csv.gz'
mut_rate_path = '/project/yuvalsim/Deep/project2/josh_full_data/error_data/simulation_results/'
mut_rate_file = 'mut_rates_gamma_smp_diff_site.csv'

output_path = mut_rate_path
output_file = 'edited_sim_data_smp_mu_diff_site.csv.gz'
os.makedirs(output_path, exist_ok=True)

der_allele_df = pd.read_csv(os.path.join(der_allele_path, der_allel_file), compression='gzip')
mut_rate_df =  pd.read_csv(os.path.join(mut_rate_path, mut_rate_file))

merged_df = mut_rate_df.merge(
    der_allele_df[['Site', 'Sampled mutation sum', 'Unique mutations']],
    on='Site',
    how='left')

merged_df[['Sampled mutation sum', 'Unique mutations']] = merged_df[
    ['Sampled mutation sum', 'Unique mutations']].fillna(0).astype(int)

merged_df = merged_df.rename(columns={'Sampled mutation sum': 'AC'})
merged_df['Site']+=1

merged_df.to_csv(os.path.join(output_path, output_file), compression='gzip', index=False)


