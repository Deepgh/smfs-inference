#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 11:15:24 2025

@author: ghosh1
"""

import numpy as np
import pandas as pd
import os
import multiprocessing
import warnings
warnings.filterwarnings('ignore')

path = '/project/yuvalsim/Deep/project2/josh_full_data/error_data/simulation_results/sampled_data/'
file = 'samp_data_10_5_5_smp_mu_diff_site.csv.gz'
df = pd.read_csv(os.path.join(path, file), compression='gzip')

def single_iteration(group):
    iter = group['Site'].iloc[0]
    unq_mut = len(group['Mutation id'].unique())
    meth_level = group['Meth level'].iloc[0]
    mut_rate = group['Mutation rate'].iloc[0]
    mut_num = group['Mutation number'].iloc[0]
    mew_gamma = group['Mean theta'].iloc[0]
    sd_gamma = group['SD theta'].iloc[0]
    smp_mew = group['Sampled mu'].iloc[0]
    
    mut_sum = group['Sampled mutation'].sum()
    
    return {'Mutation number': mut_num,
            'Site': iter,
            'Mutation rate': mut_rate,
            'Mew theta': mew_gamma,
            'SD theta': sd_gamma,
            'Sampled mu': smp_mew,
            'Meth level': meth_level,
            'Sampled mutation sum': mut_sum,
            'Unique mutations': unq_mut}


if __name__ == "__main__":
    number_of_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()))
    #number_of_cores = multiprocessing.cpu_count()  # if not on the cluster you should do this instead
    
    # Group by 'Site' to process data more efficiently
    grouped_df = df.groupby('Site')
    
    # Using multiprocessing to parallelize the iteration
    with multiprocessing.Pool(number_of_cores) as pool:
        results = pool.map(single_iteration, [group for _, group in grouped_df])
    
    # Convert the results list into a DataFrame
    df_sum_sites = pd.DataFrame(results)
    
    df_sum_sites.to_csv(os.path.join(path, 'der_allele_sum_per_site_smp_mu_diff_site.csv.gz'), compression='gzip',index=False)
