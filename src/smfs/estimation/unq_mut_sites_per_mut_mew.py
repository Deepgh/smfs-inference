#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 14:58:16 2025

@author: ghosh1
"""

import numpy as np
import pandas as pd
from scipy.stats import poisson
import multiprocessing
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def single_iteration(args):
    iter, mut_rate, lambd = args
     
    prob_sites = poisson.pmf(iter, (lambd)*mut_rate)

    return np.mean(prob_sites)

def generate_prob_sites_per_mutation(df, lam_df, unq_mut_limit=1000, number_of_cores=1):
    """
    Generate the probability of unique mutation sites per mutation number.
    """
    lambd = list(lam_df['Lambda poisson'])[0]
    unique_mutations = range(unq_mut_limit+1)  # 0 to 10
    
    df_result = pd.DataFrame({'Unique mutation': unique_mutations})

    for mut_num in df['Mutation number'].unique():
        # mut_rate = (df.loc[df['Mutation number'] == mut_num, 'Sampled mu'].iloc[0])
        mut_rate = df.loc[df['Mutation number'] == mut_num, 'Mutation rate'].iloc[0]
        
        with multiprocessing.Pool(number_of_cores) as pool:
            results = pool.map(single_iteration, [(i, mut_rate, lambd) for i in unique_mutations])

        df_result[f'Mutation number {mut_num}'] = results
    return df_result

if __name__ == "__main__":
    import os
    number_of_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()))
    #number_of_cores = multiprocessing.cpu_count()

    data_path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
    data_file = 'edited_sim_data_smp_mu_diff_site.csv'
    df = pd.read_csv(os.path.join(data_path, data_file))

    lam_path = data_path
    lam_file = 'lam_mu_sim_data_diff_site.csv'

    unq_mut_limit = 1000
    output_filename = f'prob_sites_per_mut_{unq_mut_limit}_sim_data_diff_mu_site.csv'


    lam_df = pd.read_csv(os.path.join(lam_path, lam_file))
    
    df_result = generate_prob_sites_per_mutation(df, lam_df, unq_mut_limit=unq_mut_limit, number_of_cores=number_of_cores)
    df_result.to_csv(os.path.join(lam_path, output_filename), index=False)
    
    




