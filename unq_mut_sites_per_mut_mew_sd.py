#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 14:58:16 2025

@author: ghosh1
"""

import numpy as np
import pandas as pd
import os
from scipy.special import gammaln
import multiprocessing
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
data_file = 'edited_sim_data_smp_mu_diff_site.csv'
df = pd.read_csv(os.path.join(path,data_file))

lam_path = path
lam_file = 'lam_mu_sd_sim_data_diff_site.csv'


unq_mut_limit = 1000
output_filename = 'prob_sites_per_mut_'+str(unq_mut_limit)+'_sim_data_mu_sd_diff_site.csv'


lam_df = pd.read_csv(os.path.join(lam_path, lam_file))
lambd = list(lam_df['Lambda'])[0]

df['mu_sd_rat'] = (df['Mean theta'] / df['SD theta']) ** 2
df['numerator factor'] = (lambd * df['Mean theta']) / df['mu_sd_rat']

def single_iteration(args):
    iter, mu_sd_rat, numerator_factor = args

    prob_sites = np.exp(gammaln(iter + mu_sd_rat) - gammaln(1 + iter) - gammaln(mu_sd_rat) + 
                         iter * np.log(numerator_factor) - (iter + mu_sd_rat) * 
                         np.log(1 + numerator_factor))

    return prob_sites

if __name__ == "__main__":
    number_of_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()))
      
    unique_mutations = range(unq_mut_limit+1)  # 0 to 10
    
    df_result = pd.DataFrame({'Unique mutation': unique_mutations})

    for mut_num in range(1, 145):
        mu_sd_rat = (df.loc[df['Mutation number'] == mut_num, 'mu_sd_rat'].iloc[0])
        numerator_factor = (df.loc[df['Mutation number'] == mut_num, 'numerator factor'].iloc[0])
        
        with multiprocessing.Pool(number_of_cores) as pool:
            results = pool.map(single_iteration, [(i,mu_sd_rat, numerator_factor) for i in unique_mutations])

        df_result['Mutation number ' + str(mut_num)] = results

    df_result.to_csv(os.path.join(lam_path, output_filename), index=False)