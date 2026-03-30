#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 16:59:07 2026

@author: ghosh1
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from pathlib import Path
import os
from scipy.stats import poisson
import math
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

# -------------------- Paths and files --------------------
path = '/project/yuvalsim/Deep/project2/josh_full_data/error_data/simulation_results/'
data_file = 'edited_sim_data_smp_mu_diff_site.csv.gz'
output_path = path
output_lam_filename = 'lam_smp_mu.csv'

Path(output_path).mkdir(parents=True, exist_ok=True)

# -------------------- Column names --------------------
sites_col = 'Site'
mu_col = 'Sampled mu'
ac_col = 'AC'

# -------------------- Load data --------------------
data_df = pd.read_csv(os.path.join(path, data_file), 
                      #sep='\t'
                      compression='gzip'
                      )


# -------------------- Repeat rows by number of sites ------------------
#data_df = data_df.loc[data_df.index.repeat(data_df[sites_col]), [mu_col, ac_col]].reset_index(drop=True)

# -------------------- Separate segregating and non-segregating sites --------------------
mu_seg = np.array(data_df[data_df[ac_col] != 0][mu_col])
mu_non_seg = np.array(data_df[data_df[ac_col] == 0][mu_col])


# -------------------- Log-likelihood function --------------------
def estimate(log_lam, mu_seg, mu_non_seg):
    """
    Negative log-likelihood for non-error version.
    """
    lam_val = np.exp(log_lam)
    
    # Numerical stability: avoid log(0)
    eps = 1e-12
    p_seg = np.exp(-lam_val * mu_seg)
    p_seg = np.clip(p_seg, eps, 1 - eps)
    
    p_non = np.exp(-lam_val * mu_non_seg)
    p_non = np.clip(p_non, eps, 1 - eps)
    
    log_lik = np.sum(np.log(p_non)) + np.sum(np.log(1 - p_seg))
    return -log_lik

# -------------------- Minimize negative log-likelihood --------------------
initial_guess = np.log([1e3])  # initial lambda
res = minimize(
    estimate,
    initial_guess,
    method='nelder-mead',
    args=(mu_seg, mu_non_seg),
    options={'xatol': 1e-8, 'disp': True}
)

# -------------------- Output results --------------------
lambda_est = float(np.exp(res.x[0]))
print("Estimated lambda:", lambda_est)

lam_val = pd.DataFrame([{'Lambda': lambda_est}])
lam_val.to_csv(os.path.join(output_path, output_lam_filename), index=False)
#%%

unq_mut_limit = 300
sites_unq_mut_filename = 'unq_mut_sites_'+str(unq_mut_limit)+'_sim_smp_mu.csv'
lambd = float(math.exp(res.x[0]))

df1 = pd.DataFrame()
stop_filling = False

for unq_mut in range(unq_mut_limit + 1):
    print(unq_mut)    
    if not stop_filling:
        col_name = 'Unq muts sites_' + str(unq_mut)
        data_df[col_name] = (poisson.pmf(unq_mut, lambd * data_df[mu_col]))
        poi_sum = (sum(data_df[col_name]))
        print(poi_sum)
        
        if poi_sum == 0:
            stop_filling = True
    else:
        poi_sum = 0  # After first zero, just fill 0s without any calculations
    
    data = pd.DataFrame([{'Unique mutation': unq_mut, 'Unq muts sites': poi_sum}])
    df1 = pd.concat([df1, data], ignore_index=True)

df1.to_csv(os.path.join(output_path, sites_unq_mut_filename), index=False)
