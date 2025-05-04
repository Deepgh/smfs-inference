#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 16:43:15 2024

@author: ghosh1
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import poisson
from pathlib import Path
import os


def estimate_lambda(mu_seg, mu_non_seg, log_lam0):

    #%%    
    def estimate(log_lam, mu_seg, mu_non_seg):
        return -(sum(-(np.exp(log_lam))*mu_non_seg)+sum(np.log(1-np.exp(-(np.exp(log_lam))*mu_seg)))) 

    res = minimize(estimate, log_lam0, method='nelder-mead',args=(mu_seg, mu_non_seg),
        options={'xatol': 1e-8, 'disp': True})
    return np.exp(res.x[0])


def compute_poisson_distribution(data_df, lambd, unq_mut_limit):
    df1 = pd.DataFrame()
    stop_filling = False

    for unq_mut in range(unq_mut_limit + 1):
        print(unq_mut)    
        if not stop_filling:
            col_name = 'Unq muts sites_{unq_mut}'
            data_df[col_name] = poisson.pmf(unq_mut, lambd * data_df['Mutation rate'])
            poi_sum = sum(data_df[col_name])
            print(poi_sum)
            
            if poi_sum == 0:
                stop_filling = True
        else:
            poi_sum = 0  # After first zero, just fill 0s without any calculations
        
        data = pd.DataFrame([{'Unique mutation': unq_mut, 'Unq muts sites': poi_sum}])
        df1 = pd.concat([df1, data], ignore_index=True)
    return df1


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
    from warnings import simplefilter
    simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

    path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
    data_file = 'edited_sim_data_smp_mu_diff_site.csv'

    output_path = path
    output_lam_filename = 'lam_mu_sim_data_diff_site.csv'

    unq_mut_limit = 1000
    sites_unq_mut_filename = f'unq_mut_sites_mu_{unq_mut_limit}_sim_data_diff_site.csv'


    data_df = load_data(path, data_file)
        #%%
    # df_3 = data_df[data_df['Meth level']==0]
    #data_df = data_df[data_df['Meth level']!=0]
    #data_df = data_df[data_df["Meth le"]>=10**-7]
    lam= np.array([10**8])
    mu_seg = data_df[data_df['AC']!=0]['Sampled mu'].to_numpy()
    mu_non_seg = data_df[data_df['AC']==0]['Sampled mu'].to_numpy()

    print(len(mu_seg), len(mu_non_seg), (len(mu_seg)+len(mu_non_seg))/10**6)

    lambd = estimate_lambda(mu_seg, mu_non_seg, initial_lambda=lam)
    print(lambd)
    lam_val = pd.DataFrame([{'Lambda':lambd}])
    save_csv(lam_val, output_path, output_lam_filename)

    #%%%
    # lam_df = pd.read_csv('/project/yuvalsim/Deep/project2/human_data/lam_human_data_all.csv')
    # lambd = list(lam_df['Lambda'])[0]
    df1 = compute_poisson_distribution(data_df, lambd, unq_mut_limit)    
    save_csv(df1, output_path, sites_unq_mut_filename)






















