#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 16:43:15 2024

@author: ghosh1
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from pathlib import Path
from scipy.special import gammaln
import os


def estimate_lambda(mu_seg, mu_non_seg, sig_seg, sig_non_seg, initial_lambda):
    def p0(lam,mu,sig):
        return np.exp(-(1+(mu/sig)**2)*np.log(1+(sig**2)*(lam/mu)))

    # Define the function to minimize
    def objective_function(lam, mu_seg, mu_non_seg, sig_seg, sig_non_seg):
        return -(sum(np.log(p0(lam, mu_non_seg, sig_non_seg))) + sum(np.log(1 - p0(lam, mu_seg, sig_seg))))

    # Use the Nelder-Mead method to find the minimum
    result = minimize(
        objective_function, 
        initial_lambda, 
        method='nelder-mead', 
        args = (mu_seg, mu_non_seg, sig_seg, sig_non_seg),
        options={'xatol': 1e-8, 'disp': True}
    )
    return np.exp(result.x[0])


def compute_poisson_distribution(data_df, lambd, unq_mut_limit):
    df1 = pd.DataFrame()
    stop_filling = False

    for unq_mut in range(unq_mut_limit+1):
        print(unq_mut)
        if not stop_filling:
            col_name = 'Unq muts sites_'+str(unq_mut)
            
            #data_df['numerator'] = data_df['numerator factor']**i
            
            data_df[col_name] = np.exp(gammaln(unq_mut+data_df['mu_sd_rat']) - gammaln(1+unq_mut) -
                                        gammaln(data_df['mu_sd_rat']) + unq_mut*np.log(data_df['numerator factor'])
                                        - (unq_mut+data_df['mu_sd_rat'])*np.log(1+data_df['numerator factor']))
            #print(data_df[col_name])
            
            poi_sum = sum(data_df[col_name])
            print(poi_sum)
            if poi_sum == 0:
                stop_filling = True
        else:
            poi_sum = 0  # After first zero, just fill 0s without any calculations
        data =pd.DataFrame([{'Unique mutation': unq_mut, 'Unq muts sites':poi_sum}])
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
    data_df = load_data(path, data_file)

    output_path = path
    output_lam_filename = 'lam_all_sim_data_smp_mu_diff_site.csv'

    # df_3 = data_df[data_df['Meth level']==0]
    #data_df = data_df[data_df['Meth level']!=0]
    # print(len(df_3)/10**7, len(df_2)/10**5, (len(df_3)-len(df_2))/10**7)#
    #data_df = data_df[data_df["Mew"]>=1.135*10**-8]
    lam= np.array([10**3])
    mu_seg = data_df[data_df['AC']!=0]['Mean theta'].to_numpy()
    mu_non_seg = data_df[data_df['AC']==0]['Mean theta'].to_numpy()
    sig_seg = data_df[data_df['AC']!=0]['SD theta'].to_numpy()
    sig_non_seg = data_df[data_df['AC']==0]['SD theta'].to_numpy()


    lambd = estimate_lambda(mu_seg, mu_non_seg, sig_seg, sig_non_seg, initial_lambda=lam)
    print('Estimated lambda:', lambd)   
    lam_val =pd.DataFrame([{'Lambda':lambd}])
    save_csv(lam_val, output_path, output_lam_filename)

    #%%%
    unq_mut_limit = 1000
    sites_unq_mut_filename = 'unq_mut_sites_'+str(unq_mut_limit)+'_ukbb_data_mew_sd.csv'
    # lam_file = 'lam_all.csv'
    # lam_df = pd.read_csv(os.path.join(path,lam_file))
    # lambd = list(lam_df['Lambda'])[0]

    data_df['mu_sd_rat'] =  (data_df['Mean theta'] / data_df['SD theta'])** 2
    data_df['numerator factor'] = (lambd*data_df['Mean theta'])/data_df['mu_sd_rat']
    
    df1 = compute_poisson_distribution(data_df, lambd, unq_mut_limit)
    save_csv(df1, output_path, sites_unq_mut_filename)




















