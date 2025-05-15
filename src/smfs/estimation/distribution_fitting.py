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
from scipy.special import gammaln
from pathlib import Path
from typing import Sequence, Union
import os


def estimate_lambda_poisson(mu_seg, mu_non_seg, initial_lambda=1e8):
   
    def estimate(log_lam, mu_seg, mu_non_seg):
        return -(sum(-(np.exp(log_lam))*mu_non_seg)+sum(np.log(1-np.exp(-(np.exp(log_lam))*mu_seg)))) 

    res = minimize(estimate, 
        np.log(initial_lambda), 
        method='nelder-mead',
        args=(mu_seg, mu_non_seg),
        options={'xatol': 1e-8, 'disp': True}
    )
    return np.exp(res.x[0])


def compute_poisson_distribution(
    mutation_rates: Union[Sequence[float], np.ndarray],
    lambd: float,
    max_mut: int,
    stop_on_zero: bool = True
) -> pd.DataFrame:
    """
    Compute total Poisson-distributed site probabilities for mutation counts 0..max_mut.

    Parameters
    ----------
    mutation_rates : array-like
        Per-site mutation rates.
    lambd : float
        Poisson λ scaling factor.
    max_mut : int
        Max mutation count (inclusive).
    stop_on_zero : bool
        Stop early if the total becomes zero.

    Returns
    -------
    pd.DataFrame
        With columns: ["Unique mutation", "Unq muts sites"]
    """

    rates = np.asarray(mutation_rates, dtype=float)
    result = []

    for k in range(max_mut + 1):
        p_sum = poisson.pmf(k, lambd * rates).sum()
        result.append({'Unique mutation': k, 'Unq muts sites': p_sum})
        if stop_on_zero and p_sum == 0:
            break

     # Add padding if needed // KW: Is this necessary?
    last_k = result[-1]['Unique mutation']
    result.extend(
        {'Unique mutation': j, 'Unq muts sites': 0.0}
        for j in range(last_k + 1, max_mut + 1)
    )
    return pd.DataFrame(result) # Should probably be a pd.Series


def estimate_lambda_gamma(mu_seg, mu_non_seg, sig_seg, sig_non_seg, initial_lambda):
    def p0(lam,mu,sig):
        return np.exp(-(1+(mu/sig)**2)*np.log(1+(sig**2)*(lam/mu)))

    # Define the function to minimize
    def objective_function(log_lam, mu_seg, mu_non_seg, sig_seg, sig_non_seg):
        return -(sum(np.log(p0(np.exp(log_lam),mu_non_seg,sig_non_seg)))+sum(np.log(1-p0(np.exp(log_lam),mu_seg,sig_seg))))

    # Use the Nelder-Mead method to find the minimum
    result = minimize(
        objective_function, 
        np.log(initial_lambda), 
        method='nelder-mead', 
        args = (mu_seg, mu_non_seg, sig_seg, sig_non_seg),
        options={'xatol': 1e-8, 'disp': True}
    )
    return np.exp(result.x[0])


def compute_gamma_distribution(
    mean_theta: Union[np.ndarray, list],
    sd_theta: Union[np.ndarray, list],
    lambd: float,
    unq_mut_limit: int,
    stop_on_zero: bool = True
) -> pd.DataFrame:
    """
    Compute total Gamma-distributed site probabilities for mutation counts 0..unq_mut_limit.
    
    Parameters
    ----------
    data_df : pd.DataFrame
        Must contain 'Mean theta' and 'SD theta' columns.
    lambd : float
        Lambda parameter for the Gamma distribution.
    unq_mut_limit : int
        Maximum number of unique mutations to compute.
    stop_on_zero : bool
        Whether to stop when total site probability hits 0.
    pad_zeros : bool
        Whether to fill remaining entries with 0s after stopping.
    
    Returns
    -------
    pd.DataFrame
        Columns: ['Unique mutation', 'Unq muts sites']
    """
    mu = np.asarray(mean_theta, dtype=float)
    sigma = np.asarray(sd_theta, dtype=float)

    mu_sd_rat = (mu / sigma) ** 2
    numerator = (lambd * mu) / mu_sd_rat

    results = []

    for k in range(unq_mut_limit + 1):
        # Compute the generalized Negative Binomial-like probability
        log_pmf = (
            gammaln(k + mu_sd_rat)
            - gammaln(1 + k)
            - gammaln(mu_sd_rat)
            + k * np.log(numerator)
            - (k + mu_sd_rat) * np.log(1 + numerator)
        )
        total_prob = np.exp(log_pmf).sum()
        results.append({'Unique mutation': k, 'Unq muts sites': total_prob})

        if stop_on_zero and total_prob == 0:
            break

    # Pad with zeros if needed
    last_k = results[-1]['Unique mutation']
    results.extend(
        {'Unique mutation': j, 'Unq muts sites': 0.0}
        for j in range(last_k + 1, unq_mut_limit + 1)
    )

    return pd.DataFrame(results)    


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

    lambd = estimate_lambda_poisson(mu_seg, mu_non_seg, initial_lambda=lam)
    print(lambd)
    lam_val = pd.DataFrame([{'Lambda':lambd}])
    save_csv(lam_val, output_path, output_lam_filename)

    #%%%
    # lam_df = pd.read_csv('/project/yuvalsim/Deep/project2/human_data/lam_human_data_all.csv')
    # lambd = list(lam_df['Lambda'])[0]
    df1 = compute_poisson_distribution(data_df['Mutation rate'], lambd, unq_mut_limit)    
    save_csv(df1, output_path, sites_unq_mut_filename)






















