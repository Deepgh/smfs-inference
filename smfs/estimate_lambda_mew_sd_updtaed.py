#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 14:37:33 2026

@author: ghosh1
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from pathlib import Path
import os
from scipy.stats import nbinom
from warnings import simplefilter

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

# -------------------------
# Paths and input
# -------------------------
path = '/project/yuvalsim/Deep/project2/josh_full_data/error_data/'
data_file = 'sfs_genes_dist.csv.gz'

data_df = pd.read_csv(os.path.join(path, data_file),
    compression='gzip')

output_path = path
Path(output_path).mkdir(parents=True, exist_ok=True)

# -------------------------
# Column names
# -------------------------
sites = 'Sites'
mean = 'Mean theta'
sd = 'SD theta'
ac = 'AC_nfe_down'

# -------------------------
# Initial guess for lambda
# -------------------------
lam0 = np.array([1e3])

# -------------------------
# Expand dataframe to per-site rows
# If in your data each site is row then this is not needed
# -------------------------
# data_df = (
#     data_df
#     .loc[data_df.index.repeat(data_df[sites]), [ac, mean, sd]]
#     .reset_index(drop=True)
# )

# -------------------------
# Split segregating / non-segregating sites
# -------------------------
mu_seg = data_df.loc[data_df[ac] != 0, mean].values
sig_seg = data_df.loc[data_df[ac] != 0, sd].values

mu_non = data_df.loc[data_df[ac] == 0, mean].values
sig_non = data_df.loc[data_df[ac] == 0, sd].values

# -------------------------
# Probability of zero mutational events
# (Negative binomial marginal)
# -------------------------
def p0(lam, mu, sig):
    mu = np.asarray(mu)
    sig = np.asarray(sig)

    mu_sig_sq = (mu / sig) ** 2
    return np.exp(-mu_sig_sq * np.log(1 + (lam * mu) / mu_sig_sq))

# -------------------------
# Negative log-likelihood
# -------------------------
def neg_loglik(log_lam, mu_seg, mu_non, sig_seg, sig_non):
    lam = np.exp(log_lam)

    p0_non = p0(lam, mu_non, sig_non)
    p0_seg = p0(lam, mu_seg, sig_seg)

    # Numerical safety
    eps = 1e-12
    p0_non = np.clip(p0_non, eps, 1 - eps)
    p0_seg = np.clip(p0_seg, eps, 1 - eps)

    ll = (
        np.sum(np.log(p0_non)) +
        np.sum(np.log(1 - p0_seg))
    )

    return -ll

# -------------------------
# Optimize
# -------------------------
res = minimize(
    neg_loglik,
    x0=np.log(lam0),
    args=(mu_seg, mu_non, sig_seg, sig_non),
    method='nelder-mead',
    options={'xatol': 1e-8, 'disp': True}
)

# -------------------------
# Output
# -------------------------
lam_hat = float(np.exp(res.x[0]))
print(f"Estimated lambda: {lam_hat:}")

lam_df = pd.DataFrame([{'Lambda': lam_hat}])
#lam_df.to_csv(os.path.join(output_path, 'lam_mu_sd_josh.csv'), index=False)


#%%

unq_mut_limit = 300
sites_unq_mut_filename = f'unq_mut_sites_{unq_mut_limit}_josh_mew_sd.csv'

df1 = pd.DataFrame()
stop_filling = False

for unq_mut in range(unq_mut_limit + 1):
    print(unq_mut)

    if not stop_filling:
        col_name = f'Unq muts sites_{unq_mut}'

        data_df[col_name] = nbinom.pmf(
            unq_mut,
            n=(data_df[mean] / data_df[sd]) ** 2,
            p=1 / (1 + data_df[sd] ** 2 * (lam_hat / data_df[mean])))

        poi_sum = data_df[col_name].sum()
        print(poi_sum)

        if poi_sum == 0:
            stop_filling = True
    else:
        poi_sum = 0

    df1 = pd.concat([df1, pd.DataFrame([{
            'Unique mutation': unq_mut,
            'Unq muts sites': poi_sum}])],ignore_index=True)

df1.to_csv(os.path.join(output_path, sites_unq_mut_filename),index=False)

