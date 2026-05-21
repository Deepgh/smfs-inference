#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 17:32:08 2025

@author: ghosh1
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

raw_data = '/project/yuvalsim/Deep/project2/josh_full_data/error_data/ajhg_00004094_supp_table2_mut.tsv'
df = pd.read_csv(raw_data, sep='\t')

df = df.drop(['LL_no_heterogeneity', 'LL_heterogeneity', 'CV', 'LR','ref_context', 
              'alt_context','p_value', 'p_adj'], axis=1)

df['scaled mu'] = df['mean_theta']
df['scaled sd'] = df['sd_theta']

# Compute shape and scale
df['k'] = (df['scaled mu'] / df['scaled sd']) ** 2
df['theta'] = (df['scaled sd'] ** 2) / df['scaled mu']

# Generate random seeds (e.g., between 0 and 1e9)
num_rows = len(df)
df['seed'] = np.random.randint(0, int(1e9), size=len(df)).astype(int)

# Sampling function
def sample_mean(k, theta, seed, n=100000):
    seed = int(seed)  # <- force it to be an integer
    rng = np.random.default_rng(seed)
    samples = rng.gamma(shape=k, scale=theta, size=n)
    return pd.Series({'sampled mu': np.mean(samples), 'sampled sd': np.std(samples)})

#%%
df['scaled mu'] = df['mean_theta'] 
df['scaled sd'] = df['sd_theta'] 

# Compute Gamma shape and scale
df['k'] = (df['scaled mu'] / df['scaled sd']) ** 2
df['theta'] = (df['scaled sd'] ** 2) / df['scaled mu']
df['seed'] = np.random.randint(0, int(1e9), size=len(df))
expanded_rows = []

mut_num = 1
for _, row in df.iterrows():
    print(mut_num)
    seed = int(row['seed'])
    rng = np.random.default_rng(seed)
    n = int(row['num_sites'])

    # Draw gamma samples
    samples = rng.gamma(shape=row['k'], scale=row['theta'], size=n)

    # Create a new DataFrame with repeated metadata and sampled values
    long_df = pd.DataFrame({'Mutation number':[mut_num]*n,
        'Mutation rate': [row['mu_gnomad']] * n,
        'Mean theta': [row['mean_theta']] * n,
        'SD theta': [row['sd_theta']] * n,
        'Meth_level': [row['methylation_level']]*n,
        'Sampled mu': samples
    })
    mut_num += 1
    expanded_rows.append(long_df)

# Concatenate all expanded rows into one big DataFrame
result_df = pd.concat(expanded_rows, ignore_index=True)
result_df.to_csv('/project/yuvalsim/Deep/project2/josh_full_data/error_data/scaled_mu_diff_site.csv', index=False)



#%%

# Apply per row
df[['sampled mu', 'sampled sd']] = df.apply(lambda row: sample_mean(row['k'], row['theta'], row['seed']), axis=1)

# Optional cleanup
#df.drop(columns=['k', 'theta', 'seed'], inplace=True)

df['diff'] = (df['sampled mu']- df['scaled mu'])/df['scaled mu']
plt.plot(range(len(df['diff'])), df['diff'])
plt.ylabel('relative error sampled and scaled mu')
plt.show()

plt.plot(df['scaled mu'], df['sampled mu'])
plt.xlabel('scaled mu')
plt.ylabel('mean sampled mu')
plt.show()

df['diffsd'] = (df['scaled sd']- df['scaled sd'])/df['scaled sd']
plt.plot(range(len(df['diffsd'])), df['diffsd'])
plt.ylabel('relative error sampled and scaledsd')
plt.show()

plt.plot(df['scaled sd'], df['sampled sd'])
plt.xlabel('scaled sd')
plt.ylabel('sampleld sd')
plt.show()

df.drop(columns=['k', 'theta', 'seed', 'diff', 'diffsd'], inplace=True)

df.to_csv('/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/scaled_mu_from_musd.csv', index=False)


