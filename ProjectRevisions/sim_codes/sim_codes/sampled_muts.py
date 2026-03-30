#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 10:41:22 2024

@author: ghosh1
"""

import pandas as pd
from scipy.stats import hypergeom
import os


path =  '/project/yuvalsim/Deep/project2/josh_full_data/error_data/'
file = 'der_allel_gamma_smp_mu_diff_site.csv'

output_dir = os.path.join(path,'sampled_data')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
df = pd.read_csv(os.path.join(path,file))

pop_final = 2*(10**8)

sample_size = 2*5*(10**5)
df['Sampled mutation'] = hypergeom.rvs(pop_final, df['Allele copies'], sample_size)

df_sampled =df[df['Sampled mutation'] !=0]
        
df_sampled.to_csv(os.path.join(output_dir, 'samp_data_10_5_5_smp_mu_diff_site.csv.gz'),compression='gzip', index=False)
    
