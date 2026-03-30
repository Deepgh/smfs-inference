#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 21:56:02 2026

@author: ghosh1
"""

import numpy as np
import pandas as pd
import os
from pathlib import Path

path = '/project/yuvalsim/Deep/project2/josh_full_data/error_data/simulation_results/'
data_file = 'edited_sim_data_smp_mu_diff_site.csv.gz'

data_df = pd.read_csv(os.path.join(path,data_file), 
                 #sep='\t'
                 compression='gzip'
                 )

ac = 'AC'

#%%
unq_mut_sites_path = path
unq_mut_sites_file = 'unq_mut_sites_300_sim_mew_sd.csv'

ac_limit = 300
output_path = unq_mut_sites_path 
output_filename = 'smfs_num_'+str(ac_limit)+'_sim_mew_sd.csv'

Path(output_path).mkdir(parents=True, exist_ok=True)

pmf_data_file = os.path.join(unq_mut_sites_path, unq_mut_sites_file)
df_pmf = pd.read_csv(pmf_data_file)

eps_1 = df_pmf.iloc[1]['Unq muts sites']
num_cop_1 = len(data_df[data_df[ac]==1][ac])


prop_unq_muts_1 = num_cop_1/eps_1

con_list = [prop_unq_muts_1]


f_2 = np.convolve(con_list, con_list) 
con_list.append(f_2[0])
f2 = (len(data_df[data_df[ac]==2][ac])-(f_2[0]*df_pmf.iloc[2]['Unq muts sites']))/eps_1

            
f = [prop_unq_muts_1, f2]   

for acs in range(2, ac_limit): 
    conv_lists=[]
    f_convolved = f  
    print(acs+1)
   
    for j in range(acs):
        f_convolved = np.convolve(f, f_convolved)
        conv_lists.append(f_convolved) 
    #print(a) 
    s1=0
    l=len(conv_lists)     
    
    j = acs+1
    k = 0
    while l >= 1:
        s1= s1 + conv_lists[l-1][k]*df_pmf.iloc[j]['Unq muts sites']
       
        j -= 1  # Decrease j
        k += 1  # Increase k
        l=l-1
    val= (len(data_df[data_df[ac]==acs+1][ac])-s1)/eps_1
   
    
    #print(len(data_df[data_df['AC']==i+1]['AC']), s1, val)
    f.append(val)

new_df = pd.DataFrame()

for i in range(len(f)):
    data = pd.DataFrame([{'Count': i+1, 'pred prob':f[i]}])
    new_df = pd.concat([new_df, data], ignore_index=True)

new_df.to_csv(os.path.join(output_path, output_filename), index=False)

