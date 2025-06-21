#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 00:12:49 2024

@author: ghosh1
"""

import numpy as np
import pandas as pd

def compute(data_df, df_pmf, ac_limit):
    eps_1 = df_pmf.iloc[1]['Unq muts sites']
    num_cop_1 = len(data_df[data_df['AC']==1]['AC'])
    prop_unq_muts_1 = num_cop_1/eps_1
    con_list = [prop_unq_muts_1]
    f_2 = np.convolve(con_list, con_list) 
    con_list.append(f_2[0])
    f2 = (len(data_df[data_df['AC']==2]['AC'])-(f_2[0]*df_pmf.iloc[2]['Unq muts sites']))/eps_1
    f = [prop_unq_muts_1, f2]   
    for ac in range(2, ac_limit): 
        conv_lists=[]
        f_convolved = f  
        for j in range(ac):
            f_convolved = np.convolve(f, f_convolved)
            conv_lists.append(f_convolved) 
        s1=0
        l=len(conv_lists)     
        j = ac+1
        k = 0
        while l >= 1:
            s1= s1 + conv_lists[l-1][k]*df_pmf.iloc[j]['Unq muts sites']
            j -= 1  # Decrease j
            k += 1  # Increase k
            l=l-1
        val= (len(data_df[data_df['AC']==ac+1]['AC'])-s1)/eps_1
        f.append(val)
    return f

if __name__ == "__main__":
    import os
    from pathlib import Path
    path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
    data_file = 'edited_sim_data_smp_mu_diff_site.csv'
    unq_mut_sites_path = path
    unq_mut_sites_file = 'unq_mut_sites_mu_1000_sim_data_diff_site.csv'
    ac_limit = 1000
    output_path = path
    output_filename = f'smfs_num_{ac_limit}_sim_data_mu_diff_site.csv'
    Path(output_path).mkdir(parents=True, exist_ok=True)
    data_df = pd.read_csv(os.path.join(path, data_file))
    pmf_data_file = os.path.join(unq_mut_sites_path, unq_mut_sites_file)
    df_pmf = pd.read_csv(pmf_data_file)
    f = compute(data_df, df_pmf, ac_limit)
    new_df = pd.DataFrame()
    for i in range(len(f)):
        data = pd.DataFrame([{'Count': i+1, 'pred prob':f[i]}])
        new_df = pd.concat([new_df, data], ignore_index=True)
    new_df.to_csv(os.path.join(output_path, output_filename), index=False)

