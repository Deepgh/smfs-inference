#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 16:21:56 2025

@author: ghosh1
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import os
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

data_path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
data_file = 'edited_sim_data_smp_mu_diff_site.csv'

prob_sites_path = data_path
sites_file = 'prob_sites_per_mut_1000_sim_data_mu_diff_site.csv'

ac_limit = 1000
output_path = data_path
output_filename  = 'smfs_log_lik_' + str(ac_limit) + '_sim_data_mu_diff_site.csv'

data_df = pd.read_csv(os.path.join(data_path, data_file))
#raw_df = raw_df[raw_df['Meth level']!=0]

sites_df = pd.read_csv(os.path.join(prob_sites_path, sites_file))

z_init=  np.array([0])
prob_one_unq = sites_df[sites_df['Unique mutation'] == 1].values.tolist()[0][1:]
prob_non_seg = sites_df[sites_df['Unique mutation'] == 0].values.tolist()[0][1:]
prob_i_or_more = [1-x for x in prob_non_seg]

def log_it_bounded(z, prob_ac_one_unq_mut:list):
    return (1-sum(prob_ac_one_unq_mut))/(1+np.exp(-z))

def estimate_non_singleton(z, i_count:list, noni_count:list, prob_one_unq:list, p_recurrent:list,prob_i_or_more:list):
    return -(sum(i_count*np.log((log_it_bounded(z, prob_ac_one_unq_mut)*prob_one_unq + p_recurrent)/prob_i_or_more))+ 
             sum( noni_count*np.log(1-(log_it_bounded(z, prob_ac_one_unq_mut)*prob_one_unq + p_recurrent)/prob_i_or_more)))

prob_ac_one_unq_mut =[]

for ac in range(0, ac_limit): 
    
    i_count,  noni_count = [], []                 
    for mut_typ in range(1,145):
        mut_typ_df = data_df[data_df['Mutation number'] == mut_typ]  
        #print(i, mut_typ, len(mut_typ_df))
        if len(mut_typ_df)==0:
           ac_i_count = 0
           ac_noni_count = 0
        else:
          ac_i_count = len(mut_typ_df[mut_typ_df['AC'] == ac+1])
          ac_noni_count = len(mut_typ_df[mut_typ_df['AC'] > ac+1])
        i_count.append(ac_i_count)
        noni_count.append(ac_noni_count)   
        # print(ac_noni_count+ac_i_count)
        # print('                                    ')
    
    p_recurrent = [0]*144
    if ac>0:    
        conv_lists=[]
        f_convolved = prob_ac_one_unq_mut  
        print(ac+1)
       
        for j in range(ac):
            f_convolved = np.convolve(prob_ac_one_unq_mut, f_convolved)
            conv_lists.append(f_convolved) 
        #print(a) 
        
        l=len(conv_lists)     
        
        j = ac + 1
        k = 0
        while l >= 1:
            #print(f"  j: {j}, k: {k}",l)
            prob_sites_j_copies = sites_df[sites_df['Unique mutation'] == j].values.tolist()[0][1:]
            conv_sum = [conv_lists[l-1][k]*x for x in prob_sites_j_copies]
            p_recurrent = [x+y for x,y in zip(p_recurrent, conv_sum )]
            #print(i, j,l,k)
            # print('''''''''''''''''''')
            # print(a[l-1])
            # print(a[l-1][k],j)
            j -= 1  # Decrease j
            k += 1  # Increase k
            l=l-1
        
    res = minimize(estimate_non_singleton, z_init, method='nelder-mead',
                    args=(i_count,  noni_count, prob_one_unq, p_recurrent,prob_i_or_more),
                    options={'xatol': 1e-8, 'disp': True})
    
    val = log_it_bounded(res.x[0], prob_ac_one_unq_mut)
    
    recurrent_upadate = [val*x for x in prob_one_unq]
    p_recurrent_update = [x+y  for x ,y in zip(p_recurrent,recurrent_upadate)]
    prob_i_or_more = [x-y for x,y in zip(prob_i_or_more,p_recurrent_update)]
    
    prob_ac_one_unq_mut.append(val)
    #print(prob_ac_one_unq_mut)
    

new_df = pd.DataFrame()
for i in range(len(prob_ac_one_unq_mut)):
    data = pd.DataFrame([{'Count': i+1, 'pred prob':prob_ac_one_unq_mut[i]}])
    new_df = pd.concat([new_df, data], ignore_index=True)
new_df.to_csv(os.path.join(output_path, output_filename), index=False)



