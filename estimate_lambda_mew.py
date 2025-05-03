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
import math
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
data_file = 'edited_sim_data_smp_mu_diff_site.csv'
data_df = pd.read_csv(os.path.join(path,data_file))

output_path = path
output_lam_filename = 'lam_mu_sim_data_diff_site.csv'


Path(output_path).mkdir(parents=True, exist_ok=True)

#%%
# df_3 = data_df[data_df['Meth level']==0]
#data_df = data_df[data_df['Meth level']!=0]
#data_df = data_df[data_df["Meth le"]>=10**-7]
lam= np.array([10**8])
mu_seg = list(data_df[data_df['AC']!=0]['Sampled mu'])
mu_non_seg = list(data_df[data_df['AC']==0]['Sampled mu'])

print(len(mu_seg), len(mu_non_seg), (len(mu_seg)+len(mu_non_seg))/10**6)

#%%    
x, y = [],[]
def estimate(log_lam, mu_seg, mu_non_seg):
    # x.append(log_lam)
    # y.append(sum(-(np.exp(log_lam))*mu_seg)+sum(np.log(1-np.exp(-(np.exp(log_lam))*mu_non_seg))))
    # plt.plot(x,y)
    # plt.show()
    # plt.clf()    
    return -(sum(-(np.exp(log_lam))*mu_non_seg)+sum(np.log(1-np.exp(-(np.exp(log_lam))*mu_seg)))) 

res = minimize(estimate, np.log(lam), method='nelder-mead',args=(mu_seg, mu_non_seg),
               options={'xatol': 1e-8, 'disp': True})

print(math.exp(res.x[0]))
lam_val = pd.DataFrame([{'Lambda':math.exp(res.x[0])}])
lam_val.to_csv(os.path.join(output_path, output_lam_filename), index= False)
#%%%
unq_mut_limit = 1000
sites_unq_mut_filename = 'unq_mut_sites_mu_'+str(unq_mut_limit)+'_sim_data_diff_site.csv'
lambd = (math.exp(res.x[0]))

data_df = pd.read_csv(os.path.join(path,data_file))
# lam_df = pd.read_csv('/project/yuvalsim/Deep/project2/human_data/lam_human_data_all.csv')
# lambd = list(lam_df['Lambda'])[0]
df1 = pd.DataFrame()
stop_filling = False

for unq_mut in range(unq_mut_limit + 1):
    print(unq_mut)    
    if not stop_filling:
        col_name = 'Unq muts sites_' + str(unq_mut)
        data_df[col_name] = poisson.pmf(unq_mut, lambd * data_df['Mutation rate'])
        poi_sum = sum(data_df[col_name])
        print(poi_sum)
        
        if poi_sum == 0:
            stop_filling = True
    else:
        poi_sum = 0  # After first zero, just fill 0s without any calculations
    
    data = pd.DataFrame([{'Unique mutation': unq_mut, 'Unq muts sites': poi_sum}])
    df1 = pd.concat([df1, data], ignore_index=True)

df1.to_csv(os.path.join(output_path, sites_unq_mut_filename), index=False)























