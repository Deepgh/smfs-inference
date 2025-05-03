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
import math
import matplotlib.pyplot as plt
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

path = '/project/yuvalsim/Deep/project2/human_data/sim_hum_data/final_approach/smp_mu_diff_site/seeded'
data_file = 'edited_sim_data_smp_mu_diff_site.csv'
data_df = pd.read_csv(os.path.join(path,data_file))

output_path = path
output_lam_filename = 'lam_all_sim_data_smp_mu_diff_site.csv'

Path(output_path).mkdir(parents=True, exist_ok=True)

# df_3 = data_df[data_df['Meth level']==0]
#data_df = data_df[data_df['Meth level']!=0]
# print(len(df_3)/10**7, len(df_2)/10**5, (len(df_3)-len(df_2))/10**7)#
#data_df = data_df[data_df["Mew"]>=1.135*10**-8]
lam= np.array([10**3])
mu_seg = list(data_df[data_df['AC']!=0]['Mean theta'])
mu_non_seg = list(data_df[data_df['AC']==0]['Mean theta'])
sig_seg = list(data_df[data_df['AC']!=0]['SD theta'])
sig_non_seg = list(data_df[data_df['AC']==0]['SD theta'])

x, y = [],[]
def p0(lam,mu,sig):
    sig = np.array(sig)  # Convert to NumPy array if it's a list
    mu = np.array(mu)
    return np.exp(-(1+(mu/sig)**2)*np.log(1+(sig**2)*(lam/mu)))

def estimate(log_lam, mu_seg, mu_non_seg, sig_seg, sig_non_seg):
    # x.append(log_lam)
    # y.append(sum(np.log(p0(np.exp(log_lam),mu_non_seg,sig_non_seg)))+sum(np.log(1-p0(np.exp(log_lam),mu_seg,sig_seg))))
    # plt.plot(x,y)
    # plt.show()
    # plt.clf() 
    return -(sum(np.log(p0(np.exp(log_lam),mu_non_seg,sig_non_seg)))+sum(np.log(1-p0(np.exp(log_lam),mu_seg,sig_seg)))) 
        

res = minimize(estimate, np.log(lam), method='nelder-mead',args=(mu_seg, mu_non_seg, sig_seg, sig_non_seg),
               options={'xatol': 1e-8, 'disp': True})

print(math.exp(res.x[0]))
lam_val =pd.DataFrame([{'Lambda':math.exp(res.x[0])}])
lam_val.to_csv(os.path.join(output_path, output_lam_filename), index= False)
#%%%
unq_mut_limit = 1000
sites_unq_mut_filename = 'unq_mut_sites_'+str(unq_mut_limit)+'_ukbb_data_mew_sd.csv'
lambd = (math.exp(res.x[0]))
# lam_file = 'lam_all.csv'
# lam_df = pd.read_csv(os.path.join(path,lam_file))
# lambd = list(lam_df['Lambda'])[0]

data_df['mu_sd_rat'] =  (data_df['Mean theta'] / data_df['SD theta'])** 2
data_df['numerator factor'] = (lambd*data_df['Mean theta'])/data_df['mu_sd_rat']

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

df1.to_csv(os.path.join(output_path,sites_unq_mut_filename), index=False)






















