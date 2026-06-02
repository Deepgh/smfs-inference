#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 16:09:00 2025

@author: ghosh1
"""

import numpy as np
import os
import pandas as pd
from scipy.stats import poisson, nbinom
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import gammaln


path = '/project/yuvalsim/Deep/project2/josh_full_data/error_data/'
muts_data_file = 'sfs_syn_err_ac_1M_freq_dist.csv.gz'
pred_prob_file = 'smfs_num_300_err_josh_mew_sd.csv'


pred_df =  pd.read_csv(os.path.join(path,pred_prob_file))
muts_data_df = pd.read_csv(os.path.join(path, muts_data_file), compression='gzip')

muts_data_df = muts_data_df.rename(columns={'AC_nfe_down': 'AC'})
prob_counts = list(pred_df['pred prob'])
prob_counts.insert(0,0)

conv_lists =[]
conv_lists.append(prob_counts)
conv = prob_counts

for i in range(201):
    conv = np.convolve(prob_counts, conv)
    conv_lists.append(conv)

#%%
lam_file = 'lam_mu_sd_josh_err.csv'
lam_df = pd.read_csv(os.path.join(path,lam_file))
lambd = list(lam_df['Lambda'])[0]


# output_dir_1 = path+'stopgain_results_new_data/sites_ets_data_mu/'
# if not os.path.exists(output_dir_1):
#     os.makedirs(output_dir_1)
# output_dir_2 = path+'stopgain_results_new_data/copies_mu/'
# if not os.path.exists(output_dir_2):
#     os.makedirs(output_dir_2)

muts_data_df['mu_sd_rat'] =  (muts_data_df['Mean theta'] / muts_data_df['SD theta'])** 2
muts_data_df['numerator factor'] = (lambd*muts_data_df['Mean theta'])/muts_data_df['mu_sd_rat']


mut = 1
prop_sites_df = pd.DataFrame()
all_sum_sites_data = pd.DataFrame()
for i in range(0,len(muts_data_df),((10**6)+1)):
    df = muts_data_df[i:i+((10**6)+1)]

    #mew = np.unique(df['Mew'])[0]
    gamma_mew = np.unique(df['Mean theta'])[0]
    gamma_sigma = np.unique(df['SD theta'])[0]
    mu_sigma_rat = np.unique(df['mu_sd_rat'])[0]
    mu_sig_nume = np.unique(df['numerator factor'])[0]
    meth = np.unique(df['Meth level'])[0]

    mean = np.unique(df['Mean theta'])[0]
    sd = np.unique(df['SD theta'])[0]
    # mu_sigma_rat = float(gamma_mew/gamma_sigma)**2
    # mu_sig_nume = (lambd*gamma_mew)/mu_sigma_rat

    df1= pd.DataFrame()
    tot_sites = sum(df['Sites'])

    site_0 = (df[df['AC']==0])['Sites'].iloc[0]
    pred_0 = tot_sites*(1/((1+mu_sig_nume)**mu_sigma_rat))
    #pred_0 = tot_all_sites*poisson.pmf(0, lambd*mew)


    pred_site_prop_0 = pred_0/tot_sites
    data_site_prop_0 = site_0/tot_sites

    prop_sites_data = pd.DataFrame([{'AC': 0, #"Mew": mew,
                                     "Meth level":meth,
                               'Gamma mew':gamma_mew,'Gamma sigma':gamma_sigma,
                               'Prop sites':data_site_prop_0, 'Type':'Data'}])
    prop_sites_theory = pd.DataFrame([{'AC': 0, #"Mew": mew,
                                       "Meth level":meth,
                               'Gamma mew':gamma_mew,'Gamma sigma':gamma_sigma,
                               'Prop sites':pred_site_prop_0, 'Type':'Theory'}])
    prop_sites_df = pd.concat([prop_sites_df, prop_sites_data, prop_sites_theory], ignore_index=True)


    data_data = pd.DataFrame([{'AC': 0, #"Mew": mew,
                               "Meth level":meth,
                               'Gamma mew':gamma_mew,'Gamma sigma':gamma_sigma,
                               'Num sites':site_0, 'Type':'Data'}])
    theory_data = pd.DataFrame([{'AC': 0, #"Mew": mew,
                                 "Meth level":meth,
                                 'Gamma mew':gamma_mew,'Gamma sigma':gamma_sigma,
                                 'Num sites':pred_0, 'Type':'Theory'}])
    df1 = pd.concat([df1, data_data, theory_data], ignore_index=True)


    for copy in range(1,201):

        num_site = (df[df['AC']==copy])['Sites'].iloc[0]

        pred_sum = 0
        for l in range(copy):

            #poi_coeff = poisson.pmf(l+1, lambd*mew)
            # poi_coeff = np.exp(gammaln(l+1+mu_sigma_rat) - gammaln(1+l+1) -
            #                     gammaln(mu_sigma_rat) +
            #                     (l+1)*np.log(mu_sig_nume) -
            #                     (l+1+mu_sigma_rat)*np.log(1+mu_sig_nume))
            poi_coeff= (nbinom.pmf(l+1, n=(mean / sd)**2,
                                           p=1/(1 + sd**2 *(lambd/mean ))))
            conv_prob = poi_coeff*conv_lists[l][copy]
            #print(copy,l)

            pred_sum += conv_prob



        pred_site = tot_sites*pred_sum

        data_data = pd.DataFrame([{'AC': copy, #"Mew": mew,
                                   "Meth level":meth,
                                   'Gamma mew':gamma_mew, 'Gamma sigma':gamma_sigma,
                                   'Num sites':num_site, 'Type':'Data'}])
        theory_data = pd.DataFrame([{'AC': copy, #"Mew": mew,
                                     "Meth level":meth,
                                     'Gamma mew':gamma_mew, 'Gamma sigma':gamma_sigma,
                                     'Num sites':pred_site, 'Type':'Theory'}])
        df1 = pd.concat([df1, data_data, theory_data], ignore_index=True)


        data_site_prop = num_site/tot_sites

        prop_sites_data = pd.DataFrame([{'AC': copy, #"Mew": mew,
                                         "Meth level":meth,
                                   'Gamma mew':gamma_mew,'Gamma sigma':gamma_sigma,
                                   'Prop sites':data_site_prop, 'Type':'Data'}])
        prop_sites_theory = pd.DataFrame([{'AC': copy, #"Mew": mew,
                                           "Meth level":meth,
                                   'Gamma mew':gamma_mew,'Gamma sigma':gamma_sigma,
                                   'Prop sites':pred_sum, 'Type':'Theory'}])
        prop_sites_df = pd.concat([prop_sites_df, prop_sites_data, prop_sites_theory], ignore_index=True)

    all_sum_sites_data = pd.concat([all_sum_sites_data, df1], ignore_index=True)

    #mew_title =format_scientific(mew)
    # plt.yscale('log')
    # plt.title('$\mu$ = '+str(mew_title),fontsize=16)
    # sns.barplot(df1, x="AC", y="Num sites", hue="Type", palette=['indianred', "steelblue"])
    # ytick_labels = ['${}^{{{}}}$'.format(10, i) for i in range(0,int(np.ceil(np.log10(np.max(df1["Num sites"])))))]
    # plt.yticks([10**k for k in range(0,int(np.ceil(np.log10(np.max(df1["Num sites"])))))],ytick_labels,fontsize=16)
    # plt.xticks([0,10,20,30,40, 50],fontsize=16)
    # plt.ylabel('Site count', fontsize=18)
    # plt.xlabel('Allele count', fontsize=18)
    # plt.minorticks_off()
    # sns.despine(right=True, top=True)
    # plt.gcf().subplots_adjust(bottom=0.16, left=0.16)
    # plt.savefig(os.path.join(output_dir_1, 'mut_'+str(mut)+'.png'), dpi=1000)
    #plt.show()

    print(mut)
    # df_raw_data = muts_data_df[i:i+(10**6)+1]
    # df2 = df_raw_data[df_raw_data['AC']<=200]
    # plt.yscale('log')
    # plt.title('$\mu$ = '+str(mew_title),fontsize=16)
    # sns.barplot(df2, x="AC", y="Sites",color='darkslateblue')
    # ytick_labels = ['${}^{{{}}}$'.format(10, i) for i in range(0,int(np.ceil(np.log10(np.max(df2["Sites"])))))]
    # plt.yticks([10**k for k in range(0,int(np.ceil(np.log10(np.max(df2["Sites"])))))],ytick_labels,fontsize=16)
    # plt.xticks([0,50,100,150,200],fontsize=16)
    # plt.ylabel('Site count', fontsize=18)
    # plt.xlabel('Allele count', fontsize=18)
    # plt.minorticks_off()
    # sns.despine(right=True, top=True)
    # plt.gcf().subplots_adjust(bottom=0.16, left=0.16)
    # plt.savefig(os.path.join(output_dir_2, 'copies_raw_'+str(mut)+'.png'), dpi=1000)
    # plt.show()
    mut+=1

all_sum_sites_data.to_csv(os.path.join(path, 'all_sum_sites_mu_sd_missense.csv'), index=False)        
