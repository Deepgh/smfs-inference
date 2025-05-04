#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 15:00:36 2025

@author: ghosh1
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 00:43:33 2024

@author: deepj
"""

import pandas as pd
import numpy as np

file = '/project/yuvalsim/Deep/project2/human_data/human_data_sfs.csv'

df= pd.read_csv(file)

sfs_list=[]

for i in range(0,len(df),1001):
    anc_read_list = (df['ref_context'])
    der_read_list = (df['alt_context'])
    anc_read = anc_read_list[i]
    der_read = der_read_list[i]
    
    mut_sum = sum(df['n'][i:i+1001])
    
    comb = [anc_read, der_read, mut_sum]
    
    sfs_list.append(comb)
    
#%%

file2= '/project/yuvalsim/Deep/project2/human_data/supp_table_mut.csv'    

df2 = pd.read_csv(file2) 
mut_rate_list=[]

for i in range(0,len(df2)):
    anc_read_list = (df2['ref_context'])
    der_read_list = (df2['alt_context'])
    mut_sum_list = (df2['num_sites'])
    mu_gnomad_list = (df2['mu_gnomad'])
    
    meth_list = (df2['methylation_level'])
    theta_no_heterogeneity = (df2['theta_no_heterogeneity'])
    
    appro_theta = (df2['mean_theta'])
    sd_theta = (df2['sd_theta'])
    
    anc_read = anc_read_list[i]
    der_read = der_read_list[i]   
    mut_sum = mut_sum_list[i]
    mu_gnomad =  mu_gnomad_list[i]
    
    mu_corrected = theta_no_heterogeneity[i]
    meth_level = meth_list[i]
    
    theta_mean = appro_theta[i]
    theta_sd = sd_theta[i]
    
    comb = [anc_read, der_read, mut_sum, mu_gnomad,
            theta_mean, theta_sd,
            mu_corrected, meth_level
            ]
    
    mut_rate_list.append(comb)
    
#%%
common_list_1 =[]
sfs_common = []
for mr in mut_rate_list:
    for sfs in sfs_list:
        if mr[:3]==sfs[:3]:
            common_list_1.append(mr)
            sfs_common.append(sfs)
            
print(f'Common list 1 match: {len(common_list_1)}')
#%%

set1 = set(map(tuple, sfs_list))
set2 = set(map(tuple, mut_rate_list))

unique_in_sfs = set1 - set(map(tuple,sfs_common))

unique_in_sfs_list = [list(elem) for elem in unique_in_sfs]     

unique_in_mutr = set2 - set(map(tuple,common_list_1))

# Convert back to list of lists if needed
unique_in_mutr_list = [list(elem) for elem in unique_in_mutr]     

#%%

mapping = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}

# Function to transform a string based on the mapping
def transform_string(s):
    return ''.join(mapping[char] for char in s)

# Transform the first two elements of each inner list
mutr_list_conv = [
    [transform_string(item) if isinstance(item, str) else item for item in lst[:2]] + lst[2:]
    for lst in unique_in_mutr_list]

#%%
def swap_first_third(s):
    if len(s) < 3:
        return s  # Return the string as is if its length is less than 3
    return s[2] + s[1] + s[0]  # Swap the 1st and 3rd characters

# Apply the transformation
mutr_swap_list = [
    [swap_first_third(item) if isinstance(item, str) else item for item in sublist]
    for sublist in mutr_list_conv]

swap_match = []
sfs_swap = []
for sw in mutr_swap_list:
    for j in  unique_in_sfs_list:
        if sw[:3]==j[:3]:
            swap_match.append(sw)
            sfs_swap.append(j)

print(f'Swap match: {len(swap_match)}')
#%%

all_append = common_list_1 + swap_match 

print(f'Len of all: {len(all_append)}')

matched_df = pd.DataFrame(all_append, columns=["Anc DNA", "Der DNA", "Mut sum", "Mew",
                                               "Mean theta", "SD theta",
                                               "Mew corrected","Meth level"
                                               ])

#%%
df2_subset = matched_df.iloc[:, :3]

common_elements = [sublist2 for sublist1 in all_append for sublist2 in sfs_list if sublist1[:3] == sublist2]
common_count = len(common_elements)

print("Number of common elements:", common_count)

#%%
#matched_df.to_csv('D:/project2_yuval/matched_data_mews_144_2.csv', index=False)

#%%

new_df = pd.DataFrame()
for i in range(0,len(df),1001):
    anc_read_list = (df['ref_context'])
    der_read_list = (df['alt_context'])
    anc_read = anc_read_list[i]
    der_read = der_read_list[i]
    
    mut_sum = sum(df['n'][i:i+1001])
    
    comb = [anc_read, der_read, mut_sum]
    for ml in all_append:
        if ml[:3] == comb:
            df_x = df.iloc[i:i+1001]
            df_x['Mew'] = ml[3]
            
            df_x['Mew corrected'] = ml[6]
            df_x['Meth level'] = ml[7]
            
            df_x['Mean theta'] = ml[4]
            df_x['SD theta'] = ml[5]
            
            new_df = pd.concat([new_df, df_x], ignore_index=True)
            all_append.remove(ml)
    
#new_df.to_csv('D:/project2_yuval/human_sfs_mu_corrected.csv', index=False)

# df_non_zero = df[df['AC']!=0]
# df_zero = df[df['AC']==0]


# non_seg_list = [value for index, row in df_zero.iterrows() for value in [row['Mew']] * row['n']]
# seg_list = [value for index, row in df_non_zero.iterrows() for value in [row['Mew']] * row['n']]

# new_df_1 = new_df.iloc[::1001, :].reset_index(drop=True)
#new_df.to_csv('D:/project2_yuval/matched_gamma_mu_sd_meth.csv', index=False)
#%%
# data = list(new_df)
# df = pd.DataFrame(data, columns=["Index", "Repeat", "Value"])

# Generate the transformed DataFrame
new_df['Mutation type'] = (new_df.index // 1001) + 1
cols = ['Mutation type'] + [col for col in new_df.columns if col != 'Mutation type']
new_df = new_df[cols]

rows = []
k=1
for _, row in new_df.iterrows():
    rows.extend([[#row["ref_context"], row["alt_context"],
                  row['Mutation type'], row['Meth level'],
                  row["AC"], row["Mew"],
                  row["Mean theta"], row["SD theta"]]] * int(row["n"]))

# Create a new DataFrame with the repeated rows
df_y = pd.DataFrame(rows, columns=[#"ref_context", "alt_context",
               'Mutation type',   "Meth level",                  
              "AC", "Mew",
              "Mean theta", "SD theta"
              ])

df_y.to_csv('/project/yuvalsim/Deep/project2/gamma_data/human_sfs_ets_mus_sds.csv', index=False)
    
    
    
    
    
    
    