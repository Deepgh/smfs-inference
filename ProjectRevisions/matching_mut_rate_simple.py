#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  6 23:29:20 2025

@author: ghosh1
"""

import pandas as pd


df = pd.read_csv('/project/yuvalsim/Deep/project2/josh_full_data/sfs_syn_ac_1M.csv.gz', compression='gzip')
df['methylation_level'] = df['methylation_level'].fillna(0).astype(int)

# Aggregate mutation count per unique triplet
sfs_grouped = df.groupby(['ref_context', 'alt_context', 'methylation_level'])['Sites'].sum().reset_index()
sfs_grouped.columns = ['ref_context', 'alt_context', 'methylation_level', 'mut_sum']

# Load second file
df2 = pd.read_csv('/project/yuvalsim/Deep/project2/josh_full_data/error_data/ajhg_00004094_supp_table2_mut.tsv',
                  sep='\t')

# Keep only required columns
df2 = df2[['ref_context', 'alt_context', 'methylation_level','mu_gnomad', 
           'mean_theta', 'sd_theta',  'error']]

# Merge on triplets
matched = pd.merge(df2, sfs_grouped, on=['ref_context', 'alt_context', 'methylation_level'], how='inner')
print(f"Common list 1 match: {len(matched)}")

# Identify unmatched
sfs_set = set(tuple(x) for x in sfs_grouped[['ref_context', 'alt_context', 'methylation_level']].values)
df2_set = set(tuple(x) for x in df2[['ref_context', 'alt_context', 'methylation_level']].values)
unmatched_in_sfs = sfs_set - df2_set
unmatched_in_df2 = df2_set - sfs_set

# # Map complementary bases
# mapping = str.maketrans("ACGT", "TGCA")

def complement(seq):
    return ''.join({'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}[base] for base in seq)

def swap_first_third(seq):
    return seq[2] + seq[1] + seq[0] if len(seq) == 3 else seq

# Convert df2 unmatched for swap testingon 
df2_unmatched = df2[df2.apply(lambda row: (row['ref_context'], row['alt_context'], row['methylation_level']) in unmatched_in_df2, axis=1)].copy()
df2_unmatched['ref_context'] = df2_unmatched['ref_context'].apply(complement)
df2_unmatched['alt_context'] = df2_unmatched['alt_context'].apply(complement)
df2_unmatched['ref_context'] = df2_unmatched['ref_context'].apply(swap_first_third)
df2_unmatched['alt_context'] = df2_unmatched['alt_context'].apply(swap_first_third)

# Match again
swap_match = pd.merge(df2_unmatched, sfs_grouped, on=['ref_context', 'alt_context', 'methylation_level'], how='inner')
print(f"Swap match: {len(swap_match)}")

# Combine all matches
all_matches = pd.concat([matched, swap_match], ignore_index=True)
print(f"Len of all: {len(all_matches)}")

triplet_key = ['ref_context', 'alt_context', 'methylation_level']
df = df.merge(all_matches, on=triplet_key, how='inner')

# Rename for clarity
df.rename(columns={
    'methylation_level': 'Meth level',
    'mu_gnomad': 'Mew',
    'mean_theta': 'Mean theta',
    'sd_theta': 'SD theta',
    #'theta_error': 'Error theta',
    'error': "Error"
}, inplace=True)

df = df.drop(columns=['mut_sum']) 

#%%
# Save to compressed file
df2 = pd.read_csv('/project/yuvalsim/Deep/project2/josh_full_data/error_data/sfs_syn_err_ac_1M_freq_dist.csv.gz', compression='gzip')
#df_y.to_csv('/project/yuvalsim/Deep/project2/josh_full_data/sfs_per_site_missense_ac_1M.csv.gz', compression='gzip', index=False)
