#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 00:12:49 2024

@author: ghosh1
"""

import numpy as np

class BasicSMFSEstimator:
    def __init__(self, allele_counts, expected_sites_1_mut):
        self.counts = allele_counts
        self.expected_sites_1_mut = expected_sites_1_mut
    
    def get(self, ac, expected_ac_count_from_gt1_mut):
        return (np.sum(self.counts == ac) - expected_ac_count_from_gt1_mut)/self.expected_sites_1_mut
    

def marginal_p_k_given_j_gt1(k, p_k_given_j_eq1, p_j):
    p_k_given_j = p_k_given_j_eq1[:k-1]  
    p_k_given_j_gt1 = 0
    for j in range(2, k+1):
        p_k_given_j = np.convolve(p_k_given_j_eq1[:k-1], p_k_given_j)
        p_k_given_j_gt1 += p_k_given_j[k - j]*p_j[j]
    return p_k_given_j_gt1

def compute(allele_counts, expected_sites_per_mutation_count, ac_limit):
    """
    Estimate the conditional site mutation frequency spectrum (SMFS), i.e.,
    P(AC = k | 1 mutation), for allele counts ranging from 1 to `ac_limit`.

    The method assumes that each site may have received one or more mutations,
    and that the allele count (AC) at a site is the sum of independent contributions
    from each mutation. Under this assumption, the distribution of AC under j mutations
    is the j-fold convolution of the single-mutation AC distribution.

    This function uses observed allele counts and the expected number of sites with
    j unique mutations to iteratively deconvolve and recover the underlying
    distribution P(AC = k | 1 mutation) for each k = 1, ..., ac_limit.

    Parameters
    ----------
    allele_counts : np.ndarray
        Array of integer allele counts per site.

    expected_sites_per_mutation_count : np.ndarray
        Array where the entry at index j gives the expected number of
        sites with exactly j unique mutations. Must be defined for at least
        j = 1 to ac_limit.

    ac_limit : int
        Maximum allele count (inclusive) to compute SMFS up to. Must be smaller than or 
        equal to the max number of mutations used in `expected_sites_per_mutation_count`.

    Returns
    -------
    np.ndarray
        Array of length `ac_limit`, where the k-1-th entry gives
        P(AC = k | 1 mutation) for k = 1 to ac_limit.
    """
    sfms = BasicSMFSEstimator(allele_counts, expected_sites_per_mutation_count[1])
    p_ac_given_1_mut = np.zeros(ac_limit)
    for ac in range(1, ac_limit+1): 
        expected_ac_count_from_gt1_mut = marginal_p_k_given_j_gt1(ac, p_ac_given_1_mut, expected_sites_per_mutation_count)
        p_ac_given_1_mut[ac-1] = sfms.get(ac, expected_ac_count_from_gt1_mut)
        
    return p_ac_given_1_mut

