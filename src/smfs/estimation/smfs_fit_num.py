#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 00:12:49 2024

@author: ghosh1
"""

import numpy as np

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
    p_ac_given_1_mut = np.zeros(ac_limit)
    for ac in range(1, ac_limit+1): 
        p_ac_given_j_mut = p_ac_given_1_mut[:ac-1]  
        expected_ac_count_from_gt1_mut = 0
        for j in range(2, ac+1):
            p_ac_given_j_mut = np.convolve(p_ac_given_1_mut[:ac-1], p_ac_given_j_mut)
            expected_ac_count_from_gt1_mut += p_ac_given_j_mut[ac - j]*expected_sites_per_mutation_count[j]

        p_ac_given_1_mut[ac-1] = (
            (np.sum(allele_counts == ac) - expected_ac_count_from_gt1_mut)/expected_sites_per_mutation_count[1] 
        )
        
    return p_ac_given_1_mut

