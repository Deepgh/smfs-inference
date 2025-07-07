#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 16:21:56 2025

@author: ghosh1
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from collections import defaultdict


class SMFSEstimator:
    def __init__(self, data_df, sites_df, ac_limit, num_mut_types):
        
        self.z_init = np.array([0.0])
        self.i_count, self.noni_count = self.get_ac_counts_by_type(data_df, ac_limit, num_mut_types)

        self.p_j_eq1 = np.array(sites_df[sites_df['Unique mutation'] == 1].values.tolist()[0][1:])
        self.p_k_geq_ac_and_j_neq0 = 1 - np.array(sites_df[sites_df['Unique mutation'] == 0].values.tolist()[0][1:])

    def get_ac_counts_by_type(self, data_df, ac_limit, num_mut_types):
        i_count = defaultdict(lambda: np.zeros(num_mut_types, dtype=int))
        noni_count = defaultdict(lambda: np.zeros(num_mut_types, dtype=int))
        for ac in range(ac_limit):
            for mut_typ in range(1, num_mut_types + 1):
                mut_typ_df = data_df[data_df['Mutation number'] == mut_typ]
                if len(mut_typ_df) != 0:
                    i_count[ac][mut_typ-1] = len(mut_typ_df[mut_typ_df['AC'] == ac])
                    noni_count[ac][mut_typ-1] = len(mut_typ_df[mut_typ_df['AC'] > ac])
        return i_count, noni_count

    @staticmethod
    def log_it_bounded(z, vals):
        return (1.0 - np.sum(vals))/(1 + np.exp(-z))
    

    def fit_tail(self, z, ac, prob_ac_one_unq_mut, p_k_eq_ac_and_j_gt1):
        """
        Compute the negative log-likelihood of observing allele count `ac` conditioned on the
        count being at least `ac` and having at least one mutation (j ≠ 0).

        This method models the probability of observing count K = ac as a mixture 
        of two components:
        - Single-mutation events: P(K = ac | j = 1) * P(j = 1)
        - Recurrent-mutation events: P(K = ac, j > 1)

        The likelihood is conditioned on the event (K ≥ ac, j ≠ 0), i.e., the site is 
        segregating and has at least one mutation.

        Parameters
        ----------
        z : float
            Logit-transformed scalar controlling the mixing weight of the j = 1 component.

        ac : int
            The allele count (K = ac) for which the likelihood is evaluated.

        prob_ac_one_unq_mut : array-like
            Distribution P(K = k | j = 1), used to compute the expected single-mutation
            contribution to allele count ac.

        p_k_eq_ac_and_j_gt1 : array-like
            Vector of recurrent contributions to P(K = ac, j > 1), typically computed
            via convolution.

        Returns
        -------
        float
            Negative log-likelihood of the data at allele count ac under the current model.
        """
        p_k_eq_ac_and_j_eq1 = self.log_it_bounded(z, prob_ac_one_unq_mut) * self.p_j_eq1
        p_k_eq_ac_given_tail = (
            p_k_eq_ac_and_j_eq1 + p_k_eq_ac_and_j_gt1
        ) / self.p_k_geq_ac_and_j_neq0

        return -(
            np.sum(self.i_count[ac] * np.log(p_k_eq_ac_given_tail)) +
            np.sum(self.noni_count[ac] * np.log(1 - p_k_eq_ac_given_tail))
        )
    
    def get(self, ac, p_ac_given_1_mut, p_k_eq_ac_and_j_gt1):
        res = minimize(
            self.fit_tail, 
            self.z_init, 
            method='nelder-mead',
            args=(ac, p_ac_given_1_mut[:ac], p_k_eq_ac_and_j_gt1),
            options={'xatol': 1e-8, 'disp': True}
        )
        p_k_eq_ac_given_j_eq1 = self.log_it_bounded(res.x[0], p_ac_given_1_mut[:ac])
        self.p_k_geq_ac_and_j_neq0 -= p_k_eq_ac_and_j_gt1 + p_k_eq_ac_given_j_eq1 * self.p_j_eq1
        return p_k_eq_ac_given_j_eq1


def marginal_p_k_given_j_gt1(k, p_k_given_j_eq1, p_j):
    """
    Compute the marginal probability P(K = k, 1 < j <= k).

    Parameters
    ----------
    k : int
        Target value of the total count K.

    p_k_given_j_eq1 : array-like
        Probability mass function for K given N = 1, i.e., P(K = i | N = 1) is given for i = 0 to k-1.

    p_j : np.ndarray
        Entry p_j[j] gives P(N = j) for 0 <= j <= k.

    Returns
    -------
    np.ndarray
        Marginal probabilities P(K = k, 1 < j <= k) = SUM P(K = k | N = j) * P(N = j).
        This is the contribution to K = k from j > 1 sources.
    """
    p_k_given_j = p_k_given_j_eq1[:k]  
    p_k_given_j_gt1 = np.zeros(p_j.shape[1], dtype=float)
    for j in range(2, k+1):
        p_k_given_j = np.convolve(p_k_given_j_eq1[:k], p_k_given_j)
        p_k_given_j_gt1 += p_k_given_j[k - j] * p_j[j]
    return p_k_given_j_gt1

def get_mut_count_pmf_by_type(sites_df, ac_limit, num_mut_types):

    mut_count_pmf = np.zeros((ac_limit+1, num_mut_types))
    for j in range(ac_limit + 1):
        mut_count_pmf[j] = np.array(
            sites_df[sites_df['Unique mutation'] == j].values.tolist()[0][1:]
        )
    return mut_count_pmf

def compute_sfms(data_df, sites_df, ac_limit):
    num_mut_types = data_df['Mutation number'].nunique()

    smfs = SMFSEstimator(data_df, sites_df, ac_limit, num_mut_types)
    # j is unique mutation count
    # k is allele count
    mut_count_pmf = get_mut_count_pmf_by_type(sites_df, ac_limit, num_mut_types)

    p_ac_given_1_mut = np.zeros(ac_limit)
    for ac in range(ac_limit): 
        p_ac_given_mut_count_gt1 = marginal_p_k_given_j_gt1(ac, p_ac_given_1_mut, mut_count_pmf)
        p_ac_given_1_mut[ac] = smfs.get(ac, p_ac_given_1_mut, p_ac_given_mut_count_gt1)
    return p_ac_given_1_mut
    
if __name__ == "__main__":
    
    import os
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

    prob_ac_one_unq_mut = compute_sfms(data_df, sites_df, ac_limit)

    new_df = pd.DataFrame()
    for i in range(len(prob_ac_one_unq_mut)):
        data = pd.DataFrame([{'Count': i+1, 'pred prob':prob_ac_one_unq_mut[i]}])
        new_df = pd.concat([new_df, data], ignore_index=True)
    new_df.to_csv(os.path.join(output_path, output_filename), index=False)



