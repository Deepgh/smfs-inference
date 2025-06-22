#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 16:43:15 2024

@author: ghosh1
"""

import numpy as np
from scipy.optimize import minimize
from scipy.stats import poisson, nbinom, bernoulli
from typing import Callable, Sequence, Union

def poisson_pmf_model(mutation_rates: Union[Sequence[float], np.ndarray]):
    return lambda k, lam: poisson.pmf(k, lam * mutation_rates)

def poisson_gamma_pmf_model(mu, sigma):
    return lambda k, lam: nbinom.pmf(k, n=(mu/sigma)**2, p=1/(1 + sigma**2 *(lam/mu)))

def fit_model_scaling(
    p0_func: Callable[[float], np.ndarray],
    is_seg: np.ndarray,
    initial_lambda: float = 1e8
) -> float:
    """
    Estimate the global mutation rate scaling factor λ using binary segregation status.

    This function estimates λ by maximizing the likelihood of observing
    segregating and non-segregating sites under a given probabilistic model
    of P(K = 0 | λ), where K is the number of observed mutations at a site.

    Parameters
    ----------
    p0_func : Callable[[float], np.ndarray]
        A function that takes a scalar λ and returns a vector of
        probabilities P(K = 0 | λ) for each site.
    is_seg : np.ndarray
        A boolean array where True indicates that the site is segregating
        (i.e., K > 0), and False indicates that the site is non-segregating (K = 0).
    initial_lambda : float, optional
        Initial guess for λ (default is 1e8). Optimization is performed in log space.

    Returns
    -------
    float
        The maximum likelihood estimate of λ.
    """
    def objective(log_lam: float) -> float:
        lam = np.exp(log_lam)
        p0 = p0_func(lam)
        return -bernoulli.logpmf(~is_seg, p0).sum()
    
    result = minimize(
        objective,
        np.log(initial_lambda),
        method='nelder-mead',
        options={'xatol': 1e-8, 'disp': True}
    )
    return np.exp(result.x[0])


def unique_mutation_distribution(
        site_mutation_probability: Callable[[int], np.ndarray],
        max_mutations: int,
        stop_on_zero: bool = False
) -> np.ndarray:
    result = np.zeros(max_mutations + 1, dtype=float)
    for k in range(max_mutations + 1):
        result[k] = site_mutation_probability(k).mean() # P(# unique mutations = k)
        if stop_on_zero and result[k] == 0.0: break
    return result


def expected_sites_per_mutation_count(
    site_mutation_probability: Callable[[int], np.ndarray],
    max_mutations: int, 
    stop_on_zero: bool = True
) -> np.ndarray:
    """
    Compute the expected number of sites with k unique mutations for k = 0..max_mutations.

    Parameters
    ----------
    site_mutation_probability : Callable[[int], np.ndarray]
        A function that takes an integer k and returns a NumPy array of probabilities,
        such that the i-th value is "the probability that site i has exactly k unique mutations".

    max_mutations : int
        The maximum number of unique mutations (k) to evaluate.

    stop_on_zero : bool, optional
        If True, stop computing when the total expected number of sites with k mutations becomes zero.

    Returns
    -------
    np.ndarray
        An array of length (max_mutations + 1) containing the expected number of sites
        with exactly k unique mutations for k = 0 to max_mutations.
    """

    result = np.zeros(max_mutations + 1, dtype=float)
    for k in range(max_mutations + 1):
        result[k] = site_mutation_probability(k).sum()
        if stop_on_zero and result[k] == 0.0: break
    return result



















