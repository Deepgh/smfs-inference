#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 18:15:18 2025

Author: ghosh1
"""

import numpy as np
import pandas as pd
import multiprocessing
import os
import math
from random import randrange

# constants
initial_pop =  10**4
final_pop =  10**8

# file paths
data_path = '/project/yuvalsim/Deep/project2/josh_full_data/error_data/'
muts_file = 'scaled_mu_diff_site.csv'
scaled_mu_file = os.path.join(data_path, muts_file)

output_dir = '/project/yuvalsim/Deep/project2/josh_full_data/error_data/'
output_file = 'der_allel_gamma_smp_mu_diff_site.csv'
mut_rate_file = 'mut_rates_gamma_smp_diff_site.csv'

# read mutation data
sampled_mu_df = pd.read_csv(scaled_mu_file)
pd.options.mode.chained_assignment = None


def get_sample(dist: str, lambdaa: float) -> float:
    return np.random.poisson(lambdaa) if dist == 'poisson' else np.random.normal(lambdaa, np.sqrt(lambdaa))


def extract_row_data(row):
    return {'Mutation number': row['Mutation number'],
        'Mutation rate': row['Mutation rate'],
        'Mean theta': row['Mean theta'],
        'SD theta': row['SD theta'],
        'Meth level': row['Meth_level'],
        'Sampled mu': row['Sampled mu'] }

def simulate_mutations_over_time(time_scale, lambdaa, initial_population, smp_mew, mutid_start):
    current_population = initial_population
    mutation_data = []
    results = []
    mutid = mutid_start

    for t in reversed(range(1, time_scale + 1)):
        theta = 2 * current_population * smp_mew
        m = np.random.poisson(theta)

        for _ in range(m):
            mutation_data.append([1, t, 'dist', mutid])
            mutid += 1

        for mut in mutation_data:
            rate = lambdaa * mut[0]
            dist = 'poisson' if rate < 10**9 else 'normal'
            mut[0] = get_sample(dist, rate)

        mutation_data = [m for m in mutation_data if m[0] != 0]

        if t == 1:
            for mut in mutation_data:
                results.append((mut[0], mut[3]))  # (allele copies, mutid)

        current_population *= lambdaa

    return results


def single_iteration(iter_index: int):
    
    np.random.seed((randrange(2 * 10**7) + os.getpid()) % 123456789123456)


    row = sampled_mu_df.iloc[iter_index]
    data = extract_row_data(row)

    time_scale = math.ceil(250 * math.log10(final_pop / initial_pop))
    lambdaa = (final_pop / initial_pop) ** (1 / time_scale)

    df_muts_rate = pd.DataFrame([{
        'Site': iter_index,
        **data}])

    sim_results = simulate_mutations_over_time(
        time_scale=time_scale,
        lambdaa=lambdaa,
        initial_population=initial_pop,
        smp_mew=data['Sampled mu'],
        mutid_start=0
    )

    final_data = [{
        'Site': iter_index,
        **data,
        'Allele copies': allele_copies,
        'Mutation id': mutid
    } for allele_copies, mutid in sim_results]

    return pd.DataFrame(final_data), df_muts_rate


def main():
    os.makedirs(output_dir, exist_ok=True)

    num_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', multiprocessing.cpu_count()))
    run_range = len(sampled_mu_df)

    with multiprocessing.Pool(num_cores) as pool:
        results = pool.map(single_iteration, range(run_range))

    df_results = pd.concat([r[0] for r in results], ignore_index=True).sort_values(by='Site')
    df_rates = pd.concat([r[1] for r in results], ignore_index=True).sort_values(by='Site')
    
    if 'Mutation number' in df_results.columns:
        cols = ['Mutation number'] + [col for col in df_results.columns if col != 'Mutation number']
        df_results = df_results[cols]
    if 'Mutation number' in df_results.columns:
        cols = ['Mutation number'] + [col for col in df_rates.columns if col != 'Mutation number']
        df_rates = df_rates[cols]
        
    df_results.to_csv(os.path.join(output_dir, output_file), index=False)
    df_rates.to_csv(os.path.join(output_dir, mut_rate_file), index=False)


if __name__ == "__main__":
    main()
