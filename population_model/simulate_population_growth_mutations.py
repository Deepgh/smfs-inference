#!/usr/bin/env python3

import math
import multiprocessing
import os
from random import randrange

import numpy as np
import pandas as pd

# Paths and files
INPUT_DIR = "/project/yuvalsim/Deep/project2/josh_full_data/error_data/"
INPUT_FILE = "scaled_mu_diff_site.csv"
OUTPUT_DIR = "/project/yuvalsim/Deep/project2/josh_full_data/error_data/"
OUTPUT_FILE = "der_allel_gamma_smp_mu_diff_site.csv"
MUTATION_RATE_OUTPUT_FILE = "mut_rates_gamma_smp_diff_site.csv"

# Column names
SITE_COL = "Site"
MUTATION_NUMBER_COL = "Mutation number"
MUTATION_RATE_COL = "Mutation rate"
MEAN_THETA_COL = "Mean theta"
SD_THETA_COL = "SD theta"
METHYLATION_LEVEL_COL = "Meth level"
INPUT_METHYLATION_LEVEL_COL = "Meth_level"
SAMPLED_MU_COL = "Sampled mu"
ALLELE_COPIES_COL = "Allele copies"
MUTATION_ID_COL = "Mutation id"

# Population model settings
INITIAL_POPULATION = 10**4
FINAL_POPULATION = 10**8
TIME_SCALE_MULTIPLIER = 250
NORMAL_APPROXIMATION_THRESHOLD = 10**9
RANDOM_SEED_MODULUS = 123456789123456
RANDOM_SEED_RANGE = 2 * 10**7

pd.options.mode.chained_assignment = None


def get_sample(dist: str, xi: float) -> float:
    return np.random.poisson(xi) if dist == "poisson" else np.random.normal(xi, np.sqrt(xi))


def extract_row_data(row):
    return {
        MUTATION_NUMBER_COL: row[MUTATION_NUMBER_COL],
        MUTATION_RATE_COL: row[MUTATION_RATE_COL],
        MEAN_THETA_COL: row[MEAN_THETA_COL],
        SD_THETA_COL: row[SD_THETA_COL],
        METHYLATION_LEVEL_COL: row[INPUT_METHYLATION_LEVEL_COL],
        SAMPLED_MU_COL: row[SAMPLED_MU_COL],
    }


def simulate_mutations_over_time(time_scale, xi, initial_population, sampled_mu, mutid_start):
    current_population = initial_population
    mutation_data = []
    results = []
    mutid = mutid_start

    for t in reversed(range(1, time_scale + 1)):
        theta = 2 * current_population * sampled_mu
        m = np.random.poisson(theta)

        for _ in range(m):
            mutation_data.append([1, t, "dist", mutid])
            mutid += 1

        for mut in mutation_data:
            rate = xi * mut[0]
            dist = "poisson" if rate < NORMAL_APPROXIMATION_THRESHOLD else "normal"
            mut[0] = get_sample(dist, rate)

        mutation_data = [m for m in mutation_data if m[0] != 0]

        if t == 1:
            for mut in mutation_data:
                results.append((mut[0], mut[3]))

        current_population *= xi

    return results


def single_iteration(args):
    iter_index, row = args
    np.random.seed((randrange(RANDOM_SEED_RANGE) + os.getpid()) % RANDOM_SEED_MODULUS)

    data = extract_row_data(row)
    time_scale = math.ceil(TIME_SCALE_MULTIPLIER * math.log10(FINAL_POPULATION / INITIAL_POPULATION))
    xi = (FINAL_POPULATION / INITIAL_POPULATION) ** (1 / time_scale)

    df_muts_rate = pd.DataFrame([{SITE_COL: iter_index, **data}])

    sim_results = simulate_mutations_over_time(
        time_scale=time_scale,
        xi=xi,
        initial_population=INITIAL_POPULATION,
        sampled_mu=data[SAMPLED_MU_COL],
        mutid_start=0,
    )

    final_data = [
        {
            SITE_COL: iter_index,
            **data,
            ALLELE_COPIES_COL: allele_copies,
            MUTATION_ID_COL: mutid,
        }
        for allele_copies, mutid in sim_results
    ]

    return pd.DataFrame(final_data), df_muts_rate


def move_column_first(df, column):
    if column not in df.columns:
        return df

    cols = [column] + [col for col in df.columns if col != column]
    return df[cols]


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    sampled_mu_df = pd.read_csv(os.path.join(INPUT_DIR, INPUT_FILE))
    run_items = list(sampled_mu_df.iterrows())
    num_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", multiprocessing.cpu_count()))

    with multiprocessing.Pool(num_cores) as pool:
        results = pool.map(single_iteration, run_items)

    df_results = pd.concat([r[0] for r in results], ignore_index=True).sort_values(by=SITE_COL)
    df_rates = pd.concat([r[1] for r in results], ignore_index=True).sort_values(by=SITE_COL)

    df_results = move_column_first(df_results, MUTATION_NUMBER_COL)
    df_rates = move_column_first(df_rates, MUTATION_NUMBER_COL)

    df_results.to_csv(os.path.join(OUTPUT_DIR, OUTPUT_FILE), index=False)
    df_rates.to_csv(os.path.join(OUTPUT_DIR, MUTATION_RATE_OUTPUT_FILE), index=False)


if __name__ == "__main__":
    main()
