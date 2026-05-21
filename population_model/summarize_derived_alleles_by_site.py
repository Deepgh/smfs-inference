#!/usr/bin/env python3

import multiprocessing
import os
import sys
import warnings
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow_config import apply_config

CONFIG_SECTION = "population_model.summarize_derived_alleles_by_site"

# Paths and files
INPUT_DIR = "/project/yuvalsim/Deep/project2/josh_full_data/error_data/simulation_results/sampled_data/"
INPUT_FILE = "samp_data_10_5_5_smp_mu_diff_site.csv.gz"
OUTPUT_DIR = INPUT_DIR
OUTPUT_FILE = "der_allele_sum_per_site_smp_mu_diff_site.csv.gz"

# Column names
SITE_COL = "Site"
MUTATION_ID_COL = "Mutation id"
METHYLATION_LEVEL_COL = "Meth level"
MUTATION_RATE_COL = "Mutation rate"
MUTATION_NUMBER_COL = "Mutation number"
MEAN_THETA_COL = "Mean theta"
SD_THETA_COL = "SD theta"
SAMPLED_MU_COL = "Sampled mu"
SAMPLED_MUTATION_COL = "Sampled mutation"

warnings.filterwarnings("ignore")


def single_iteration(group):
    site = group[SITE_COL].iloc[0]
    unq_mut = len(group[MUTATION_ID_COL].unique())
    methylation_level = group[METHYLATION_LEVEL_COL].iloc[0]
    mut_rate = group[MUTATION_RATE_COL].iloc[0]
    mut_num = group[MUTATION_NUMBER_COL].iloc[0]
    mean_theta = group[MEAN_THETA_COL].iloc[0]
    sd_gamma = group[SD_THETA_COL].iloc[0]
    sampled_mu = group[SAMPLED_MU_COL].iloc[0]
    mut_sum = group[SAMPLED_MUTATION_COL].sum()

    return {
        "Mutation number": mut_num,
        "Site": site,
        "Mutation rate": mut_rate,
        "Mew theta": mean_theta,
        "SD theta": sd_gamma,
        "Sampled mu": sampled_mu,
        "Meth level": methylation_level,
        "Sampled mutation sum": mut_sum,
        "Unique mutations": unq_mut,
    }


def main():
    apply_config(CONFIG_SECTION, globals())

    df = pd.read_csv(os.path.join(INPUT_DIR, INPUT_FILE), compression="gzip")
    number_of_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", multiprocessing.cpu_count()))
    grouped_df = df.groupby(SITE_COL)

    with multiprocessing.Pool(number_of_cores) as pool:
        results = pool.map(single_iteration, [group for _, group in grouped_df])

    df_sum_sites = pd.DataFrame(results)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df_sum_sites.to_csv(
        os.path.join(OUTPUT_DIR, OUTPUT_FILE),
        compression="gzip",
        index=False,
    )


if __name__ == "__main__":
    main()
