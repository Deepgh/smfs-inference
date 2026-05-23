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
INPUT_DIR = "results/sampled_data"
INPUT_FILE = "sampled_observed_mutations.csv.gz"
OUTPUT_DIR = INPUT_DIR
OUTPUT_FILE = "observed_mutations_by_site.csv.gz"

# Column names
SITE_COL = "site"
MUTATION_ID_COL = "mutation_id"
METHYLATION_LEVEL_COL = "methylation_level"
MUTATION_RATE_COL = "mutation_rate"
MUTATION_NUMBER_COL = "mutation_number"
MEAN_THETA_COL = "mean_theta"
SD_THETA_COL = "sd_theta"
SAMPLED_MU_COL = "sampled_mu"
SAMPLED_MUTATION_COL = "sampled_mutation"

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
        "mutation_number": mut_num,
        "site": site,
        "mutation_rate": mut_rate,
        "mean_theta": mean_theta,
        "sd_theta": sd_gamma,
        "sampled_mu": sampled_mu,
        "methylation_level": methylation_level,
        "sampled_mutation_sum": mut_sum,
        "unique_mutations": unq_mut,
    }


def main(df, number_of_cores):
    grouped_df = df.groupby(SITE_COL)

    with multiprocessing.Pool(number_of_cores) as pool:
        results = pool.map(single_iteration, [group for _, group in grouped_df])

    return pd.DataFrame(results)


if __name__ == "__main__":
    apply_config(CONFIG_SECTION, globals())

    input_df = pd.read_csv(os.path.join(INPUT_DIR, INPUT_FILE), compression="gzip")
    cores = int(os.environ.get("SLURM_CPUS_PER_TASK", multiprocessing.cpu_count()))

    df_sum_sites = main(input_df, cores)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df_sum_sites.to_csv(
        os.path.join(OUTPUT_DIR, OUTPUT_FILE),
        compression="gzip",
        index=False,
    )
