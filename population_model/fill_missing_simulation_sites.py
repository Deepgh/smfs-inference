#!/usr/bin/env python3

import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow_config import apply_config

CONFIG_SECTION = "population_model.fill_missing_simulation_sites"

# Paths and files
DERIVED_ALLELE_DIR = (
    "/project/yuvalsim/Deep/project2/josh_full_data/error_data/"
    "simulation_results/sampled_data/"
)
DERIVED_ALLELE_FILE = "der_allele_sum_per_site_smp_mu_diff_site.csv.gz"
MUTATION_RATE_DIR = "/project/yuvalsim/Deep/project2/josh_full_data/error_data/simulation_results/"
MUTATION_RATE_FILE = "mut_rates_gamma_smp_diff_site.csv"
OUTPUT_DIR = MUTATION_RATE_DIR
OUTPUT_FILE = "edited_sim_data_smp_mu_diff_site.csv.gz"

# Column names
SITE_COL = "Site"
SAMPLED_MUTATION_SUM_COL = "Sampled mutation sum"
UNIQUE_MUTATIONS_COL = "Unique mutations"
AC_COL = "AC"


def main():
    apply_config(CONFIG_SECTION, globals())

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    derived_allele_df = pd.read_csv(
        os.path.join(DERIVED_ALLELE_DIR, DERIVED_ALLELE_FILE),
        compression="gzip",
    )
    mutation_rate_df = pd.read_csv(os.path.join(MUTATION_RATE_DIR, MUTATION_RATE_FILE))

    merge_cols = [SITE_COL, SAMPLED_MUTATION_SUM_COL, UNIQUE_MUTATIONS_COL]
    merged_df = mutation_rate_df.merge(
        derived_allele_df[merge_cols],
        on=SITE_COL,
        how="left",
    )

    fill_cols = [SAMPLED_MUTATION_SUM_COL, UNIQUE_MUTATIONS_COL]
    merged_df[fill_cols] = merged_df[fill_cols].fillna(0).astype(int)

    merged_df = merged_df.rename(columns={SAMPLED_MUTATION_SUM_COL: AC_COL})
    merged_df[SITE_COL] += 1

    merged_df.to_csv(
        os.path.join(OUTPUT_DIR, OUTPUT_FILE),
        compression="gzip",
        index=False,
    )


if __name__ == "__main__":
    main()
