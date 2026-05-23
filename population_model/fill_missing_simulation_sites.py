#!/usr/bin/env python3

import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow_config import apply_config

CONFIG_SECTION = "population_model.fill_missing_simulation_sites"

# Paths and files
DERIVED_ALLELE_DIR = "results/sampled_data"
DERIVED_ALLELE_FILE = "observed_mutations_by_site.csv.gz"
MUTATION_RATE_DIR = "results"
MUTATION_RATE_FILE = "simulated_site_mutation_rates.csv"
OUTPUT_DIR = MUTATION_RATE_DIR
OUTPUT_FILE = "simulation_sites_with_observed_counts.csv.gz"

# Column names
SITE_COL = "site"
SAMPLED_MUTATION_SUM_COL = "sampled_mutation_sum"
UNIQUE_MUTATIONS_COL = "unique_mutations"
AC_COL = "allele_count"


def main(derived_allele_df, mutation_rate_df):
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

    return merged_df


if __name__ == "__main__":
    apply_config(CONFIG_SECTION, globals())

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    derived_allele_input_df = pd.read_csv(
        os.path.join(DERIVED_ALLELE_DIR, DERIVED_ALLELE_FILE),
        compression="gzip",
    )
    mutation_rate_input_df = pd.read_csv(os.path.join(MUTATION_RATE_DIR, MUTATION_RATE_FILE))

    output_df = main(derived_allele_input_df, mutation_rate_input_df)
    output_df.to_csv(
        os.path.join(OUTPUT_DIR, OUTPUT_FILE),
        compression="gzip",
        index=False,
    )
