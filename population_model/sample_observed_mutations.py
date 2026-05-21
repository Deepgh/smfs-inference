#!/usr/bin/env python3

import os

import pandas as pd
from scipy.stats import hypergeom

# Paths and files
INPUT_DIR = "/project/yuvalsim/Deep/project2/josh_full_data/error_data/"
INPUT_FILE = "der_allel_gamma_smp_mu_diff_site.csv"
OUTPUT_DIR = os.path.join(INPUT_DIR, "sampled_data")
OUTPUT_FILE = "samp_data_10_5_5_smp_mu_diff_site.csv.gz"

# Column names
ALLELE_COPIES_COL = "Allele copies"
SAMPLED_MUTATION_COL = "Sampled mutation"

# Sampling settings
FINAL_POPULATION = 2 * (10**8)
SAMPLE_SIZE = 2 * 5 * (10**5)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(os.path.join(INPUT_DIR, INPUT_FILE))
    df[SAMPLED_MUTATION_COL] = hypergeom.rvs(
        FINAL_POPULATION,
        df[ALLELE_COPIES_COL],
        SAMPLE_SIZE,
    )

    df_sampled = df[df[SAMPLED_MUTATION_COL] != 0]
    df_sampled.to_csv(
        os.path.join(OUTPUT_DIR, OUTPUT_FILE),
        compression="gzip",
        index=False,
    )


if __name__ == "__main__":
    main()
