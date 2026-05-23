#!/usr/bin/env python3

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow_config import apply_config

CONFIG_SECTION = "smfs.estimate_smfs_from_unique_mutations"

# Paths and files
INPUT_DIR = "data"
INPUT_FILE = "observed_allele_counts.csv.gz"
UNIQUE_MUTATION_SITES_DIR = "results"
UNIQUE_MUTATION_SITES_FILE = "expected_unique_mutation_sites_negative_binomial_300.csv"
OUTPUT_DIR = "results"

# Column names
AC_COL = "allele_count"
SITES_COL = "sites"
UNIQUE_MUTATION_SITES_COL = "unique_mutation_sites"

# Model settings
AC_LIMIT = 300
OUTPUT_FILE = f"smfs_observed_{AC_LIMIT}.csv"

# Optional error correction settings
APPLY_ERROR_CORRECTION = False
ERROR_FILE = "data/mutation_rates.tsv"
ERROR_NUM_SITES_COL = "num_sites"
ERROR_RATE_COL = "error"
ERROR_SCALE = 10**6


def main(data_df, df_pmf, err_df=None):
    eps_1 = df_pmf.iloc[1][UNIQUE_MUTATION_SITES_COL]
    num_cop_1 = sum(data_df[data_df[AC_COL] == 1][SITES_COL])
    prop_unq_muts_1 = num_cop_1 / eps_1

    if APPLY_ERROR_CORRECTION:
        if err_df is None:
            raise ValueError("err_df is required when APPLY_ERROR_CORRECTION is true")

        wt_sc_error = sum(err_df[ERROR_NUM_SITES_COL] * err_df[ERROR_RATE_COL] * ERROR_SCALE)
        prop_unq_muts_1 = (num_cop_1 - wt_sc_error) / eps_1

    con_list = [prop_unq_muts_1]

    f_2 = np.convolve(con_list, con_list)
    con_list.append(f_2[0])
    f2 = (
        sum(data_df[data_df[AC_COL] == 2][SITES_COL])
        - (f_2[0] * df_pmf.iloc[2][UNIQUE_MUTATION_SITES_COL])
    ) / eps_1

    f = [prop_unq_muts_1, f2]

    for acs in range(2, AC_LIMIT):
        conv_lists = []
        f_convolved = f

        for _ in range(acs):
            f_convolved = np.convolve(f, f_convolved)
            conv_lists.append(f_convolved)

        s1 = 0
        l = len(conv_lists)
        j = acs + 1
        k = 0

        while l >= 1:
            s1 = s1 + conv_lists[l - 1][k] * df_pmf.iloc[j][UNIQUE_MUTATION_SITES_COL]
            j -= 1
            k += 1
            l -= 1

        val = (sum(data_df[data_df[AC_COL] == acs + 1][SITES_COL]) - s1) / eps_1
        f.append(val)

    return pd.DataFrame(
        [{"allele_count": i + 1, "predicted_probability": value} for i, value in enumerate(f)]
    )


if __name__ == "__main__":
    apply_config(CONFIG_SECTION, globals())

    input_df = pd.read_csv(
        os.path.join(INPUT_DIR, INPUT_FILE),
        compression="gzip",
    )

    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    pmf_csv = os.path.join(UNIQUE_MUTATION_SITES_DIR, UNIQUE_MUTATION_SITES_FILE)
    pmf_df = pd.read_csv(pmf_csv)

    error_df = None
    if APPLY_ERROR_CORRECTION:
        error_df = pd.read_csv(ERROR_FILE, sep="\t")

    output_df = main(input_df, pmf_df, error_df)
    output_df.to_csv(os.path.join(OUTPUT_DIR, OUTPUT_FILE), index=False)
