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
INPUT_DIR = "/project/yuvalsim/Deep/project2/josh_full_data/error_data/simulation_results/"
INPUT_FILE = "edited_sim_data_smp_mu_diff_site.csv.gz"
UNIQUE_MUTATION_SITES_DIR = INPUT_DIR
UNIQUE_MUTATION_SITES_FILE = "unq_mut_sites_300_josh_mu.csv"
OUTPUT_DIR = UNIQUE_MUTATION_SITES_DIR

# Column names
AC_COL = "AC_"
SITES_COL = "Sites"
UNIQUE_MUTATION_SITES_COL = "Unq muts sites"

# Model settings
AC_LIMIT = 300
OUTPUT_FILE = f"smfs_num_{AC_LIMIT}_sim_smp_mu.csv"

# Optional error correction settings
APPLY_ERROR_CORRECTION = False
ERROR_FILE = "ajhg_00004094_supp_table2_mut.tsv"
ERROR_NUM_SITES_COL = "num_sites"
ERROR_RATE_COL = "error"
ERROR_SCALE = 10**6


def main():
    apply_config(CONFIG_SECTION, globals())

    data_df = pd.read_csv(
        os.path.join(INPUT_DIR, INPUT_FILE),
        compression="gzip",
    )

    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    pmf_csv = os.path.join(UNIQUE_MUTATION_SITES_DIR, UNIQUE_MUTATION_SITES_FILE)
    df_pmf = pd.read_csv(pmf_csv)

    eps_1 = df_pmf.iloc[1][UNIQUE_MUTATION_SITES_COL]
    num_cop_1 = sum(data_df[data_df[AC_COL] == 1][SITES_COL])
    prop_unq_muts_1 = num_cop_1 / eps_1

    if APPLY_ERROR_CORRECTION:
        err_df = pd.read_csv(os.path.join(UNIQUE_MUTATION_SITES_DIR, ERROR_FILE), sep="\t")
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
        print(acs + 1)

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

    new_df = pd.DataFrame()

    for i in range(len(f)):
        data = pd.DataFrame([{"Count": i + 1, "pred prob": f[i]}])
        new_df = pd.concat([new_df, data], ignore_index=True)

    new_df.to_csv(os.path.join(OUTPUT_DIR, OUTPUT_FILE), index=False)


if __name__ == "__main__":
    main()
