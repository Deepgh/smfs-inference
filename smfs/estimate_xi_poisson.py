#!/usr/bin/env python3

import math
import os
import sys
from pathlib import Path
from warnings import simplefilter

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import poisson

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow_config import apply_config

CONFIG_SECTION = "smfs.estimate_xi_poisson"

# Paths and files
INPUT_DIR = "/project/yuvalsim/Deep/project2/josh_full_data/error_data/simulation_results/"
INPUT_FILE = "edited_sim_data_smp_mu_diff_site.csv.gz"
OUTPUT_DIR = INPUT_DIR
OUTPUT_XI_FILE = "lam_smp_mu.csv"

# Column names
SITE_COL = "Site"
MU_COL = "Sampled mu"
AC_COL = "AC"

# Model settings
INITIAL_XI = 1e3
UNIQUE_MUTATION_LIMIT = 300
OUTPUT_UNIQUE_MUTATION_FILE = f"unq_mut_sites_{UNIQUE_MUTATION_LIMIT}_sim_smp_mu.csv"
OPTIMIZER_METHOD = "nelder-mead"
OPTIMIZER_OPTIONS = {"xatol": 1e-8, "disp": True}
EPS = 1e-12

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def estimate(log_xi, mu_seg, mu_non_seg):
    """Negative log-likelihood for the Poisson mutation-rate model."""
    xi_val = np.exp(log_xi)

    p_seg = np.exp(-xi_val * mu_seg)
    p_seg = np.clip(p_seg, EPS, 1 - EPS)

    p_non = np.exp(-xi_val * mu_non_seg)
    p_non = np.clip(p_non, EPS, 1 - EPS)

    log_lik = np.sum(np.log(p_non)) + np.sum(np.log(1 - p_seg))
    return -log_lik


def main(data_df):
    data_df = data_df.copy()

    # If rows contain site counts rather than one row per site, expand here.
    # data_df = data_df.loc[
    #     data_df.index.repeat(data_df[SITE_COL]), [MU_COL, AC_COL]
    # ].reset_index(drop=True)

    mu_seg = np.array(data_df[data_df[AC_COL] != 0][MU_COL])
    mu_non_seg = np.array(data_df[data_df[AC_COL] == 0][MU_COL])

    res = minimize(
        estimate,
        np.log([INITIAL_XI]),
        method=OPTIMIZER_METHOD,
        args=(mu_seg, mu_non_seg),
        options={**OPTIMIZER_OPTIONS, "disp": False},
    )

    xi_est = float(np.exp(res.x[0]))

    xi_df = pd.DataFrame([{"Lambda": xi_est}])

    xi = float(math.exp(res.x[0]))
    unique_mutation_rows = []
    stop_filling = False

    for unq_mut in range(UNIQUE_MUTATION_LIMIT + 1):
        if not stop_filling:
            col_name = f"Unq muts sites_{unq_mut}"
            data_df[col_name] = poisson.pmf(unq_mut, xi * data_df[MU_COL])
            poi_sum = sum(data_df[col_name])

            if poi_sum == 0:
                stop_filling = True
        else:
            poi_sum = 0

        unique_mutation_rows.append(
            {"Unique mutation": unq_mut, "Unq muts sites": poi_sum}
        )

    unique_mutation_df = pd.DataFrame(unique_mutation_rows)
    return xi_est, xi_df, unique_mutation_df


if __name__ == "__main__":
    apply_config(CONFIG_SECTION, globals())

    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    input_df = pd.read_csv(
        os.path.join(INPUT_DIR, INPUT_FILE),
        compression="gzip",
    )

    estimated_xi, xi_output_df, unique_mutation_output_df = main(input_df)

    print("Estimated xi:", estimated_xi)

    xi_output_df.to_csv(os.path.join(OUTPUT_DIR, OUTPUT_XI_FILE), index=False)
    unique_mutation_output_df.to_csv(
        os.path.join(OUTPUT_DIR, OUTPUT_UNIQUE_MUTATION_FILE),
        index=False,
    )
