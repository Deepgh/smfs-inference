#!/usr/bin/env python3

import math
import os
from pathlib import Path
from warnings import simplefilter

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import poisson

# Paths and files
INPUT_DIR = "/project/yuvalsim/Deep/project2/josh_full_data/error_data/simulation_results/"
INPUT_FILE = "edited_sim_data_smp_mu_diff_site.csv.gz"
OUTPUT_DIR = INPUT_DIR
OUTPUT_LAMBDA_FILE = "lam_smp_mu.csv"

# Column names
SITE_COL = "Site"
MU_COL = "Sampled mu"
AC_COL = "AC"

# Model settings
INITIAL_LAMBDA = 1e3
UNIQUE_MUTATION_LIMIT = 300
OUTPUT_UNIQUE_MUTATION_FILE = f"unq_mut_sites_{UNIQUE_MUTATION_LIMIT}_sim_smp_mu.csv"
OPTIMIZER_METHOD = "nelder-mead"
OPTIMIZER_OPTIONS = {"xatol": 1e-8, "disp": True}
EPS = 1e-12

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def estimate(log_lam, mu_seg, mu_non_seg):
    """Negative log-likelihood for the Poisson mutation-rate model."""
    lam_val = np.exp(log_lam)

    p_seg = np.exp(-lam_val * mu_seg)
    p_seg = np.clip(p_seg, EPS, 1 - EPS)

    p_non = np.exp(-lam_val * mu_non_seg)
    p_non = np.clip(p_non, EPS, 1 - EPS)

    log_lik = np.sum(np.log(p_non)) + np.sum(np.log(1 - p_seg))
    return -log_lik


def main():
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    data_df = pd.read_csv(
        os.path.join(INPUT_DIR, INPUT_FILE),
        compression="gzip",
    )

    # If rows contain site counts rather than one row per site, expand here.
    # data_df = data_df.loc[
    #     data_df.index.repeat(data_df[SITE_COL]), [MU_COL, AC_COL]
    # ].reset_index(drop=True)

    mu_seg = np.array(data_df[data_df[AC_COL] != 0][MU_COL])
    mu_non_seg = np.array(data_df[data_df[AC_COL] == 0][MU_COL])

    res = minimize(
        estimate,
        np.log([INITIAL_LAMBDA]),
        method=OPTIMIZER_METHOD,
        args=(mu_seg, mu_non_seg),
        options=OPTIMIZER_OPTIONS,
    )

    lambda_est = float(np.exp(res.x[0]))
    print("Estimated lambda:", lambda_est)

    lam_val = pd.DataFrame([{"Lambda": lambda_est}])
    lam_val.to_csv(os.path.join(OUTPUT_DIR, OUTPUT_LAMBDA_FILE), index=False)

    lambd = float(math.exp(res.x[0]))
    df1 = pd.DataFrame()
    stop_filling = False

    for unq_mut in range(UNIQUE_MUTATION_LIMIT + 1):
        print(unq_mut)
        if not stop_filling:
            col_name = f"Unq muts sites_{unq_mut}"
            data_df[col_name] = poisson.pmf(unq_mut, lambd * data_df[MU_COL])
            poi_sum = sum(data_df[col_name])
            print(poi_sum)

            if poi_sum == 0:
                stop_filling = True
        else:
            poi_sum = 0

        data = pd.DataFrame([{"Unique mutation": unq_mut, "Unq muts sites": poi_sum}])
        df1 = pd.concat([df1, data], ignore_index=True)

    df1.to_csv(os.path.join(OUTPUT_DIR, OUTPUT_UNIQUE_MUTATION_FILE), index=False)


if __name__ == "__main__":
    main()
