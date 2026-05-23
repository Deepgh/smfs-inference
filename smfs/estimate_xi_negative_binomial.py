#!/usr/bin/env python3

import os
import sys
from pathlib import Path
from warnings import simplefilter

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import nbinom

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow_config import apply_config

CONFIG_SECTION = "smfs.estimate_xi_negative_binomial"

# Paths and files
INPUT_DIR = "/project/yuvalsim/Deep/project2/josh_full_data/error_data/"
INPUT_FILE = "sfs_genes_dist.csv.gz"
OUTPUT_DIR = INPUT_DIR
OUTPUT_XI_FILE = "lam_mu_sd_josh.csv"
SAVE_XI_OUTPUT = False

# Column names
SITES_COL = "Sites"
MEAN_COL = "Mean theta"
SD_COL = "SD theta"
AC_COL = "AC_nfe_down"

# Model settings
INITIAL_XI = 1e3
UNIQUE_MUTATION_LIMIT = 300
OUTPUT_UNIQUE_MUTATION_FILE = f"unq_mut_sites_{UNIQUE_MUTATION_LIMIT}_josh_mew_sd.csv"
OPTIMIZER_METHOD = "nelder-mead"
OPTIMIZER_OPTIONS = {"xatol": 1e-8, "disp": True}
EPS = 1e-12

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def p0(xi, mu, sig):
    """Probability of zero mutational events under the negative-binomial marginal."""
    mu = np.asarray(mu)
    sig = np.asarray(sig)

    mu_sig_sq = (mu / sig) ** 2
    return np.exp(-mu_sig_sq * np.log(1 + (xi * mu) / mu_sig_sq))


def neg_loglik(log_xi, mu_seg, mu_non, sig_seg, sig_non):
    xi = np.exp(log_xi)

    p0_non = p0(xi, mu_non, sig_non)
    p0_seg = p0(xi, mu_seg, sig_seg)

    p0_non = np.clip(p0_non, EPS, 1 - EPS)
    p0_seg = np.clip(p0_seg, EPS, 1 - EPS)

    ll = (
        np.sum(np.log(p0_non)) +
        np.sum(np.log(1 - p0_seg))
    )

    return -ll


def main(data_df):
    data_df = data_df.copy()

    # If rows contain site counts rather than one row per site, expand here.
    # data_df = (
    #     data_df
    #     .loc[data_df.index.repeat(data_df[SITES_COL]), [AC_COL, MEAN_COL, SD_COL]]
    #     .reset_index(drop=True)
    # )

    mu_seg = data_df.loc[data_df[AC_COL] != 0, MEAN_COL].values
    sig_seg = data_df.loc[data_df[AC_COL] != 0, SD_COL].values
    mu_non = data_df.loc[data_df[AC_COL] == 0, MEAN_COL].values
    sig_non = data_df.loc[data_df[AC_COL] == 0, SD_COL].values

    res = minimize(
        neg_loglik,
        x0=np.log([INITIAL_XI]),
        args=(mu_seg, mu_non, sig_seg, sig_non),
        method=OPTIMIZER_METHOD,
        options={**OPTIMIZER_OPTIONS, "disp": False},
    )

    xi_hat = float(np.exp(res.x[0]))

    xi_df = pd.DataFrame([{"Lambda": xi_hat}])

    unique_mutation_rows = []
    stop_filling = False

    for unq_mut in range(UNIQUE_MUTATION_LIMIT + 1):
        if not stop_filling:
            col_name = f"Unq muts sites_{unq_mut}"

            data_df[col_name] = nbinom.pmf(
                unq_mut,
                n=(data_df[MEAN_COL] / data_df[SD_COL]) ** 2,
                p=1 / (1 + data_df[SD_COL] ** 2 * (xi_hat / data_df[MEAN_COL])),
            )

            poi_sum = data_df[col_name].sum()

            if poi_sum == 0:
                stop_filling = True
        else:
            poi_sum = 0

        unique_mutation_rows.append(
            {"Unique mutation": unq_mut, "Unq muts sites": poi_sum}
        )

    unique_mutation_df = pd.DataFrame(unique_mutation_rows)
    return xi_hat, xi_df, unique_mutation_df


if __name__ == "__main__":
    apply_config(CONFIG_SECTION, globals())

    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    input_df = pd.read_csv(
        os.path.join(INPUT_DIR, INPUT_FILE),
        compression="gzip",
    )

    estimated_xi, xi_output_df, unique_mutation_output_df = main(input_df)

    print(f"Estimated xi: {estimated_xi:}")

    if SAVE_XI_OUTPUT:
        xi_output_df.to_csv(os.path.join(OUTPUT_DIR, OUTPUT_XI_FILE), index=False)

    unique_mutation_output_df.to_csv(
        os.path.join(OUTPUT_DIR, OUTPUT_UNIQUE_MUTATION_FILE),
        index=False,
    )
