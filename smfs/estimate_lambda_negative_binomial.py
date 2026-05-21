#!/usr/bin/env python3

import os
from pathlib import Path
from warnings import simplefilter

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import nbinom

# Paths and files
INPUT_DIR = "/project/yuvalsim/Deep/project2/josh_full_data/error_data/"
INPUT_FILE = "sfs_genes_dist.csv.gz"
OUTPUT_DIR = INPUT_DIR
OUTPUT_LAMBDA_FILE = "lam_mu_sd_josh.csv"
SAVE_LAMBDA_OUTPUT = False

# Column names
SITES_COL = "Sites"
MEAN_COL = "Mean theta"
SD_COL = "SD theta"
AC_COL = "AC_nfe_down"

# Model settings
INITIAL_LAMBDA = 1e3
UNIQUE_MUTATION_LIMIT = 300
OUTPUT_UNIQUE_MUTATION_FILE = f"unq_mut_sites_{UNIQUE_MUTATION_LIMIT}_josh_mew_sd.csv"
OPTIMIZER_METHOD = "nelder-mead"
OPTIMIZER_OPTIONS = {"xatol": 1e-8, "disp": True}
EPS = 1e-12

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def p0(lam, mu, sig):
    """Probability of zero mutational events under the negative-binomial marginal."""
    mu = np.asarray(mu)
    sig = np.asarray(sig)

    mu_sig_sq = (mu / sig) ** 2
    return np.exp(-mu_sig_sq * np.log(1 + (lam * mu) / mu_sig_sq))


def neg_loglik(log_lam, mu_seg, mu_non, sig_seg, sig_non):
    lam = np.exp(log_lam)

    p0_non = p0(lam, mu_non, sig_non)
    p0_seg = p0(lam, mu_seg, sig_seg)

    p0_non = np.clip(p0_non, EPS, 1 - EPS)
    p0_seg = np.clip(p0_seg, EPS, 1 - EPS)

    ll = (
        np.sum(np.log(p0_non)) +
        np.sum(np.log(1 - p0_seg))
    )

    return -ll


def main():
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    data_df = pd.read_csv(
        os.path.join(INPUT_DIR, INPUT_FILE),
        compression="gzip",
    )

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
        x0=np.log([INITIAL_LAMBDA]),
        args=(mu_seg, mu_non, sig_seg, sig_non),
        method=OPTIMIZER_METHOD,
        options=OPTIMIZER_OPTIONS,
    )

    lam_hat = float(np.exp(res.x[0]))
    print(f"Estimated lambda: {lam_hat:}")

    lam_df = pd.DataFrame([{"Lambda": lam_hat}])
    if SAVE_LAMBDA_OUTPUT:
        lam_df.to_csv(os.path.join(OUTPUT_DIR, OUTPUT_LAMBDA_FILE), index=False)

    df1 = pd.DataFrame()
    stop_filling = False

    for unq_mut in range(UNIQUE_MUTATION_LIMIT + 1):
        print(unq_mut)

        if not stop_filling:
            col_name = f"Unq muts sites_{unq_mut}"

            data_df[col_name] = nbinom.pmf(
                unq_mut,
                n=(data_df[MEAN_COL] / data_df[SD_COL]) ** 2,
                p=1 / (1 + data_df[SD_COL] ** 2 * (lam_hat / data_df[MEAN_COL])),
            )

            poi_sum = data_df[col_name].sum()
            print(poi_sum)

            if poi_sum == 0:
                stop_filling = True
        else:
            poi_sum = 0

        df1 = pd.concat(
            [df1, pd.DataFrame([{"Unique mutation": unq_mut, "Unq muts sites": poi_sum}])],
            ignore_index=True,
        )

    df1.to_csv(os.path.join(OUTPUT_DIR, OUTPUT_UNIQUE_MUTATION_FILE), index=False)


if __name__ == "__main__":
    main()
