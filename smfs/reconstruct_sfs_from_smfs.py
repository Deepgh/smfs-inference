#!/usr/bin/env python3

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import nbinom

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow_config import apply_config

CONFIG_SECTION = "smfs.reconstruct_sfs_from_smfs"

# Paths and files
INPUT_DIR = "data"
INPUT_FILE = "observed_site_frequency_spectrum.csv.gz"
PREDICTED_PROBABILITY_DIR = "results"
PREDICTED_PROBABILITY_FILE = "smfs_observed_300.csv"
XI_DIR = "results"
XI_FILE = "xi_negative_binomial.csv"
OUTPUT_DIR = "results"
OUTPUT_FILE = "reconstructed_sfs_by_context.csv"

# Input column names
AC_COL = "allele_count"
SITES_COL = "sites"
MEAN_COL = "mean_theta"
SD_COL = "sd_theta"
METHYLATION_LEVEL_COL = "methylation_level"
PREDICTED_PROBABILITY_COL = "predicted_probability"
XI_COL = "xi"

# Internal column names
MU_SD_RATIO_COL = "mu_sd_ratio"
NUMERATOR_FACTOR_COL = "numerator_factor"

# Output column names
AC_OUT_COL = "allele_count"
METHYLATION_LEVEL_OUT_COL = "methylation_level"
MEAN_OUT_COL = "mean_theta"
SD_OUT_COL = "sd_theta"
NUM_SITES_OUT_COL = "num_sites"
PROP_SITES_OUT_COL = "proportion_sites"
TYPE_OUT_COL = "type"
DATA_TYPE_LABEL = "data"
THEORY_TYPE_LABEL = "theory"

# Model settings
AC_LIMIT = 200
ROWS_PER_CONTEXT = 10**6 + 1


def main(pred_df, muts_data_df, xi_df):
    prob_counts = list(pred_df[PREDICTED_PROBABILITY_COL])
    prob_counts.insert(0, 0)

    conv_lists = [prob_counts]
    conv = prob_counts

    for _ in range(AC_LIMIT + 1):
        conv = np.convolve(prob_counts, conv)
        conv_lists.append(conv)

    xi = list(xi_df[XI_COL])[0]

    muts_data_df[MU_SD_RATIO_COL] = (
        muts_data_df[MEAN_COL] / muts_data_df[SD_COL]
    ) ** 2
    muts_data_df[NUMERATOR_FACTOR_COL] = (
        xi * muts_data_df[MEAN_COL]
    ) / muts_data_df[MU_SD_RATIO_COL]

    mut = 1
    prop_sites_df = pd.DataFrame()
    all_sum_sites_data = pd.DataFrame()

    for i in range(0, len(muts_data_df), ROWS_PER_CONTEXT):
        df = muts_data_df[i:i + ROWS_PER_CONTEXT]

        gamma_mew = np.unique(df[MEAN_COL])[0]
        gamma_sigma = np.unique(df[SD_COL])[0]
        mu_sigma_rat = np.unique(df[MU_SD_RATIO_COL])[0]
        mu_sig_nume = np.unique(df[NUMERATOR_FACTOR_COL])[0]
        meth = np.unique(df[METHYLATION_LEVEL_COL])[0]

        mean = np.unique(df[MEAN_COL])[0]
        sd = np.unique(df[SD_COL])[0]

        df1 = pd.DataFrame()
        tot_sites = sum(df[SITES_COL])

        site_0 = (df[df[AC_COL] == 0])[SITES_COL].iloc[0]
        pred_0 = tot_sites * (1 / ((1 + mu_sig_nume) ** mu_sigma_rat))

        pred_site_prop_0 = pred_0 / tot_sites
        data_site_prop_0 = site_0 / tot_sites

        prop_sites_data = pd.DataFrame(
            [
                {
                    AC_OUT_COL: 0,
                    METHYLATION_LEVEL_OUT_COL: meth,
                    MEAN_OUT_COL: gamma_mew,
                    SD_OUT_COL: gamma_sigma,
                    PROP_SITES_OUT_COL: data_site_prop_0,
                    TYPE_OUT_COL: DATA_TYPE_LABEL,
                }
            ]
        )
        prop_sites_theory = pd.DataFrame(
            [
                {
                    AC_OUT_COL: 0,
                    METHYLATION_LEVEL_OUT_COL: meth,
                    MEAN_OUT_COL: gamma_mew,
                    SD_OUT_COL: gamma_sigma,
                    PROP_SITES_OUT_COL: pred_site_prop_0,
                    TYPE_OUT_COL: THEORY_TYPE_LABEL,
                }
            ]
        )
        prop_sites_df = pd.concat(
            [prop_sites_df, prop_sites_data, prop_sites_theory],
            ignore_index=True,
        )

        data_data = pd.DataFrame(
            [
                {
                    AC_OUT_COL: 0,
                    METHYLATION_LEVEL_OUT_COL: meth,
                    MEAN_OUT_COL: gamma_mew,
                    SD_OUT_COL: gamma_sigma,
                    NUM_SITES_OUT_COL: site_0,
                    TYPE_OUT_COL: DATA_TYPE_LABEL,
                }
            ]
        )
        theory_data = pd.DataFrame(
            [
                {
                    AC_OUT_COL: 0,
                    METHYLATION_LEVEL_OUT_COL: meth,
                    MEAN_OUT_COL: gamma_mew,
                    SD_OUT_COL: gamma_sigma,
                    NUM_SITES_OUT_COL: pred_0,
                    TYPE_OUT_COL: THEORY_TYPE_LABEL,
                }
            ]
        )
        df1 = pd.concat([df1, data_data, theory_data], ignore_index=True)

        for copy in range(1, AC_LIMIT + 1):
            num_site = (df[df[AC_COL] == copy])[SITES_COL].iloc[0]

            pred_sum = 0
            for l in range(copy):
                poi_coeff = nbinom.pmf(
                    l + 1,
                    n=(mean / sd) ** 2,
                    p=1 / (1 + sd**2 * (xi / mean)),
                )
                conv_prob = poi_coeff * conv_lists[l][copy]
                pred_sum += conv_prob

            pred_site = tot_sites * pred_sum

            data_data = pd.DataFrame(
                [
                    {
                        AC_OUT_COL: copy,
                        METHYLATION_LEVEL_OUT_COL: meth,
                        MEAN_OUT_COL: gamma_mew,
                        SD_OUT_COL: gamma_sigma,
                        NUM_SITES_OUT_COL: num_site,
                        TYPE_OUT_COL: DATA_TYPE_LABEL,
                    }
                ]
            )
            theory_data = pd.DataFrame(
                [
                    {
                        AC_OUT_COL: copy,
                        METHYLATION_LEVEL_OUT_COL: meth,
                        MEAN_OUT_COL: gamma_mew,
                        SD_OUT_COL: gamma_sigma,
                        NUM_SITES_OUT_COL: pred_site,
                        TYPE_OUT_COL: THEORY_TYPE_LABEL,
                    }
                ]
            )
            df1 = pd.concat([df1, data_data, theory_data], ignore_index=True)

            data_site_prop = num_site / tot_sites

            prop_sites_data = pd.DataFrame(
                [
                    {
                        AC_OUT_COL: copy,
                        METHYLATION_LEVEL_OUT_COL: meth,
                        MEAN_OUT_COL: gamma_mew,
                        SD_OUT_COL: gamma_sigma,
                        PROP_SITES_OUT_COL: data_site_prop,
                        TYPE_OUT_COL: DATA_TYPE_LABEL,
                    }
                ]
            )
            prop_sites_theory = pd.DataFrame(
                [
                    {
                        AC_OUT_COL: copy,
                        METHYLATION_LEVEL_OUT_COL: meth,
                        MEAN_OUT_COL: gamma_mew,
                        SD_OUT_COL: gamma_sigma,
                        PROP_SITES_OUT_COL: pred_sum,
                        TYPE_OUT_COL: THEORY_TYPE_LABEL,
                    }
                ]
            )
            prop_sites_df = pd.concat(
                [prop_sites_df, prop_sites_data, prop_sites_theory],
                ignore_index=True,
            )

        all_sum_sites_data = pd.concat([all_sum_sites_data, df1], ignore_index=True)
        print(mut)
        mut += 1

    return all_sum_sites_data, prop_sites_df


if __name__ == "__main__":
    apply_config(CONFIG_SECTION, globals())

    predicted_probability_df = pd.read_csv(
        os.path.join(PREDICTED_PROBABILITY_DIR, PREDICTED_PROBABILITY_FILE)
    )
    mutation_data_df = pd.read_csv(
        os.path.join(INPUT_DIR, INPUT_FILE),
        compression="infer",
    )
    xi_df = pd.read_csv(os.path.join(XI_DIR, XI_FILE))

    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    output_df, _ = main(predicted_probability_df, mutation_data_df, xi_df)
    output_df.to_csv(os.path.join(OUTPUT_DIR, OUTPUT_FILE), index=False)
