#!/usr/bin/env python3

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow_config import apply_config

CONFIG_SECTION = "population_model.sample_site_mutation_rates"

# Paths and files
INPUT_FILE = (
    "/project/yuvalsim/Deep/project2/josh_full_data/error_data/"
    "ajhg_00004094_supp_table2_mut.tsv"
)
EXPANDED_OUTPUT_FILE = (
    "/project/yuvalsim/Deep/project2/josh_full_data/error_data/"
    "scaled_mu_diff_site.csv"
)
DIAGNOSTIC_OUTPUT_FILE = (
    "/project/yuvalsim/Deep/project2/human_data/sim_hum_data/"
    "final_approach/scaled_mu_from_musd.csv"
)

# Column names
MEAN_THETA_COL = "mean_theta"
SD_THETA_COL = "sd_theta"
SCALED_MU_COL = "scaled mu"
SCALED_SD_COL = "scaled sd"
GAMMA_SHAPE_COL = "k"
GAMMA_SCALE_COL = "theta"
SEED_COL = "seed"
NUM_SITES_COL = "num_sites"
MU_COL = "mu"
METHYLATION_LEVEL_COL = "methylation_level"
SAMPLED_MU_COL = "sampled mu"
SAMPLED_SD_COL = "sampled sd"
DIFF_COL = "diff"
DIFF_SD_COL = "diffsd"

# Output column names
MUTATION_NUMBER_OUT_COL = "Mutation number"
MUTATION_RATE_OUT_COL = "Mutation rate"
MEAN_THETA_OUT_COL = "Mean theta"
SD_THETA_OUT_COL = "SD theta"
METHYLATION_LEVEL_OUT_COL = "Meth_level"
SAMPLED_MU_OUT_COL = "Sampled mu"

# Sampling and diagnostic settings
RANDOM_SEED_MAX = int(1e9)
DIAGNOSTIC_SAMPLE_SIZE = 100000
SHOW_DIAGNOSTIC_PLOTS = False
WRITE_DIAGNOSTIC_OUTPUT = True

DROP_COLUMNS = [
    "LL_no_heterogeneity",
    "LL_heterogeneity",
    "CV",
    "LR",
    "ref_context",
    "alt_context",
    "p_value",
    "p_adj",
]


def sample_mean(k, theta, seed, n=DIAGNOSTIC_SAMPLE_SIZE):
    seed = int(seed)
    rng = np.random.default_rng(seed)
    samples = rng.gamma(shape=k, scale=theta, size=n)
    return pd.Series({SAMPLED_MU_COL: np.mean(samples), SAMPLED_SD_COL: np.std(samples)})


def add_gamma_parameters(df):
    df[SCALED_MU_COL] = df[MEAN_THETA_COL]
    df[SCALED_SD_COL] = df[SD_THETA_COL]
    df[GAMMA_SHAPE_COL] = (df[SCALED_MU_COL] / df[SCALED_SD_COL]) ** 2
    df[GAMMA_SCALE_COL] = (df[SCALED_SD_COL] ** 2) / df[SCALED_MU_COL]
    df[SEED_COL] = np.random.randint(0, RANDOM_SEED_MAX, size=len(df)).astype(int)
    return df


def write_expanded_site_rates(df):
    expanded_rows = []
    mut_num = 1

    for _, row in df.iterrows():
        print(mut_num)
        rng = np.random.default_rng(int(row[SEED_COL]))
        n = int(row[NUM_SITES_COL])
        samples = rng.gamma(shape=row[GAMMA_SHAPE_COL], scale=row[GAMMA_SCALE_COL], size=n)

        long_df = pd.DataFrame(
            {
                MUTATION_NUMBER_OUT_COL: [mut_num] * n,
                MUTATION_RATE_OUT_COL: [row[MU_COL]] * n,
                MEAN_THETA_OUT_COL: [row[MEAN_THETA_COL]] * n,
                SD_THETA_OUT_COL: [row[SD_THETA_COL]] * n,
                METHYLATION_LEVEL_OUT_COL: [row[METHYLATION_LEVEL_COL]] * n,
                SAMPLED_MU_OUT_COL: samples,
            }
        )
        mut_num += 1
        expanded_rows.append(long_df)

    result_df = pd.concat(expanded_rows, ignore_index=True)
    result_df.to_csv(EXPANDED_OUTPUT_FILE, index=False)


def plot_diagnostics(df):
    df[DIFF_COL] = (df[SAMPLED_MU_COL] - df[SCALED_MU_COL]) / df[SCALED_MU_COL]
    plt.plot(range(len(df[DIFF_COL])), df[DIFF_COL])
    plt.ylabel("relative error sampled and scaled mu")
    plt.show()

    plt.plot(df[SCALED_MU_COL], df[SAMPLED_MU_COL])
    plt.xlabel("scaled mu")
    plt.ylabel("mean sampled mu")
    plt.show()

    df[DIFF_SD_COL] = (df[SCALED_SD_COL] - df[SCALED_SD_COL]) / df[SCALED_SD_COL]
    plt.plot(range(len(df[DIFF_SD_COL])), df[DIFF_SD_COL])
    plt.ylabel("relative error sampled and scaledsd")
    plt.show()

    plt.plot(df[SCALED_SD_COL], df[SAMPLED_SD_COL])
    plt.xlabel("scaled sd")
    plt.ylabel("sampleld sd")
    plt.show()


def main():
    apply_config(CONFIG_SECTION, globals())

    df = pd.read_csv(INPUT_FILE, sep="\t")
    df = df.drop(DROP_COLUMNS, axis=1)
    df = add_gamma_parameters(df)

    write_expanded_site_rates(df)

    df[[SAMPLED_MU_COL, SAMPLED_SD_COL]] = df.apply(
        lambda row: sample_mean(
            row[GAMMA_SHAPE_COL],
            row[GAMMA_SCALE_COL],
            row[SEED_COL],
            DIAGNOSTIC_SAMPLE_SIZE,
        ),
        axis=1,
    )

    if SHOW_DIAGNOSTIC_PLOTS:
        plot_diagnostics(df)

    cleanup_cols = [GAMMA_SHAPE_COL, GAMMA_SCALE_COL, SEED_COL, DIFF_COL, DIFF_SD_COL]
    df = df.drop(columns=[col for col in cleanup_cols if col in df.columns])

    if WRITE_DIAGNOSTIC_OUTPUT:
        df.to_csv(DIAGNOSTIC_OUTPUT_FILE, index=False)


if __name__ == "__main__":
    main()
