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
INPUT_FILE = "data/mutation_rates.tsv"
EXPANDED_OUTPUT_FILE = "results/sampled_site_mutation_rates.csv"
DIAGNOSTIC_OUTPUT_FILE = "results/mutation_rate_sampling_diagnostics.csv"

# Column names
MEAN_THETA_COL = "mean_theta"
SD_THETA_COL = "sd_theta"
SCALED_MU_COL = "scaled_mu"
SCALED_SD_COL = "scaled_sd"
GAMMA_SHAPE_COL = "gamma_shape"
GAMMA_SCALE_COL = "gamma_scale"
SEED_COL = "seed"
NUM_SITES_COL = "num_sites"
MU_GNOMAD_COL = "mu_gnomad"
METHYLATION_LEVEL_COL = "methylation_level"
SAMPLED_MU_COL = "sampled_mu"
SAMPLED_SD_COL = "sampled_sd"
DIFF_COL = "relative_mu_error"
DIFF_SD_COL = "relative_sd_error"

# Output column names
MUTATION_NUMBER_OUT_COL = "mutation_number"
MUTATION_RATE_OUT_COL = "mutation_rate"
MEAN_THETA_OUT_COL = "mean_theta"
SD_THETA_OUT_COL = "sd_theta"
METHYLATION_LEVEL_OUT_COL = "methylation_level"
SAMPLED_MU_OUT_COL = "sampled_mu"

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


def build_expanded_site_rates(df):
    expanded_rows = []
    mut_num = 1

    for _, row in df.iterrows():
        rng = np.random.default_rng(int(row[SEED_COL]))
        n = int(row[NUM_SITES_COL])
        samples = rng.gamma(shape=row[GAMMA_SHAPE_COL], scale=row[GAMMA_SCALE_COL], size=n)

        long_df = pd.DataFrame(
            {
                MUTATION_NUMBER_OUT_COL: [mut_num] * n,
                MUTATION_RATE_OUT_COL: [row[MU_GNOMAD_COL]] * n,
                MEAN_THETA_OUT_COL: [row[MEAN_THETA_COL]] * n,
                SD_THETA_OUT_COL: [row[SD_THETA_COL]] * n,
                METHYLATION_LEVEL_OUT_COL: [row[METHYLATION_LEVEL_COL]] * n,
                SAMPLED_MU_OUT_COL: samples,
            }
        )
        mut_num += 1
        expanded_rows.append(long_df)

    return pd.concat(expanded_rows, ignore_index=True)


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


def main(df):
    df = df.copy()
    df = df.drop(DROP_COLUMNS, axis=1)
    df = add_gamma_parameters(df)

    expanded_df = build_expanded_site_rates(df)

    df[[SAMPLED_MU_COL, SAMPLED_SD_COL]] = df.apply(
        lambda row: sample_mean(
            row[GAMMA_SHAPE_COL],
            row[GAMMA_SCALE_COL],
            row[SEED_COL],
            DIAGNOSTIC_SAMPLE_SIZE,
        ),
        axis=1,
    )

    cleanup_cols = [GAMMA_SHAPE_COL, GAMMA_SCALE_COL, SEED_COL, DIFF_COL, DIFF_SD_COL]
    diagnostic_df = df.drop(columns=[col for col in cleanup_cols if col in df.columns])
    return expanded_df, diagnostic_df


if __name__ == "__main__":
    apply_config(CONFIG_SECTION, globals())

    input_df = pd.read_csv(INPUT_FILE, sep="\t")
    expanded_output_df, diagnostic_output_df = main(input_df)

    Path(EXPANDED_OUTPUT_FILE).parent.mkdir(parents=True, exist_ok=True)
    expanded_output_df.to_csv(EXPANDED_OUTPUT_FILE, index=False)

    if SHOW_DIAGNOSTIC_PLOTS:
        plot_diagnostics(diagnostic_output_df.copy())

    if WRITE_DIAGNOSTIC_OUTPUT:
        Path(DIAGNOSTIC_OUTPUT_FILE).parent.mkdir(parents=True, exist_ok=True)
        diagnostic_output_df.to_csv(DIAGNOSTIC_OUTPUT_FILE, index=False)
