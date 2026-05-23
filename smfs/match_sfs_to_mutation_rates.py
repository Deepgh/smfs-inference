#!/usr/bin/env python3

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow_config import apply_config

CONFIG_SECTION = "smfs.match_sfs_to_mutation_rates"

# Paths and files
SFS_FILE = "data/site_frequency_spectrum.csv.gz"
MUTATION_RATE_FILE = "data/mutation_rates.tsv"
OUTPUT_FILE = "results/matched_site_frequency_spectrum.csv.gz"
SAVE_OUTPUT = False

# Column names
REF_CONTEXT_COL = "ref_context"
ALT_CONTEXT_COL = "alt_context"
METHYLATION_LEVEL_COL = "methylation_level"
SITES_COL = "sites"
MUTATION_SUM_COL = "mut_sum"

MU_GNOMAD_COL = "mu_gnomad"
MEAN_THETA_COL = "mean_theta"
SD_THETA_COL = "sd_theta"
ERROR_COL = "error"
METHYLATION_LEVEL_OUT_COL = "methylation_level"
MU_OUT_COL = "mutation_rate"
MEAN_THETA_OUT_COL = "mean_theta"
SD_THETA_OUT_COL = "sd_theta"
ERROR_OUT_COL = "error"


def complement(seq):
    return "".join({"C": "G", "G": "C", "T": "A", "A": "T"}[base] for base in seq)


def swap_first_third(seq):
    return seq[2] + seq[1] + seq[0] if len(seq) == 3 else seq


def main(sfs_df, mutation_rate_df):
    mutation_rate_columns = [
        REF_CONTEXT_COL,
        ALT_CONTEXT_COL,
        METHYLATION_LEVEL_COL,
        MU_GNOMAD_COL,
        MEAN_THETA_COL,
        SD_THETA_COL,
        ERROR_COL,
    ]
    triplet_key = [REF_CONTEXT_COL, ALT_CONTEXT_COL, METHYLATION_LEVEL_COL]
    rename_columns = {
        METHYLATION_LEVEL_COL: METHYLATION_LEVEL_OUT_COL,
        MU_GNOMAD_COL: MU_OUT_COL,
        MEAN_THETA_COL: MEAN_THETA_OUT_COL,
        SD_THETA_COL: SD_THETA_OUT_COL,
        ERROR_COL: ERROR_OUT_COL,
    }

    sfs_df = sfs_df.copy()
    sfs_df[METHYLATION_LEVEL_COL] = sfs_df[METHYLATION_LEVEL_COL].fillna(0).astype(int)

    sfs_grouped = (
        sfs_df
        .groupby(triplet_key)[SITES_COL]
        .sum()
        .reset_index()
        .rename(columns={SITES_COL: MUTATION_SUM_COL})
    )

    mutation_rate_df = mutation_rate_df[mutation_rate_columns].copy()

    matched = pd.merge(mutation_rate_df, sfs_grouped, on=triplet_key, how="inner")

    sfs_set = set(tuple(x) for x in sfs_grouped[triplet_key].values)
    mutation_rate_set = set(tuple(x) for x in mutation_rate_df[triplet_key].values)
    unmatched_in_mutation_rates = mutation_rate_set - sfs_set

    mutation_rate_unmatched = mutation_rate_df[
        mutation_rate_df.apply(
            lambda row: tuple(row[triplet_key]) in unmatched_in_mutation_rates,
            axis=1,
        )
    ].copy()

    mutation_rate_unmatched[REF_CONTEXT_COL] = mutation_rate_unmatched[REF_CONTEXT_COL].apply(complement)
    mutation_rate_unmatched[ALT_CONTEXT_COL] = mutation_rate_unmatched[ALT_CONTEXT_COL].apply(complement)
    mutation_rate_unmatched[REF_CONTEXT_COL] = mutation_rate_unmatched[REF_CONTEXT_COL].apply(swap_first_third)
    mutation_rate_unmatched[ALT_CONTEXT_COL] = mutation_rate_unmatched[ALT_CONTEXT_COL].apply(swap_first_third)

    swap_match = pd.merge(mutation_rate_unmatched, sfs_grouped, on=triplet_key, how="inner")

    all_matches = pd.concat([matched, swap_match], ignore_index=True)

    matched_sfs = sfs_df.merge(all_matches, on=triplet_key, how="inner")
    matched_sfs = matched_sfs.rename(columns=rename_columns)
    matched_sfs = matched_sfs.drop(columns=[MUTATION_SUM_COL])

    match_counts = {
        "direct_matches": len(matched),
        "reverse_complement_matches": len(swap_match),
        "total_matches": len(all_matches),
    }
    return matched_sfs, match_counts


if __name__ == "__main__":
    apply_config(CONFIG_SECTION, globals())

    sfs_data = pd.read_csv(SFS_FILE, compression="gzip")
    mutation_rate_data = pd.read_csv(MUTATION_RATE_FILE, sep="\t")

    output_df, counts = main(sfs_data, mutation_rate_data)

    print(f"Common list 1 match: {counts['direct_matches']}")
    print(f"Swap match: {counts['reverse_complement_matches']}")
    print(f"Len of all: {counts['total_matches']}")

    if SAVE_OUTPUT:
        Path(OUTPUT_FILE).parent.mkdir(parents=True, exist_ok=True)
        output_df.to_csv(OUTPUT_FILE, compression="gzip", index=False)
