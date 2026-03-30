#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd


def complement(seq):
    return ''.join({'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}[base] for base in seq)


def swap_first_third(seq):
    return seq[2] + seq[1] + seq[0] if len(seq) == 3 else seq


def match_mutation_metadata(
    sfs_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    triplet_key=("ref_context", "alt_context", "methylation_level"),
) -> pd.DataFrame:
    sfs_grouped = (
        sfs_df
        .assign(methylation_level=sfs_df["methylation_level"].fillna(0).astype(int))
        .groupby(list(triplet_key))["Sites"]
        .sum()
        .reset_index(name="mut_sum")
    )

    metadata_cols = list(triplet_key) + ["mu_gnomad", "mean_theta", "sd_theta", "error"]
    metadata_trimmed = metadata_df[metadata_cols].copy()

    matched = pd.merge(metadata_trimmed, sfs_grouped, on=list(triplet_key), how="inner")

    sfs_set = set(tuple(x) for x in sfs_grouped[list(triplet_key)].values)
    metadata_set = set(tuple(x) for x in metadata_trimmed[list(triplet_key)].values)
    unmatched_in_metadata = metadata_set - sfs_set

    metadata_unmatched = metadata_trimmed[
        metadata_trimmed.apply(lambda row: tuple(row[list(triplet_key)]) in unmatched_in_metadata, axis=1)
    ].copy()
    metadata_unmatched["ref_context"] = metadata_unmatched["ref_context"].apply(complement).apply(swap_first_third)
    metadata_unmatched["alt_context"] = metadata_unmatched["alt_context"].apply(complement).apply(swap_first_third)

    swap_match = pd.merge(metadata_unmatched, sfs_grouped, on=list(triplet_key), how="inner")
    all_matches = pd.concat([matched, swap_match], ignore_index=True)

    result = sfs_df.merge(all_matches, on=list(triplet_key), how="inner")
    return result.rename(columns={
        "methylation_level": "Meth level",
        "mu_gnomad": "Mew",
        "mean_theta": "Mean theta",
        "sd_theta": "SD theta",
        "error": "Error",
    }).drop(columns=["mut_sum"])

