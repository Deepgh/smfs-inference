#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import pandas as pd


def summarize_site_mutations(
    sampled_df: pd.DataFrame,
    site_col: str = "Site",
    mutation_id_col: str = "Mutation id",
    sampled_mutation_col: str = "Sampled mutation",
) -> pd.DataFrame:
    grouped = sampled_df.groupby(site_col, sort=True)
    result = grouped.agg({
        "Mutation number": "first",
        "Mutation rate": "first",
        "Mean theta": "first",
        "SD theta": "first",
        "Sampled mu": "first",
        "Meth level": "first",
        sampled_mutation_col: "sum",
    }).reset_index()
    result["Unique mutations"] = grouped[mutation_id_col].nunique().to_numpy()
    return result.rename(columns={sampled_mutation_col: "Sampled mutation sum"})


def fill_missing_site_rows(
    site_summary_df: pd.DataFrame,
    site_metadata_df: pd.DataFrame,
    site_col: str = "Site",
) -> pd.DataFrame:
    merged = site_metadata_df.merge(
        site_summary_df[[site_col, "Sampled mutation sum", "Unique mutations"]],
        on=site_col,
        how="left",
    )
    merged[["Sampled mutation sum", "Unique mutations"]] = (
        merged[["Sampled mutation sum", "Unique mutations"]]
        .fillna(0)
        .astype(int)
    )
    return merged.rename(columns={"Sampled mutation sum": "AC"})
