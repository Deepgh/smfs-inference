#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd


def gamma_shape_scale(mean: pd.Series, sd: pd.Series) -> tuple[pd.Series, pd.Series]:
    shape = (mean / sd) ** 2
    scale = (sd ** 2) / mean
    return shape, scale


def summarize_gamma_draws(
    metadata_df: pd.DataFrame,
    mean_col: str = "mean_theta",
    sd_col: str = "sd_theta",
    num_samples: int = 100_000,
    seed_col: str = "seed",
    random_seed: Optional[int] = None,
) -> pd.DataFrame:
    result = metadata_df.copy()
    shape, scale = gamma_shape_scale(result[mean_col], result[sd_col])
    result["k"] = shape
    result["theta"] = scale

    rng = np.random.default_rng(random_seed)
    if seed_col not in result.columns:
        result[seed_col] = rng.integers(0, int(1e9), size=len(result), dtype=np.int64)

    sampled_means = []
    sampled_sds = []
    for row in result.itertuples(index=False):
        row_rng = np.random.default_rng(int(getattr(row, seed_col)))
        samples = row_rng.gamma(shape=getattr(row, "k"), scale=getattr(row, "theta"), size=num_samples)
        sampled_means.append(float(np.mean(samples)))
        sampled_sds.append(float(np.std(samples)))

    result["sampled mu"] = sampled_means
    result["sampled sd"] = sampled_sds
    return result


def sample_site_rates(
    metadata_df: pd.DataFrame,
    num_sites_col: str = "num_sites",
    mean_col: str = "mean_theta",
    sd_col: str = "sd_theta",
    mutation_number_start: int = 1,
    mutation_rate_col: str = "mu_gnomad",
    methylation_col: str = "methylation_level",
    random_seed: Optional[int] = None,
) -> pd.DataFrame:
    metadata = metadata_df.copy()
    shape, scale = gamma_shape_scale(metadata[mean_col], metadata[sd_col])
    metadata["k"] = shape
    metadata["theta"] = scale

    rng = np.random.default_rng(random_seed)
    metadata["seed"] = rng.integers(0, int(1e9), size=len(metadata), dtype=np.int64)

    expanded_rows = []
    for mutation_number, row in enumerate(metadata.itertuples(index=False), start=mutation_number_start):
        row_rng = np.random.default_rng(int(row.seed))
        num_sites = int(getattr(row, num_sites_col))
        samples = row_rng.gamma(shape=row.k, scale=row.theta, size=num_sites)
        expanded_rows.append(pd.DataFrame({
            "Mutation number": np.full(num_sites, mutation_number, dtype=int),
            "Mutation rate": np.full(num_sites, getattr(row, mutation_rate_col), dtype=float),
            "Mean theta": np.full(num_sites, getattr(row, mean_col), dtype=float),
            "SD theta": np.full(num_sites, getattr(row, sd_col), dtype=float),
            "Meth level": np.full(num_sites, getattr(row, methylation_col), dtype=int),
            "Sampled mu": samples,
        }))

    if not expanded_rows:
        return pd.DataFrame(columns=[
            "Mutation number",
            "Mutation rate",
            "Mean theta",
            "SD theta",
            "Meth level",
            "Sampled mu",
        ])
    return pd.concat(expanded_rows, ignore_index=True)
