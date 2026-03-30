#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd


def sample_observed_mutations(
    simulation_df: pd.DataFrame,
    population_size: int,
    sample_size: int,
    allele_copies_col: str = "Allele copies",
    output_col: str = "Sampled mutation",
    random_seed: Optional[int] = None,
) -> pd.DataFrame:
    rng = np.random.default_rng(random_seed)
    result = simulation_df.copy()
    draws = rng.hypergeometric(
        ngood=result[allele_copies_col].to_numpy(dtype=np.int64),
        nbad=population_size - result[allele_copies_col].to_numpy(dtype=np.int64),
        nsample=sample_size,
    )
    result[output_col] = draws
    return result[result[output_col] != 0].copy()
