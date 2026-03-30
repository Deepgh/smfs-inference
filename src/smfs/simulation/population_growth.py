#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import math
import multiprocessing
from typing import Optional

import numpy as np
import pandas as pd


def growth_time_scale(initial_population: int, final_population: int, scale_factor: int = 250) -> int:
    return math.ceil(scale_factor * math.log10(final_population / initial_population))


def _offspring_sample(mean_count: float, rng: np.random.Generator) -> int:
    if mean_count < 1e9:
        return int(rng.poisson(mean_count))
    return int(rng.normal(mean_count, np.sqrt(mean_count)))


def simulate_mutations_over_time(
    time_scale: int,
    growth_factor: float,
    initial_population: int,
    sampled_mu: float,
    rng: np.random.Generator,
    mutation_id_start: int = 0,
) -> list[tuple[int, int]]:
    current_population = initial_population
    mutation_data: list[list[int]] = []
    results: list[tuple[int, int]] = []
    mutation_id = mutation_id_start

    for t in reversed(range(1, time_scale + 1)):
        theta = 2 * current_population * sampled_mu
        new_mutations = int(rng.poisson(theta))

        for _ in range(new_mutations):
            mutation_data.append([1, mutation_id])
            mutation_id += 1

        for mutation in mutation_data:
            mutation[0] = _offspring_sample(growth_factor * mutation[0], rng)

        mutation_data = [mutation for mutation in mutation_data if mutation[0] != 0]

        if t == 1:
            results.extend((allele_copies, mut_id) for allele_copies, mut_id in mutation_data)

        current_population *= growth_factor

    return results


def _simulate_single_site(args) -> tuple[pd.DataFrame, pd.DataFrame]:
    (
        site_index,
        row_dict,
        initial_population,
        final_population,
        base_seed,
    ) = args
    time_scale = growth_time_scale(initial_population, final_population)
    growth_factor = (final_population / initial_population) ** (1 / time_scale)
    rng = np.random.default_rng(None if base_seed is None else base_seed + site_index)

    metadata = {
        "Mutation number": row_dict["Mutation number"],
        "Mutation rate": row_dict["Mutation rate"],
        "Mean theta": row_dict["Mean theta"],
        "SD theta": row_dict["SD theta"],
        "Meth level": row_dict["Meth level"],
        "Sampled mu": row_dict["Sampled mu"],
    }

    site_metadata = pd.DataFrame([{"Site": site_index, **metadata}])
    sim_results = simulate_mutations_over_time(
        time_scale=time_scale,
        growth_factor=growth_factor,
        initial_population=initial_population,
        sampled_mu=row_dict["Sampled mu"],
        rng=rng,
        mutation_id_start=0,
    )
    site_results = pd.DataFrame([
        {
            "Site": site_index,
            **metadata,
            "Allele copies": allele_copies,
            "Mutation id": mutation_id,
        }
        for allele_copies, mutation_id in sim_results
    ])
    return site_results, site_metadata


def simulate_population_growth_dataframe(
    sampled_mu_df: pd.DataFrame,
    initial_population: int = 10**4,
    final_population: int = 10**8,
    processes: Optional[int] = None,
    random_seed: Optional[int] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    args = [
        (
            site_index,
            {
                "Mutation number": row["Mutation number"],
                "Mutation rate": row["Mutation rate"],
                "Mean theta": row["Mean theta"],
                "SD theta": row["SD theta"],
                "Meth level": row.get("Meth level", row.get("Meth_level")),
                "Sampled mu": row["Sampled mu"],
            },
            initial_population,
            final_population,
            random_seed,
        )
        for site_index, row in enumerate(sampled_mu_df.to_dict("records"))
    ]

    if processes == 1:
        results = [_simulate_single_site(arg) for arg in args]
    else:
        worker_count = processes or multiprocessing.cpu_count()
        with multiprocessing.Pool(worker_count) as pool:
            results = pool.map(_simulate_single_site, args)

    site_result_frames = [result[0] for result in results if not result[0].empty]
    if site_result_frames:
        site_results = pd.concat(site_result_frames, ignore_index=True).sort_values("Site")
    else:
        site_results = pd.DataFrame(columns=[
            "Site",
            "Mutation number",
            "Mutation rate",
            "Mean theta",
            "SD theta",
            "Meth level",
            "Sampled mu",
            "Allele copies",
            "Mutation id",
        ])
    site_metadata = pd.concat([result[1] for result in results], ignore_index=True).sort_values("Site")
    return site_results, site_metadata
