#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

from smfs.simulation.population_growth import simulate_population_growth_dataframe


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path")
    parser.add_argument("site_results_path")
    parser.add_argument("site_metadata_path")
    parser.add_argument("--initial-population", type=int, default=10**4)
    parser.add_argument("--final-population", type=int, default=10**8)
    parser.add_argument("--processes", type=int, default=None)
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    sampled_mu_df = pd.read_csv(args.input_path, compression="infer")
    site_results, site_metadata = simulate_population_growth_dataframe(
        sampled_mu_df,
        initial_population=args.initial_population,
        final_population=args.final_population,
        processes=args.processes,
        random_seed=args.seed,
    )
    site_results.to_csv(args.site_results_path, index=False)
    site_metadata.to_csv(args.site_metadata_path, index=False)


if __name__ == "__main__":
    main()
