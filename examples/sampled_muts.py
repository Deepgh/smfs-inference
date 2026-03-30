#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

from smfs.simulation.downsampling import sample_observed_mutations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path")
    parser.add_argument("output_path")
    parser.add_argument("--population-size", type=int, required=True)
    parser.add_argument("--sample-size", type=int, required=True)
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    df = pd.read_csv(args.input_path, compression="infer")
    sampled = sample_observed_mutations(
        df,
        population_size=args.population_size,
        sample_size=args.sample_size,
        random_seed=args.seed,
    )
    sampled.to_csv(args.output_path, compression="infer", index=False)


if __name__ == "__main__":
    main()
