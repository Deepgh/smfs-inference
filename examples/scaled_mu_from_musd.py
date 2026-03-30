#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

from smfs.data_prep.mutation_rate_sampling import sample_site_rates


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path")
    parser.add_argument("output_path")
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    df = pd.read_csv(args.input_path, sep="\t" if args.input_path.endswith(".tsv") else ",")
    result = sample_site_rates(df, random_seed=args.seed)
    result.to_csv(args.output_path, index=False)


if __name__ == "__main__":
    main()
