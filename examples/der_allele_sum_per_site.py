#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

from smfs.data_prep.site_postprocessing import summarize_site_mutations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path")
    parser.add_argument("output_path")
    args = parser.parse_args()

    df = pd.read_csv(args.input_path, compression="infer")
    result = summarize_site_mutations(df)
    result.to_csv(args.output_path, compression="infer", index=False)


if __name__ == "__main__":
    main()
