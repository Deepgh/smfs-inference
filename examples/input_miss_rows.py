#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

from smfs.data_prep.site_postprocessing import fill_missing_site_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("site_summary_path")
    parser.add_argument("site_metadata_path")
    parser.add_argument("output_path")
    parser.add_argument("--site-offset", type=int, default=0)
    args = parser.parse_args()

    site_summary_df = pd.read_csv(args.site_summary_path, compression="infer")
    site_metadata_df = pd.read_csv(args.site_metadata_path, compression="infer")
    result = fill_missing_site_rows(site_summary_df, site_metadata_df)
    if args.site_offset:
        result["Site"] = result["Site"] + args.site_offset
    result.to_csv(args.output_path, compression="infer", index=False)


if __name__ == "__main__":
    main()
