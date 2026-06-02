# Reference Data

`ajhg_00004094_supp_table2_mut.tsv` is the mutation-rate reference table
distributed with this repository. `configs/template.toml` points to this file
for mutation-rate matching, optional error correction, and per-site mutation-rate
sampling.

The source table retains its original `mu_gnomad` column name. The template
maps that source column through `mu_col = "mu_gnomad"`.
