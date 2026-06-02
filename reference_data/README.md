# Reference Data

`mut_rates.csv` is the normalized mutation-rate reference table
distributed with this repository. `configs/template.toml` points to this file
for mutation-rate matching, optional error correction, and per-site mutation-rate
sampling.

The table uses the normalized `mu` column name expected by the default scripts.
It is derived from the AJHG supplementary mutation-rate table with unused
source metadata columns removed and `mu_gnomad` renamed to `mu`.
