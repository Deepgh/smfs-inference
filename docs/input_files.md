# Input File Requirements

This repository contains standalone research scripts. Each script can use either the constants defined near the top of the script or a TOML config file passed with `--config`.

The tables below describe the *meaning* of each required input column, the config variable that controls the column name, and the current default column name. Edit the config variable when your data use a different column name.

For example:

```toml
[smfs.estimate_xi_poisson]
mu_col = "my_sampled_mu_column"
ac_col = "allele_count"
```

The original project data are not included in this repository.

## SMFS Scripts

### `smfs/match_sfs_to_mutation_rates.py`

Matches SFS rows to mutation-rate metadata by reference context, alternate context, and methylation level.

Input files:

| Config variable | Description |
| --- | --- |
| `sfs_file` | SFS table, read as gzipped CSV. |
| `mutation_rate_file` | Mutation-rate metadata table, read as tab-separated text. |

Required `sfs_file` columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Reference sequence context | `ref_context_col` | `ref_context` |
| Alternate sequence context | `alt_context_col` | `alt_context` |
| Methylation level | `methylation_level_col` | `methylation_level` |
| Number of sites for the SFS row | `sites_col` | `Sites` |

Required `mutation_rate_file` columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Reference sequence context | `ref_context_col` | `ref_context` |
| Alternate sequence context | `alt_context_col` | `alt_context` |
| Methylation level | `methylation_level_col` | `methylation_level` |
| gnomAD mutation rate | `mu_gnomad_col` | `mu_gnomad` |
| Mean theta | `mean_theta_col` | `mean_theta` |
| Standard deviation of theta | `sd_theta_col` | `sd_theta` |
| Error estimate | `error_col` | `error` |

Output:

| Config variable | Description |
| --- | --- |
| `output_file` | Matched SFS table with mutation-rate metadata merged in. Written only when `save_output = true`. |

Output column names are controlled by `methylation_level_out_col`, `mu_out_col`, `mean_theta_out_col`, `sd_theta_out_col`, and `error_out_col`.

### `smfs/estimate_xi_poisson.py`

Estimates xi using a Poisson model with site-specific sampled mutation rates.

Input file:

| Config variable | Description |
| --- | --- |
| `input_dir` / `input_file` | Simulation or site table, read as gzipped CSV. |

Required columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Site identifier | `site_col` | `Site` |
| Site-specific sampled mutation rate | `mu_col` | `Sampled mu` |
| Allele count | `ac_col` | `AC` |

Outputs:

| Config variable | Description |
| --- | --- |
| `output_xi_file` | CSV with the estimated xi value. The output column is currently named `Lambda` for compatibility. |
| `output_unique_mutation_file` | CSV with expected site counts by unique-mutation count. Columns are currently `Unique mutation` and `Unq muts sites`. |

### `smfs/estimate_xi_negative_binomial.py`

Estimates xi while accounting for mean and standard deviation uncertainty using a negative-binomial marginal model.

Input file:

| Config variable | Description |
| --- | --- |
| `input_dir` / `input_file` | Site or SFS table, read as gzipped CSV. |

Required columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Number of sites represented by the row | `sites_col` | `Sites` |
| Mean theta | `mean_col` | `Mean theta` |
| Standard deviation of theta | `sd_col` | `SD theta` |
| Allele count | `ac_col` | `AC_nfe_down` |

Outputs:

| Config variable | Description |
| --- | --- |
| `output_xi_file` | CSV with the estimated xi value. Written only when `save_xi_output = true`. The output column is currently named `Lambda` for compatibility. |
| `output_unique_mutation_file` | CSV with expected site counts by unique-mutation count. Columns are currently `Unique mutation` and `Unq muts sites`. |

### `smfs/estimate_smfs_from_unique_mutations.py`

Numerically estimates predicted SMFS probabilities using observed allele counts and an expected unique-mutation count distribution.

Input files:

| Config variable | Description |
| --- | --- |
| `input_dir` / `input_file` | Observed allele-count table, read as gzipped CSV. |
| `unique_mutation_sites_dir` / `unique_mutation_sites_file` | Expected unique-mutation count distribution. |

Required `input_file` columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Allele count | `ac_col` | `AC_` |
| Number of sites at that allele count | `sites_col` | `Sites` |

Required `unique_mutation_sites_file` columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Expected number of sites with a given number of unique mutations | `unique_mutation_sites_col` | `Unq muts sites` |

Optional error-correction input:

| Config variable | Description |
| --- | --- |
| `error_file` | Error table used only when `apply_error_correction = true`. |

Required `error_file` columns when error correction is enabled:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Number of sites | `error_num_sites_col` | `num_sites` |
| Error rate | `error_rate_col` | `error` |

Output:

| Config variable | Description |
| --- | --- |
| `output_file` | Predicted SMFS probabilities. Columns are currently `Count` and `pred prob`. |

## Population Model Scripts

### `population_model/sample_site_mutation_rates.py`

Samples per-site mutation rates from mutation-rate metadata with mean and standard deviation values.

Input file:

| Config variable | Description |
| --- | --- |
| `input_file` | Mutation-rate metadata table, read as tab-separated text. |

Required columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Mean theta | `mean_theta_col` | `mean_theta` |
| Standard deviation of theta | `sd_theta_col` | `sd_theta` |
| Number of sites represented by the row | `num_sites_col` | `num_sites` |
| gnomAD mutation rate | `mu_gnomad_col` | `mu_gnomad` |
| Methylation level | `methylation_level_col` | `methylation_level` |

The script also drops these columns if present in the input: `LL_no_heterogeneity`, `LL_heterogeneity`, `CV`, `LR`, `ref_context`, `alt_context`, `p_value`, `p_adj`.

Outputs:

| Config variable | Description |
| --- | --- |
| `expanded_output_file` | Expanded per-site sampled mutation-rate table. |
| `diagnostic_output_file` | Metadata with sampled diagnostic means and standard deviations. Written only when `write_diagnostic_output = true`. |

Expanded output column names are controlled by `mutation_number_out_col`, `mutation_rate_out_col`, `mean_theta_out_col`, `sd_theta_out_col`, `methylation_level_out_col`, and `sampled_mu_out_col`.

### `population_model/simulate_population_growth_mutations.py`

Simulates mutation accumulation through a growing population using sampled site-specific mutation rates.

Input file:

| Config variable | Description |
| --- | --- |
| `input_dir` / `input_file` | Expanded per-site mutation-rate table, usually generated by `sample_site_mutation_rates.py`. |

Required columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Mutation type/row identifier | `mutation_number_col` | `Mutation number` |
| Mutation rate metadata | `mutation_rate_col` | `Mutation rate` |
| Mean theta | `mean_theta_col` | `Mean theta` |
| Standard deviation of theta | `sd_theta_col` | `SD theta` |
| Methylation level | `input_methylation_level_col` | `Meth_level` |
| Site-specific sampled mutation rate | `sampled_mu_col` | `Sampled mu` |

Outputs:

| Config variable | Description |
| --- | --- |
| `output_file` | Simulated mutation observations. |
| `mutation_rate_output_file` | Per-site mutation-rate table used in the simulation. |

Generated output column names are controlled by `site_col`, `methylation_level_col`, `allele_copies_col`, and `mutation_id_col`, along with the metadata column variables listed above.

### `population_model/sample_observed_mutations.py`

Samples observed mutations from simulated allele-copy counts using a hypergeometric model.

Input file:

| Config variable | Description |
| --- | --- |
| `input_dir` / `input_file` | Simulated mutation observations, usually generated by `simulate_population_growth_mutations.py`. |

Required columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Number of allele copies in the simulated population | `allele_copies_col` | `Allele copies` |

Output:

| Config variable | Description |
| --- | --- |
| `output_dir` / `output_file` | Sampled mutation table filtered to rows where the sampled mutation count is nonzero. |

The sampled mutation output column is controlled by `sampled_mutation_col`.

### `population_model/summarize_derived_alleles_by_site.py`

Aggregates sampled mutation records by site.

Input file:

| Config variable | Description |
| --- | --- |
| `input_dir` / `input_file` | Sampled mutation table, usually generated by `sample_observed_mutations.py`. |

Required columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Site identifier | `site_col` | `Site` |
| Mutation identifier | `mutation_id_col` | `Mutation id` |
| Methylation level | `methylation_level_col` | `Meth level` |
| Mutation rate metadata | `mutation_rate_col` | `Mutation rate` |
| Mutation type/row identifier | `mutation_number_col` | `Mutation number` |
| Mean theta | `mean_theta_col` | `Mean theta` |
| Standard deviation of theta | `sd_theta_col` | `SD theta` |
| Site-specific sampled mutation rate | `sampled_mu_col` | `Sampled mu` |
| Sampled mutation count | `sampled_mutation_col` | `Sampled mutation` |

Output:

| Config variable | Description |
| --- | --- |
| `output_dir` / `output_file` | Per-site derived-allele summary. |

The output currently includes compatibility column names such as `Mew theta`, `Meth level`, `Sampled mutation sum`, and `Unique mutations`.

### `population_model/fill_missing_simulation_sites.py`

Merges per-site derived-allele summaries back into the full mutation-rate table so sites with no observed sampled mutations are retained.

Input files:

| Config variable | Description |
| --- | --- |
| `derived_allele_dir` / `derived_allele_file` | Per-site summary, usually generated by `summarize_derived_alleles_by_site.py`. |
| `mutation_rate_dir` / `mutation_rate_file` | Full mutation-rate table, usually generated by `simulate_population_growth_mutations.py`. |

Required `derived_allele_file` columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Site identifier | `site_col` | `Site` |
| Total sampled mutation count for the site | `sampled_mutation_sum_col` | `Sampled mutation sum` |
| Number of unique mutations observed at the site | `unique_mutations_col` | `Unique mutations` |

Required `mutation_rate_file` columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Site identifier | `site_col` | `Site` |

Output:

| Config variable | Description |
| --- | --- |
| `output_dir` / `output_file` | Merged simulation table with missing sampled values filled with zero. |

The sampled mutation sum column is renamed to the value of `ac_col`, which defaults to `AC`.

### `population_model/estimate_simulated_smfs_from_unique_mutations.py`

Numerically estimates predicted SMFS probabilities from simulated data and an expected unique-mutation count distribution.

Input files:

| Config variable | Description |
| --- | --- |
| `input_dir` / `input_file` | Edited simulation data, usually generated by `fill_missing_simulation_sites.py`. |
| `unique_mutation_sites_dir` / `unique_mutation_sites_file` | Expected unique-mutation count distribution. |

Required `input_file` columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Allele count | `ac_col` | `AC` |

Required `unique_mutation_sites_file` columns:

| Meaning | Config variable | Default column name |
| --- | --- | --- |
| Expected number of sites with a given number of unique mutations | `unique_mutation_sites_col` | `Unq muts sites` |

Output:

| Config variable | Description |
| --- | --- |
| `output_dir` / `output_file` | Predicted simulated SMFS probabilities. Columns are currently `Count` and `pred prob`. |
