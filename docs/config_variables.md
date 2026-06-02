# Config Variables

Scripts can be run with an optional TOML config file:

```bash
python smfs/estimate_xi_poisson.py --config configs/template.toml
```

If `--config` is not provided, each script uses the constants defined at the top of the script. The file `configs/current_behavior.toml` preserves the original hard-coded behavior. The file `configs/template.toml` is intended as a user-editable starting point.

Config keys are lowercase versions of the uppercase constants in each script. For example, `INPUT_FILE` becomes `input_file`, and `SAVE_OUTPUT` becomes `save_output`.

## Shared Notes

- Paths may be absolute or relative to the directory where you run the command.
- Column-name settings are optional. Only change them if your input files use different column names.
- Existing output column names and generated filenames are preserved for compatibility unless explicitly changed in the config.
- Boolean values use TOML syntax: `true` or `false`.

## `smfs.match_sfs_to_mutation_rates`

| Variable | Description |
| --- | --- |
| `sfs_file` | Input SFS table. |
| `mutation_rate_file` | Input mutation-rate metadata table. |
| `output_file` | Output matched SFS file. |
| `save_output` | Whether to write `output_file`. |
| `ref_context_col` | Reference-context column. |
| `alt_context_col` | Alternate-context column. |
| `methylation_level_col` | Methylation-level column. |
| `sites_col` | Site-count column. |
| `mutation_sum_col` | Temporary grouped mutation-count column. |
| `mu_col` | Input gnomAD mutation-rate column. |
| `mean_theta_col` | Input mean-theta column. |
| `sd_theta_col` | Input standard-deviation theta column. |
| `error_col` | Input error column. |
| `methylation_level_out_col` | Output methylation-level column. |
| `mu_out_col` | Output mutation-rate column. |
| `mean_theta_out_col` | Output mean-theta column. |
| `sd_theta_out_col` | Output standard-deviation theta column. |
| `error_out_col` | Output error column. |

## `smfs.estimate_xi_poisson`

| Variable | Description |
| --- | --- |
| `input_dir` | Directory containing `input_file`. |
| `input_file` | Input simulation/site table. |
| `output_dir` | Output directory. |
| `output_xi_file` | Output file for the estimated xi value. |
| `site_col` | Site identifier column. |
| `mu_col` | Mutation-rate column. |
| `ac_col` | Allele-count column. |
| `initial_xi` | Initial optimizer value for xi. |
| `unique_mutation_limit` | Maximum unique-mutation count to evaluate. |
| `output_unique_mutation_file` | Output file for expected unique-mutation site counts. |
| `optimizer_method` | SciPy optimizer method. |
| `eps` | Numerical clipping value used for likelihood stability. |

## `smfs.estimate_xi_negative_binomial`

| Variable | Description |
| --- | --- |
| `input_dir` | Directory containing `input_file`. |
| `input_file` | Input site/SFS table. |
| `output_dir` | Output directory. |
| `output_xi_file` | Output file for the estimated xi value. |
| `save_xi_output` | Whether to write `output_xi_file`. |
| `sites_col` | Site-count column. |
| `mean_col` | Mean-theta column. |
| `sd_col` | Standard-deviation theta column. |
| `ac_col` | Allele-count column. |
| `initial_xi` | Initial optimizer value for xi. |
| `unique_mutation_limit` | Maximum unique-mutation count to evaluate. |
| `output_unique_mutation_file` | Output file for expected unique-mutation site counts. |
| `optimizer_method` | SciPy optimizer method. |
| `eps` | Numerical clipping value used for likelihood stability. |

## `smfs.estimate_smfs_from_unique_mutations`

| Variable | Description |
| --- | --- |
| `input_dir` | Directory containing `input_file`. |
| `input_file` | Observed allele-count table. |
| `unique_mutation_sites_dir` | Directory containing `unique_mutation_sites_file`. |
| `unique_mutation_sites_file` | Expected unique-mutation site-count distribution. |
| `output_dir` | Output directory. |
| `ac_col` | Allele-count column. |
| `sites_col` | Site-count column. |
| `unique_mutation_sites_col` | Expected unique-mutation site-count column. |
| `ac_limit` | Maximum allele count to evaluate. |
| `output_file` | Output SMFS probability file. |
| `apply_error_correction` | Whether to use the optional error correction. |
| `error_file` | Error table used when `apply_error_correction = true`. |
| `error_num_sites_col` | Error-table site-count column. |
| `error_rate_col` | Error-rate column. |
| `error_scale` | Scale factor applied to the error estimate. |

## `smfs.reconstruct_sfs_from_smfs`

| Variable | Description |
| --- | --- |
| `input_dir` | Directory containing `input_file`. |
| `input_file` | Observed site-frequency spectrum grouped by context. |
| `predicted_probability_dir` | Directory containing `predicted_probability_file`. |
| `predicted_probability_file` | Inferred SMFS probability table. |
| `xi_dir` | Directory containing `xi_file`. |
| `xi_file` | Estimated xi table. |
| `output_dir` | Output directory. |
| `output_file` | Output reconstructed SFS comparison table. |
| `ac_col` | Input allele-count column. |
| `sites_col` | Input site-count column. |
| `mean_col` | Input mean-theta column. |
| `sd_col` | Input standard-deviation theta column. |
| `methylation_level_col` | Input methylation-level column. |
| `predicted_probability_col` | Inferred SMFS probability column. |
| `xi_col` | Estimated xi column. |
| `ac_out_col` | Output allele-count column. |
| `methylation_level_out_col` | Output methylation-level column. |
| `mean_out_col` | Output mean-theta column. |
| `sd_out_col` | Output standard-deviation theta column. |
| `num_sites_out_col` | Output site-count column. |
| `prop_sites_out_col` | Internal proportional site-count column. |
| `type_out_col` | Output comparison-type column. |
| `data_type_label` | Output label for observed-data rows. |
| `theory_type_label` | Output label for reconstructed-SFS rows. |
| `ac_limit` | Maximum allele count included in the reconstruction. |
| `rows_per_context` | Number of contiguous input rows per context group. |

## `population_model.sample_site_mutation_rates`

| Variable | Description |
| --- | --- |
| `input_file` | Input mutation-rate metadata table. |
| `expanded_output_file` | Output expanded per-site sampled mutation-rate table. |
| `diagnostic_output_file` | Output diagnostic table. |
| `mean_theta_col` | Input mean-theta column. |
| `sd_theta_col` | Input standard-deviation theta column. |
| `scaled_mu_col` | Internal scaled-mu column. |
| `scaled_sd_col` | Internal scaled-standard-deviation column. |
| `gamma_shape_col` | Gamma shape column. |
| `gamma_scale_col` | Gamma scale column. |
| `seed_col` | Random seed column. |
| `num_sites_col` | Number-of-sites column. |
| `mu_col` | Input gnomAD mutation-rate column. |
| `methylation_level_col` | Input methylation-level column. |
| `sampled_mu_col` | Diagnostic sampled-mu column. |
| `sampled_sd_col` | Diagnostic sampled-standard-deviation column. |
| `diff_col` | Diagnostic relative-error column for mu. |
| `diff_sd_col` | Diagnostic relative-error column for standard deviation. |
| `mutation_number_out_col` | Output mutation-number column. |
| `mutation_rate_out_col` | Output mutation-rate column. |
| `mean_theta_out_col` | Output mean-theta column. |
| `sd_theta_out_col` | Output standard-deviation theta column. |
| `methylation_level_out_col` | Output methylation-level column. |
| `sampled_mu_out_col` | Output sampled-mu column. |
| `random_seed_max` | Upper bound for generated random seeds. |
| `diagnostic_sample_size` | Number of gamma samples used for diagnostics. |
| `show_diagnostic_plots` | Whether to show diagnostic plots. |
| `write_diagnostic_output` | Whether to write `diagnostic_output_file`. |

## `population_model.simulate_population_growth_mutations`

| Variable | Description |
| --- | --- |
| `input_dir` | Directory containing `input_file`. |
| `input_file` | Expanded per-site sampled mutation-rate table. |
| `output_dir` | Output directory. |
| `output_file` | Simulated mutation observations. |
| `mutation_rate_output_file` | Per-site mutation-rate table used by the simulation. |
| `site_col` | Site identifier column. |
| `mutation_number_col` | Mutation-number column. |
| `mutation_rate_col` | Mutation-rate column. |
| `mean_theta_col` | Mean-theta column. |
| `sd_theta_col` | Standard-deviation theta column. |
| `methylation_level_col` | Output methylation-level column. |
| `input_methylation_level_col` | Input methylation-level column. |
| `sampled_mu_col` | Sampled-mu column. |
| `allele_copies_col` | Output allele-copy column. |
| `mutation_id_col` | Output mutation-id column. |
| `initial_population` | Starting population size. |
| `final_population` | Final population size. |
| `time_scale_multiplier` | Multiplier used to compute the simulation time scale. |
| `normal_approximation_threshold` | Threshold where normal sampling replaces Poisson sampling. |
| `random_seed_modulus` | Modulus used when creating worker seeds. |
| `random_seed_range` | Random offset range used when creating worker seeds. |

## `population_model.sample_observed_mutations`

| Variable | Description |
| --- | --- |
| `input_dir` | Directory containing `input_file`. |
| `input_file` | Simulated mutation observations. |
| `output_dir` | Output directory. |
| `output_file` | Sampled observed-mutation output file. |
| `allele_copies_col` | Allele-copy column. |
| `sampled_mutation_col` | Output sampled-mutation column. |
| `final_population` | Final population size used by the hypergeometric sampling model. |
| `sample_size` | Sample size used by the hypergeometric sampling model. |

## `population_model.summarize_derived_alleles_by_site`

| Variable | Description |
| --- | --- |
| `input_dir` | Directory containing `input_file`. |
| `input_file` | Sampled mutation table. |
| `output_dir` | Output directory. |
| `output_file` | Per-site derived-allele summary. |
| `site_col` | Site identifier column. |
| `mutation_id_col` | Mutation-id column. |
| `methylation_level_col` | Methylation-level column. |
| `mutation_rate_col` | Mutation-rate column. |
| `mutation_number_col` | Mutation-number column. |
| `mean_theta_col` | Mean-theta column. |
| `sd_theta_col` | Standard-deviation theta column. |
| `sampled_mu_col` | Sampled-mu column. |
| `sampled_mutation_col` | Sampled-mutation column. |

## `population_model.fill_missing_simulation_sites`

| Variable | Description |
| --- | --- |
| `derived_allele_dir` | Directory containing `derived_allele_file`. |
| `derived_allele_file` | Per-site derived-allele summary. |
| `mutation_rate_dir` | Directory containing `mutation_rate_file`. |
| `mutation_rate_file` | Full mutation-rate table. |
| `output_dir` | Output directory. |
| `output_file` | Merged simulation table. |
| `site_col` | Site identifier column. |
| `sampled_mutation_sum_col` | Sampled mutation sum column. |
| `unique_mutations_col` | Unique mutations column. |
| `ac_col` | Output allele-count column. |

## `population_model.estimate_simulated_smfs_from_unique_mutations`

| Variable | Description |
| --- | --- |
| `input_dir` | Directory containing `input_file`. |
| `input_file` | Edited simulation data. |
| `unique_mutation_sites_dir` | Directory containing `unique_mutation_sites_file`. |
| `unique_mutation_sites_file` | Expected unique-mutation site-count distribution. |
| `output_dir` | Output directory. |
| `ac_col` | Allele-count column. |
| `unique_mutation_sites_col` | Expected unique-mutation site-count column. |
| `ac_limit` | Maximum allele count to evaluate. |
| `output_file` | Output simulated SMFS probability file. |
