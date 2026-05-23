# Accounting for Recurrent Mutation in the Frequency Spectrum of Rare Alleles

This repository contains Python scripts associated with the manuscript *Accounting for recurrent mutation in the frequency spectrum of rare alleles*. The code analyzes the human site frequency spectrum (SFS), estimates the single mutation frequency spectrum (SMFS), models recurrent mutation, and simulates allele-count data under population growth. The code is organized as two independent script collections:

* `smfs/`: scripts for matching mutation-rate metadata, estimating xi, estimating expected numbers of unique mutations, and deriving predicted SMFS curves.
* `population_model/`: scripts for generating and processing simulated population-genetic data used by the SMFS analyses.

The scripts can be configured with TOML files passed through `--config`. `configs/template.toml` is a user-editable starting point, and `configs/paper_config.toml` records the paths and settings used for the submitted manuscript.

## Overview

The code supports two related workflows:

* Simulation of mutation counts across sites under a growing population model.
* Estimation of recurrent mutation quantities and predicted site-frequency spectra from simulation or observed-data summaries.

The analysis uses tabular input data containing mutation rates, methylation levels, allele counts, site counts, and mutation-specific uncertainty parameters such as mean and standard deviation of theta.

The manuscript defines the SMFS as the frequency spectrum of alleles descending from a single mutational event. Under recurrent mutation, the SFS can be represented as a weighted sum of convolutions of the SMFS with itself, separating recurrent mutation from the population-genetic processes that shape the SMFS.

## File Descriptions

### `smfs/`

* `smfs/match_sfs_to_mutation_rates.py`: Matches SFS rows to mutation-rate metadata by reference context, alternate context, and methylation level. It also attempts a reverse-complement context match for rows that do not match directly, then merges matched mutation-rate fields such as `mu_gnomad`, `mean_theta`, `sd_theta`, and error estimates into the SFS data.

* `smfs/estimate_xi_poisson.py`: Estimates a global xi parameter using a Poisson model with site-specific sampled mutation rates. It separates segregating and non-segregating sites, minimizes the negative log-likelihood for xi, writes the xi estimate, and computes the expected number of sites with each possible number of unique mutations.

* `smfs/estimate_xi_negative_binomial.py`: Estimates xi while accounting for uncertainty in mutation-rate/theta estimates using a negative-binomial marginal model. It uses per-site mean and standard deviation columns, estimates xi by likelihood optimization, and writes expected counts of sites with each number of unique mutations.

* `smfs/estimate_smfs_from_unique_mutations.py`: Uses the expected unique-mutation count distribution and observed allele-count data to numerically derive predicted SMFS probabilities up to a specified allele-count limit. This version is configured for data where allele-count observations are stored with site-count weights.

### `population_model/`

* `population_model/sample_site_mutation_rates.py`: Reads mutation-rate metadata with mean and standard deviation values, converts those values into gamma-distribution parameters, samples per-site mutation rates, and writes an expanded table with one sampled mutation-rate row per site. It also includes diagnostic plots comparing sampled and input values.

* `population_model/simulate_population_growth_mutations.py`: Simulates mutation accumulation through a growing population using sampled site-specific mutation rates. It outputs both simulated derived-allele observations and the mutation-rate table used for each simulated site.

* `population_model/sample_observed_mutations.py`: Samples observed mutations from simulated allele-copy counts using a hypergeometric model. It filters out mutations not observed in the sample and writes the sampled mutation data to a compressed CSV file.

* `population_model/summarize_derived_alleles_by_site.py`: Aggregates sampled mutation records by site. For each site, it computes the total sampled derived-allele count, the number of unique mutations, and carries forward metadata such as mutation rate, methylation level, sampled mu, mean theta, and standard deviation of theta.

* `population_model/fill_missing_simulation_sites.py`: Merges per-site derived-allele summaries back into the full mutation-rate table so sites with no observed sampled mutations are retained. Missing allele-count and unique-mutation values are filled with zero, producing an edited simulation dataset for downstream SMFS estimation.

* `population_model/estimate_simulated_smfs_from_unique_mutations.py`: Numerically estimates predicted SMFS probabilities from simulated data and an expected unique-mutation count distribution. This script is similar in purpose to `smfs/estimate_smfs_from_unique_mutations.py`, but is configured for simulation data where allele counts are counted by rows rather than by a separate site-count column.

## Typical Workflow

Copy the template config and edit paths for your environment:

```bash
cp configs/template.toml configs/my_analysis.toml
```

Each script reads only its own section of the config file. For example, `smfs/estimate_xi_poisson.py` reads the `[smfs.estimate_xi_poisson]` section:

```toml
[smfs.estimate_xi_poisson]
input_dir = "results"
input_file = "simulation_sites_with_observed_counts.csv.gz"
output_dir = "results"
output_xi_file = "xi_poisson.csv"
initial_xi = 1000.0
unique_mutation_limit = 300
```

Run scripts with the edited config:

```bash
python smfs/estimate_xi_poisson.py --config configs/my_analysis.toml
```

If `--config` is omitted, scripts fall back to the constants defined near the top of each script. See `docs/config_variables.md` for all supported config variables and `docs/input_files.md` for required input columns.

A typical simulation workflow is:

```bash
python population_model/sample_site_mutation_rates.py --config configs/my_analysis.toml
python population_model/simulate_population_growth_mutations.py --config configs/my_analysis.toml
python population_model/sample_observed_mutations.py --config configs/my_analysis.toml
python population_model/summarize_derived_alleles_by_site.py --config configs/my_analysis.toml
python population_model/fill_missing_simulation_sites.py --config configs/my_analysis.toml
```

A typical SMFS estimation workflow is:

```bash
python smfs/estimate_xi_poisson.py --config configs/my_analysis.toml
python smfs/estimate_xi_negative_binomial.py --config configs/my_analysis.toml
python smfs/estimate_smfs_from_unique_mutations.py --config configs/my_analysis.toml
```

Use `smfs/match_sfs_to_mutation_rates.py` when preparing SFS data that needs to be matched to mutation-rate metadata before estimation.

## Input Data

The analysis relies on:

* Simulation outputs from the population-model scripts.
* Tabulated mutation-specific metadata, including mutation rates, methylation level, sequence context, mean theta, and standard deviation of theta.
* SFS or allele-count summaries with site-count information.

The original data files are not included in this repository due to size and project-specific storage constraints.

See `docs/input_files.md` for the expected input columns and generated outputs for each script.
See `docs/config_variables.md` for the available TOML config variables.

## Reproducibility Notes

The original input datasets are not distributed with this repository. To reproduce an analysis, provide equivalent input tables with the columns described in `docs/input_files.md` and configure paths with a TOML file such as `configs/template.toml` or `configs/paper_config.toml`.

Several population-model scripts use random sampling. Re-running those scripts may produce different simulated outputs unless random seeds and input files are controlled. Default generated tables use normalized snake_case column names; use `configs/paper_config.toml` when working with the original publication data layout.

## Requirements

The scripts use Python 3 and the following third-party libraries:

```bash
pip install numpy pandas scipy matplotlib
```

The project dependencies are also listed in `pyproject.toml`.

## Output

Depending on which scripts are run, outputs include:

* Simulated mutation-rate tables.
* Sampled mutation records.
* Per-site derived-allele summaries.
* Xi estimates.
* Expected counts of sites with each number of unique mutations.
* Predicted SMFS probability tables.
* Diagnostic plots for sampled mutation-rate checks.

## Citation

If you use this code for your research, please cite:

> Deepjyoti Ghosh, Kyle Williams, Joshua Schraiber, and Yuval Simons. *Accounting for recurrent mutation in the frequency spectrum of rare alleles.* Manuscript submitted.

This repository also includes `CITATION.cff` for software citation metadata.

## License

This code is released under the MIT License.
