# Allele Mutation Spectrum Modeling and Recurrent Mutation Estimation

This repository contains Python scripts used to analyze the human site frequency spectrum (SFS), estimate recurrent mutation probabilities, and simulate allele-count data under population growth. The code is organized as two independent script collections:

* `smfs/`: scripts for matching mutation-rate metadata, estimating lambda, estimating expected numbers of unique mutations, and deriving predicted SMFS curves.
* `population_model/`: scripts for generating and processing simulated population-genetic data used by the SMFS analyses.

The scripts currently use hard-coded input and output paths from the original analysis environment. Before running them in a new environment, update the path variables near the top of each script.

## Overview

The code supports two related workflows:

* Simulation of mutation counts across sites under a growing population model.
* Estimation of recurrent mutation quantities and predicted site-frequency spectra from simulation or observed-data summaries.

The analysis uses tabular input data containing mutation rates, methylation levels, allele counts, site counts, and mutation-specific uncertainty parameters such as mean and standard deviation of theta.

## File Descriptions

### `smfs/`

* `smfs/match_sfs_to_mutation_rates.py`: Matches SFS rows to mutation-rate metadata by reference context, alternate context, and methylation level. It also attempts a reverse-complement context match for rows that do not match directly, then merges matched mutation-rate fields such as `mu_gnomad`, `mean_theta`, `sd_theta`, and error estimates into the SFS data.

* `smfs/estimate_lambda_poisson.py`: Estimates a global lambda parameter using a Poisson model with site-specific sampled mutation rates. It separates segregating and non-segregating sites, minimizes the negative log-likelihood for lambda, writes the lambda estimate, and computes the expected number of sites with each possible number of unique mutations.

* `smfs/estimate_lambda_negative_binomial.py`: Estimates lambda while accounting for uncertainty in mutation-rate/theta estimates using a negative-binomial marginal model. It uses per-site mean and standard deviation columns, estimates lambda by likelihood optimization, and writes expected counts of sites with each number of unique mutations.

* `smfs/estimate_smfs_from_unique_mutations.py`: Uses the expected unique-mutation count distribution and observed allele-count data to numerically derive predicted SMFS probabilities up to a specified allele-count limit. This version is configured for data where allele-count observations are stored with site-count weights.

### `population_model/`

* `population_model/sample_site_mutation_rates.py`: Reads mutation-rate metadata with mean and standard deviation values, converts those values into gamma-distribution parameters, samples per-site mutation rates, and writes an expanded table with one sampled mutation-rate row per site. It also includes diagnostic plots comparing sampled and input values.

* `population_model/simulate_population_growth_mutations.py`: Simulates mutation accumulation through a growing population using sampled site-specific mutation rates. It outputs both simulated derived-allele observations and the mutation-rate table used for each simulated site.

* `population_model/sample_observed_mutations.py`: Samples observed mutations from simulated allele-copy counts using a hypergeometric model. It filters out mutations not observed in the sample and writes the sampled mutation data to a compressed CSV file.

* `population_model/summarize_derived_alleles_by_site.py`: Aggregates sampled mutation records by site. For each site, it computes the total sampled derived-allele count, the number of unique mutations, and carries forward metadata such as mutation rate, methylation level, sampled mu, mean theta, and standard deviation of theta.

* `population_model/fill_missing_simulation_sites.py`: Merges per-site derived-allele summaries back into the full mutation-rate table so sites with no observed sampled mutations are retained. Missing allele-count and unique-mutation values are filled with zero, producing an edited simulation dataset for downstream SMFS estimation.

* `population_model/estimate_simulated_smfs_from_unique_mutations.py`: Numerically estimates predicted SMFS probabilities from simulated data and an expected unique-mutation count distribution. This script is similar in purpose to `smfs/estimate_smfs_from_unique_mutations.py`, but is configured for simulation data where allele counts are counted by rows rather than by a separate site-count column.

## Typical Workflow

A typical simulation workflow is:

```bash
python population_model/sample_site_mutation_rates.py
python population_model/simulate_population_growth_mutations.py
python population_model/sample_observed_mutations.py
python population_model/summarize_derived_alleles_by_site.py
python population_model/fill_missing_simulation_sites.py
```

A typical SMFS estimation workflow is:

```bash
python smfs/estimate_lambda_poisson.py
python smfs/estimate_lambda_negative_binomial.py
python smfs/estimate_smfs_from_unique_mutations.py
```

Use `smfs/match_sfs_to_mutation_rates.py` when preparing SFS data that needs to be matched to mutation-rate metadata before estimation.

## Input Data

The analysis relies on:

* Simulation outputs from the population-model scripts.
* Tabulated mutation-specific metadata, including mutation rates, methylation level, sequence context, mean theta, and standard deviation of theta.
* SFS or allele-count summaries with site-count information.

The original data files are not included in this repository due to size and project-specific storage constraints.

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
* Lambda estimates.
* Expected counts of sites with each number of unique mutations.
* Predicted SMFS probability tables.
* Diagnostic plots for sampled mutation-rate checks.

## Citation

If you use this code for your research, please cite:

> [Author Names]. *Estimating recurrent mutation probabilities in allele frequency spectra from human population data.* [Journal], [Year], [DOI].

## License

This code is released under the MIT License.
