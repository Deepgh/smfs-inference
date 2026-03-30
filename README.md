# Allele Mutation Spectrum Modeling and Recurrent Mutation Estimation

This repository contains Python scripts used to analyze the human site frequency spectrum (SFS) to infer the probability of recurrent mutations and mutation-specific rates in population genomic data. The work contributes to understanding how often mutations at the same site arise independently, a critical factor in interpreting rare variant patterns, molecular evolution, and disease genetics.

## Overview

The code estimates the expected number and distribution of unique mutations at each allele count using two approaches:

* **Maximum likelihood estimation (MLE)** to fit a Poisson or negative binomial model of mutation recurrence.
* **Convolution-based modeling** to account for recurrence due to shared site-specific mutation rates.

It processes simulation-based datasets and applies computational models to infer:

* Per-site recurrence probabilities.
* Mutation-specific parameters (`μ`, `θ`, etc.).
* Xi (`ξ`): the effective mutation flux or compound rate parameter of interest.
* Adjusted SFS prediction curves for publication or regression testing.

## Project Structure

* `matching_mut_rate_simple.py`: Aligns mutation rate and methylation-level datasets to simulation-based SFS entries by mutation-context triplet.
* `estimate_lambda_mew.py` / `estimate_lambda_mew_sd.py`: Estimate the global mutation flux parameter (`ξ`) using observed segregating/non-segregating counts.
* `unq_mut_sites_per_mut_mew.py`: Predicts the expected number of sites with a given number of unique mutations per mutation type using Poisson processes.
* `unq_mut_sites_per_mut_mew_sd.py`: Extends the above with variance modeling using negative binomial distributions.
* `est_smfs_num.py`: Uses convolution to numerically estimate the expected SFS under the inferred model, including simulation tables where rows can represent multiple sites.
* `est_smfs_log_lik.py`: Applies likelihood maximization to derive the predicted SFS and estimate the probability of recurrent mutation events.
* `scaled_mu_from_musd.py`: Samples site-level mutation rates from mutation-type-level gamma summaries.
* `pop_growth_sampled_mu_diff_sites.py` / `sampled_muts.py`: Simulate recurrent mutation growth and downsample to an observed cohort.
* `der_allele_sum_per_site.py` / `input_miss_rows.py`: Collapse simulated mutations back to one row per site and reinsert zero-count sites.

## Input Data

The analysis relies on:

* Simulation outputs mimicking human genome variation under varying mutation rates.
* Tabulated mutation-specific metadata including sampled mutation rates, methylation level, and context.

> **Note**: The original data files are not included here due to size constraints. You may contact the authors or refer to the associated publication to access the input data, or use mock datasets for testing.

## Usage

These scripts are standalone and may be run with Python 3.8+ and the following libraries:

```bash
pip install numpy pandas scipy matplotlib
```

Each script can be run as a top-level program. Paths to data files must be updated or mocked as appropriate for your environment.

Example (estimating xi):

```bash
python estimate_lambda_mew.py
```

Example (inferring SMFS numerically):

```bash
python est_smfs_num.py
```

Example (MLE approach to SFS fitting):

```bash
python est_smfs_log_lik.py
```

## Output

* CSV tables of predicted allele count frequencies.
* Tabulated xi estimates per dataset.
* Intermediate files with predicted number of unique mutations per count.

These are useful for:

* Regression testing simulation pipelines.
* Comparing observed versus expected SFS.
* Plotting model fits in publication figures.

## Citation

If you use this code for your research, please cite:

> \[Author Names]. *Estimating recurrent mutation probabilities in allele frequency spectra from human population data.* \[Journal], \[Year], \[DOI].

## License

This code is released under the MIT License.
