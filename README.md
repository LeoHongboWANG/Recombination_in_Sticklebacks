# Recombination rate and efficiency of linked selection in small and large stickleback populations

This repository contains the scripts used for linkage-map construction, recombination-rate estimation, crossover analyses, and statistical modelling for the manuscript:

**Recombination rate and efficiency of linked selection in small and large stickleback populations**

## Repository contents

| File | Description |
| --- | --- |
| `Lepmap3_pipeline.smk` | Snakemake workflow for running the Lep-MAP3 linkage-map pipeline. |
| `OrderMarkers.sh` | Shell script used to run Lep-MAP3 `OrderMarkers` and calculate recombination maps/rates. |
| `model.R` | R script for GLMM analyses, model selection, R2 calculation, and variance partitioning. |
| `sex_specific_recombination_CO.R` | R script for sex-specific recombination-rate and crossover-number analyses. |
| `window_recombination_diversity_analysis.R` | R script for window-based recombination, nucleotide diversity, CpG, gene density, repeat density, GWAS summary, and path analyses. |

## Data availability

The 10 sex-specific linkage maps generated in this study have been deposited in Zenodo:

**DOI:** [10.5281/zenodo.17197823](https://doi.org/10.5281/zenodo.17197823)

These linkage maps are used as input for downstream recombination-rate estimation and sex-specific recombination analyses.

