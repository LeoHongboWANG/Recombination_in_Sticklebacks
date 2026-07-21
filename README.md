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

## Dependencies

### Command-line tools

The linkage-map pipeline requires:

- `bash`
- `snakemake`
- `java`
- `Lep-MAP3`
The Lep-MAP3 path should be edited in the scripts if it differs from the local installation.

