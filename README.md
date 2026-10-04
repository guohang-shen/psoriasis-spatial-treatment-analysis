# Psoriasis spatial and treatment-convergence analysis

This repository contains the analysis and figure-generation code for the psoriasis manuscript.

## Scope

The project evaluates locked psoriasis gene programs across two spatial platforms, treatment cohorts, an independent etanercept time course, and a randomized placebo-controlled trial. The primary analyses use patient-level summaries, matched-gene-set controls, signed Stouffer meta-analysis, and Benjamini–Hochberg correction.

## Data access

Raw expression matrices are not redistributed here. They are available from GEO:

- GSE206391
- GSE314158
- GSE117239
- GSE85034
- GSE136757
- GSE183047
- GSE278330
- GSE228421
- GSE11903

## Repository contents

- `scripts/`: core analysis and figure-generation scripts
- `derived_tables/`: selected derived tables used in the manuscript
- `analysis_environment_R4.6.0.txt`: R session information
- `figure_manifest.tsv`: figure-to-file mapping
- `CODE_AVAILABILITY.md`: release and reproducibility notes

## Reproducibility

Run scripts from the project root with R 4.6.0. Input files and accession-specific preparation steps are described in the script headers. The treatment convergence analysis uses arm-level signed changes and weighted Stouffer statistics with weights proportional to the square root of the contributing patient count. GSE117239 ustekinumab and etanercept are treated as two arms from one cohort.

This repository is a versioned code release for the manuscript. Raw GEO data remain in their original public repositories.
