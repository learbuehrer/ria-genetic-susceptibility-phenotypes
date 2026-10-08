# ria-genetic-susceptibility-phenotypes

This repository contains the analysis code and computational workflow accompanying:

**Conditional Phenotypic Relationships of Genetic Susceptibility in Ruptured Intracranial Aneurysm: Evidence from Mixed-Effects Regression and Additive Bayesian Networks**

**Authors:** [Authors]

**Journal:** [Journal]

**DOI:** [DOI, if available]

------------------------------------------------------------------------

## Overview

This repository contains the code used for data processing, descriptive analyses, statistical modelling, and computational analyses reported in the publication. The repository is intended to support **transparency, reproducibility, and reuse of the analytical workflow**. Patient-level study data are not included in this repository due to data protection and/or data-sharing restrictions. A synthetic dataset is provided to illustrate the structure of the analysis data and to facilitate testing of the analysis code.

------------------------------------------------------------------------

## Repository structure

``` text
.
├── R/
│   ├── functions/
│   │   └── 00_functions.R
│   └── scripts/
│       ├── 01_descriptives.R
│       ├── 02_regressions.R
│       └── 03_sensitivity.R
│
├── hpc/
│   ├── functions/
│   │   └── 00_mcmc_functions.R
│   ├── scripts/
│       ├── 01_RIA_me.R       
│   │   └── 02_RIA_me_nonFins.R
│   ├── load_env.sh
│   ├── submit.sh
│   ├── renv.lock
│   └── README.md
│
├── data/
│   └── synthetic/
│       ├── synthetic_data.rds
│       └── README.md
│
├── figures/          # generated locally; not tracked by Git
├── tables/           # generated locally; not tracked by Git
│
│
├── renv.lock
├── .gitignore
└── README.md
```

### `R/`

Contains the R code used for the analyses.

-   `functions/` contains reusable functions for data processing, statistical modelling, and visualisation.
-   `analysis/` contains the analysis scripts executed for the publication.

### `hpc/`

Contains scripts required to run computationally intensive analyses on the high-performance computing (HPC) environment. The `submit.sh` script provides the corresponding job submission workflow.

### `data/synthetic/`

Contains a synthetic dataset with the same general structure as the analysis data. The synthetic data **do not contain patient-level information** and cannot be used to reproduce the reported results.

### `figures/` and `tables/`

Figures and tables generated during the analysis are intentionally **not tracked in this repository**. They can be regenerated locally from the analysis scripts.

------------------------------------------------------------------------

## Software environment

The analyses were conducted in **R**.

Package versions are managed using [`renv`](https://rstudio.github.io/renv/). The `renv.lock` file records the R package versions used for the analysis. To initialise the project environment after cloning the repository:

``` r
install.packages("renv")
renv::restore()
```

This installs the package versions recorded in `renv.lock`. To inspect the current project environment:

``` r
renv::status()
```

------------------------------------------------------------------------

## Reproducing the analyses

### 1. Clone the repository

``` bash
git clone [REPOSITORY-URL]
cd [REPOSITORY-NAME]
```

### 2. Open the R project

Open:

``` text
[REPOSITORY-NAME].Rproj
```

in RStudio.

### 3. Restore the R environment

``` r
renv::restore()
```

### 4. Run the analyses

The analysis scripts are located in:

``` text
R/scripts/
```

and are intended to be executed in the indicated order. The scripts use the study analysis data, which are not distributed with this repository. The synthetic dataset can be used to inspect and test the analysis workflow where applicable.

------------------------------------------------------------------------

## High-performance computing

Computationally intensive analyses were performed using a high-performance computing environment.

The corresponding scripts are located in:

``` text
hpc/
```

The hpc/ directory contains template scripts used for the original HPC execution environment. These scripts illustrate the computational workflow but may require adaptation to local infrastructure.

The main submission script is:

``` text
hpc/submit.sh
```

HPC-specific instructions are provided in:

``` text
hpc/README.md
```

The HPC workflow may require adaptation of paths, resource requests, environment modules, or scheduler settings to the local computing infrastructure.

------------------------------------------------------------------------

## Data availability

The underlying study data are not publicly available because they contain sensitive patient-level information and/or are subject to data-sharing restrictions. The synthetic dataset included in this repository is provided solely for demonstration and testing purposes. It does not contain real patient-level observations and cannot be used to reproduce the numerical results reported in the publication.

------------------------------------------------------------------------

## Citation

If you use this repository or its code, please cite:

> [Full citation of the publication]

------------------------------------------------------------------------

## License

This project is licensed under the [MIT License](LICENSE).

Copyright (c) [YEAR] [AUTHOR / COPYRIGHT HOLDER]
