# ria-genetic-susceptibility-phenotypes

This repository contains the analysis code and computational workflow accompanying:

**Conditional Phenotypic Relationships of Genetic Susceptibility in Ruptured Intracranial Aneurysm: Evidence from Mixed-Effects Regression and Additive Bayesian Networks**

**Authors:** [Authors]

**Journal:** [Journal]

**DOI:** [DOI, if available]

------------------------------------------------------------------------

## Overview

This repository contains the code used for data preparation, descriptive analyses, mixed-effects regression modelling, sensitivity analyses, and additive Bayesian network (ABN) modelling investigating the association between genetic susceptibility and phenotypic characteristics of ruptured intracranial aneurysms (RIA). The repository is intended to support **transparency, reproducibility, and reuse of the analytical workflow**.

As the original study dataset contains sensitive patient-level information, it cannot be publicly distributed. Instead, a synthetic dataset generated using the `synthpop` package is provided to demonstrate the structure of the data and enable execution of the analysis scripts.

------------------------------------------------------------------------

## Quick Start

Cloning the repository:

``` bash
git clone https://github.com/learbuehrer/ria-genetic-susceptibility-phenotypes.git
cd ria-genetic-susceptibility-phenotypes
```

Choose the workflow of interest:

-   Local analyses: run the scripts in `R/scripts/`:\
    Restore the R environment:

    ``` r
    install.packages("renv") 
    renv::restore()
    ```

    Run the analysis scripts in the following order:

    ``` text
    R/scripts/

    00_setup.R
    01_descriptives.R
    02_regressions.R
    03_sensitivity.R
    ```

-   HPC analyses (ABN/MCMC): see `hpc/README.md`

------------------------------------------------------------------------

## Repository structure

``` text
.
├── R/
│   ├── functions/
│   │   └── A_functions.R
│   └── scripts/
│       ├── 00_setup.R
│       ├── 01_descriptives.R
│       ├── 02_regressions.R
│       └── 03_sensitivity.R
│
├── hpc/
│   ├── functions/
│   │   └── B_mcmc_functions.R
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

Ia a synthetic dataset generated from the original study data using the R package `synthpop`.

The synthetic dataset preserves the overall structure and statistical characteristics of the original dataset while containing no real study participants.

The dataset is provided solely to:

-   demonstrate the analytical workflow,
-   enable execution of the code,
-   facilitate reproducibility of the computational analyses.

The synthetic dataset should **not** be used for scientific inference or interpretation of effect estimates.

### `figures/` and `tables/`

Figures and tables generated during the analysis are intentionally not tracked in this repository. They can be regenerated locally from the analysis scripts.

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

## High-performance computing

Computationally intensive analyses, including additive Bayesian network estimation using Markov Chain Monte Carlo (MCMC) methods, were executed on a high-performance computing environment.

The corresponding scripts are located in:

``` text
hpc/
```

These scripts reflect the original computational workflow and may require adaptation to local computing infrastructure, including:

-   file paths
-   scheduler settings
-   resource allocation
-   environment modules

Additional information is provided in:

``` text
hpc/README.md
```

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
