###############################################################################-
#### File Content ####
# Script: 00_setup.R
#
# Purpose:
#   This script prepares the analysis environment and should be sourced
#   before running any other scripts in this repository.
#
# Tasks:
#   1. Load required R packages
#   2. Define color palettes
#   3. Load project-specific helper functions
#   4. Load the synthetic study dataset
#
# Input:
#   - data/synthetic/synthetic_data.rds
#   - R/functions/A_functions.R
#
# Output:
#   - Analysis environment
#   - dat_s object containing the synthetic dataset
#
# Notes:
#   - The repository contains synthetic data generated from the original
#     study cohort.
#   - Variables ending in "_Z" denote z-standardised variables:
#
#         Z = (x - mean(x)) / sd(x)
#
#   - All scripts assume that this setup file has been sourced first.
#
###############################################################################-

###############################################################################-
# Reproducibility ####
###############################################################################-
set.seed(0815)

###############################################################################-
## Required Packages ####
###############################################################################-
library(xtable)
library(tableone)
library(stringr)
library(kableExtra)
library(ggh4x)
library(biostatUZH)
library(ggplot2)
library(ggpubr)
library(purrr)
library(scales)
library(gt)
library(broom.mixed)
library(plotly)
library(performance)
library(ggcorrplot)
library(lme4)
library(lmerTest)
library(tidyr)
library(patchwork)

# # load as last package:
library(dplyr)

###############################################################################-
# Colour Palettes ####
###############################################################################-
darkColorPalette <- c("#B07C00", "#3A8BCD",
                      "#007F5F", "#BBAF00",  
                      "#005A8C", "#A30000",  
                      "#994C80", "grey50")


## color palette for the source sites:
sourceColorPalette <- c(
  "Barcelona"        = "#F6E7B4",
  "Geneva"           = "#F1D38A",
  "Oxford"           = "#E9B96E",
  "Helsinki"         = "#4C9A8A",
  "Kuopio"           = "#3B7F73",
  "ICAN"             = "#8A7F4F",
  "Utrecht"          = "#7FB3D5",
  "Dutch"            = "#4F81A8",
  "Poland"           = "#8A8F95",
  "GOSH"             = "#8C6FA8",
  "GERFHS"           = "#C06C84",
  "UCSF"             = "#A95D75")

###############################################################################-
# Project-Specific Functions ####
#
# Custom helper functions used throughout the repository.
###############################################################################-
source("R/functions/A_functions.R")

###############################################################################-
# Synthetic Study Dataset ####
#
# The dataset contains synthetic participant-level observations that
# preserve the overall structure and statistical properties of the
# original study data while containing no real participants.
###############################################################################-

dat_s <- readRDS("data/synthetic/synthetic_data.rds")

###############################################################################-
# Quick Dataset Overview
###############################################################################-

message(
  "Loaded synthetic dataset with ",
  nrow(dat_s),
  " observations and ",
  ncol(dat_s),
  " variables.")


###############################################################################
# Exposure Variables
## Needed for all analysis, so predefine the default here.
###############################################################################

exp_variables <- c(
  "PRStotal_Z",
  "Sex")