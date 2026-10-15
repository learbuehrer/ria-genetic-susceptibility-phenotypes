################################################################################-
# File to run on the HPC ####
# Script: 01_RIA_me.R
#
# Purpose:
#   Runs additive Bayesian network analyses and MCMC validation on a
#   SLURM-based HPC environment.
#
# Input:
#   - data/dat_s.rds
#   - Banmatrix_ISGC_genetic_scores.xlsx
#
# Output:
#   - Bayesian network models
#   - MCMC diagnostics
#   - DAG figures
#   - Validation reports
#   - Run metadata
#
################################################################################-
## ============================ -
## User-defined run parameters ####
## ============================ -
run_name <- Sys.getenv("RUN_NAME")
if (run_name == "") {
  stop("RUN_NAME environment variable is not set. Did you submit via submit.sh?")
}

base_dir  <- "results"                  # output directory for figures


## Create run-specific output directory
out_dir <- file.path(base_dir, run_name)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

## load the validation parameters from environment variables
## Validation parameters (from environment)
MCMC.SEEDS <- strsplit(Sys.getenv("MCMC_SEEDS"), ",")[[1]]
MCMC.SCHEME <- as.numeric(
  strsplit(Sys.getenv("MCMC_SCHEME"), ",")[[1]])
burnin.length <- as.integer(Sys.getenv("BURNIN_LENGTH"))
thinningsteps <- as.integer(Sys.getenv("THINNING_STEPS"))

## set the seed for reproducible results: 
set.seed(0815)

# load the libraries
library(abn)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

## libraries for the MCMC validation:
library(foreach)
library(doParallel)
library(mcmcabn)
library(coda)


n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))

if (.Platform$OS.type == "unix") {
  cl <- parallel::makeForkCluster(n_cores)
} else {
  cl <- parallel::makeCluster(n_cores)
}

doParallel::registerDoParallel(cl)


## ============================ -
## Source MCMC validation functions ####
## ============================ -

source("hpc/scripts/B_mcmc_functions.R")


# load the data:
dat <- readRDS("data/synthetic/synthetic_data.rds")


################################################################################-
# ABN Setup ####
# Alternative analysis scenarios
#
# The HPC workflow can be executed under various analysis scenarios:
#
# Main analysis: complete cohort including all contributing sources, 
# using the overall PRS or AC together with the first five ancestry principal components.
#
# Non-Finnish sensitivity analysis: identical model specification but excluding the Finnish cohorts 
# (Helsinki and Kuopio) to evaluate the influence of the Finnish genetic background on network structure.
#
# PRS-cluster analysis: complete cohort using genetic cluster-specific scores (e.g. PRS2) 
# instead of the overall PRS to assess whether cluster-specific genetic architectures 
# lead to different network structures.
#
# To run a different scenario, users only need to modify the variable selection (keep), 
# the corresponding distribution specification (dists), 
# and, where applicable, the cohort filtering step prior to Bayesian network estimation (line 132).
################################################################################-

# select variables to include in the BN:
keep <- c("Sex",
          "Age.at.time.of.SAH_Z",
          "Positive.famillial.history",
          "Hypertension",
          "Smoking_Current_Former_No",
          "IAlocation_cow",
          "IA_size.at.time.of.SAH.log_Z",
          "Multiple.IAs",
          "PRStotal_Z",
          #"allele_count_Z",
          #"PRS1_Z",
          #"PRS2_Z",
          #"PRS3_Z",
          #"PRS4_Z",
          #"PRS5_Z",
          "PC1_AVG_Z",
          "PC2_AVG_Z",
          "PC3_AVG_Z",
          "PC4_AVG_Z",
          "PC5_AVG_Z",
          "Source")


datBN <- dat %>%
  #filter(!(Source %in% c("Kuopio", "Helsinki"))) %>%  ## remove the Fins by filtering the data
  select(all_of(keep)) %>% 
  na.omit() %>% 
  mutate(Source = droplevels(as.factor(Source))) 

datBN_red <- datBN %>%
  select(-Source)


# distributions
dists <- list(
  Sex ="binomial",
  Age.at.time.of.SAH_Z ="gaussian",  
  Positive.famillial.history ="binomial",
  Hypertension ="binomial",
  Smoking_Current_Former_No ="multinomial",
  IAlocation_cow ="multinomial",
  IA_size.at.time.of.SAH.log_Z ="gaussian",
  Multiple.IAs ="binomial",
  PRStotal_Z ="gaussian",
  #allele_count_Z ="gaussian",
  #PRS1_Z ="gaussian",
  PRS2_Z ="gaussian",
  #PRS3_Z ="gaussian",
  #PRS4_Z ="gaussian",
  #PRS5_Z ="gaussian",
  PC1_AVG_Z ="gaussian",
  PC2_AVG_Z ="gaussian",
  PC3_AVG_Z ="gaussian",
  PC4_AVG_Z ="gaussian",
  PC5_AVG_Z ="gaussian")

# banmatrix:
keep <- keep[keep != "Source"] # remove Source from ban matrix
allowmat <- readxl::read_excel("data/Banmatrix_ISGC_genetic_scores.xlsx",
                               sheet="allowmatrix_RIA") %>% 
  tibble::column_to_rownames('ChildParent') %>%
  as.matrix() 

## creat banmat: 1-allowmat
banmat <- 1-allowmat

# only use ban matrix for variables fitted:
banmat <- banmat[c(keep),
                 c(keep)]


# retainmatrix
retainmat <- readxl::read_excel("data/synthetic/Banmatrix_ISGC_genetic_scores.xlsx",
                                sheet="retainmatrix_RIA") %>% 
  tibble::column_to_rownames('ChildParent') %>%
  as.matrix() 

# only use retainmat for variables fitted:
retainmat <- retainmat[c(keep),
                       c(keep)]


## ============================ -
## Run meta-information ####
## ============================ -
run_meta <- list(
  run_name   = run_name,
  start_time = Sys.time(),
  steps      = character(),
  settings   = list(
    method       = "mle",
    score        = "bic",
    max.parents  = ncol(datBN) - 1,
    centre       = TRUE,
    validation   = list(
      type              = "MCMC",
      mcmc_seeds        = MCMC.SEEDS,
      mcmc_scheme       = MCMC.SCHEME,
      mcmc_iterations   = MCMC.SCHEME[1],
      burnin_length     = burnin.length,
      thinning_steps   = thinningsteps
    )
  ),
  data_info = list(
    n_obs  = nrow(datBN),
    n_vars = ncol(datBN),
    vars   = colnames(datBN)
  ),
  slurm = list(
    job_id    = Sys.getenv("SLURM_JOB_ID"),
    job_name  = Sys.getenv("SLURM_JOB_NAME"),
    node_list = Sys.getenv("SLURM_NODELIST"),
    n_tasks   = Sys.getenv("SLURM_NTASKS"),
    cpus_task = Sys.getenv("SLURM_CPUS_PER_TASK")
  )
)



run_meta$start_time <- Sys.time()

log_step <- function(meta, msg) {
  meta$steps <- c(
    meta$steps,
    paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  )
  meta
}

## ============================ -
## Fit of the BN ####
## ============================ -

cacheMLE <- buildScoreCache(
  data.df     = datBN,
  data.dists  = dists,
  method      = "mle",
  centre      = TRUE,
  max.parents = length(datBN)-1,
  dag.banned  = banmat,
  dag.retained = retainmat,
  group.var   = "Source"
)

run_meta <- log_step(run_meta, "Score cache built")

dagMP <- mostProbable(
  score.cache = cacheMLE,
  score       = "bic"
)
run_meta <- log_step(run_meta, "Most probable DAG identified")

## ---------------------------- -
## Plot: Most Probable DAG ####
## ---------------------------- -

pdf(
  file   = file.path(out_dir, paste0(run_name, "_dagMP.pdf")),
  width  = 8,
  height = 8
)
plot(dagMP)
dev.off()

run_meta <- log_step(run_meta, "Most probable DAG plotted")

## ============================ -
## Fit the DAG ####
## ============================ -

fit <- fitAbn(dagMP, method = "mle", group.var = "Source")
fit_file <- file.path(out_dir, paste0(run_name, "_fit.txt"))

sink(fit_file)
print(fit)
sink()

run_meta <- log_step(run_meta, "DAG fitted and fit saved")

## ============================ -
## Save the resulting model to JSON format ####
## ============================ -
json_export <- export_abnFit(
  fit,
  file = file.path(out_dir, paste0(run_name, "_fit.json")),
  scenario_id = run_name,
  label = paste0("ABN model - ", run_name))


## ============================ -
## Convert to bnlearn DAG ####
## ============================ -
nodes <- colnames(dagMP$dag)
dag <- bnlearn::empty.graph(nodes)
bnlearn::amat(dag) <- t(dagMP$dag)

run_meta <- log_step(run_meta, "Converted DAG to bnlearn object")

## ============================ -
## CPDAG construction  ####
## ============================ -

cons.cpdag <- bnlearn::cpdag(dag)
cons.amat  <- bnlearn::amat(cons.cpdag)

## Remove banned edges
cons.amat[cons.amat == 1 & t(banmat) == 1] <- 0

cpdag.cons.amat <- bnlearn::empty.graph(nodes)
bnlearn::amat(cpdag.cons.amat) <- cons.amat

run_meta <- log_step(run_meta, "CPDAG constructed and constrained")

## ============================ -
## DAG and CPDAG side-by-side  ####
## ============================ -

pdf(
  file   = file.path(out_dir, paste0(run_name, "_cpdag.pdf")),
  width  = 14,
  height = 7
)

par(mfrow = c(1, 2))
plot(dag, main = "DAG")
plot(cpdag.cons.amat, main = "CPDAG II")
par(mfrow = c(1, 1))
dev.off()

run_meta <- log_step(run_meta, "DAG and CPDAG plotted")

## ============================ -
# Validation: MCMC runs  ####
## ============================ -

run_meta <- log_step(run_meta, "Starting MCMC-based validation")

## Output directory for validation
validation_dir <- file.path(out_dir, "mcmc")
if (!dir.exists(validation_dir)) {
  dir.create(validation_dir, recursive = TRUE)
}


## Store validation paths in metadata
run_meta$validation <- list(
  output_dir = validation_dir
)

## ---------------------------- -
## Max parent search  ####
## ---------------------------- -

modelmaxparents <- findMaxParents(
  dat    = datBN,
  dists  = dists,
  banmat = banmat,
  retainmat = retainmat,
  SCORE  = "bic",
  PRIOR  = 1,
  group.var   = "Source"
)


run_meta <- log_step(run_meta, "Max-parent model search completed")

## ---------------------------- -
## Score plot  ####
## ---------------------------- -

ScorePlotMaxParent(
  modelmaxparents$net.scores,
  FILENAMEbase = paste0(validation_dir, "/")
)

run_meta <- log_step(run_meta, "Score-vs-parents plot created")

## ---------------------------- -
## DAG from max-parent models
## ---------------------------- -

DAGMaxParent(
  net.scores.dags = modelmaxparents$net.scores.dags,
  net.scores      = modelmaxparents$net.scores,
  par1            = 4,
  par2            = 4,
  FILENAMEbase    = paste0(validation_dir, "/"),
  filename        = "DAGMaxParents",
  retainmat       = retainmat)

run_meta <- log_step(run_meta, "Max-parent DAGs generated")

## ---------------------------- -
## Run MCMC ####
## ---------------------------- -
# record MCMC start time
run_meta <- log_step(run_meta, "MCMC start")

run_meta$mcmc_start_time <- Sys.time()

# run MCMC
modelMCMC <- mcmcmABN(
  dat = datBN,
  dists = dists,
  banmat = banmat,
  retainmat = retainmat,
  FILENAMEbase = out_dir,
  FILENAME = run_name,
  max.par = modelmaxparents$max.par,
  MCMC.SEEDS = MCMC.SEEDS,
  MCMC.SCHEME = MCMC.SCHEME,
  PROB.REV = 0,
  PROB.MBR = 0,
  MCMC.PRIOR = 1,
  SCORE = "bic",
  group.var   = "Source"
)

# record MCMC end time
run_meta$mcmc_end_time <- Sys.time()

run_meta <- log_step(run_meta, "MCMC finished")

# save to validation dir
saveRDS(
  modelMCMC,
  file = file.path(validation_dir, paste0("modelMCMC_", run_name, ".rds"))
)

## ============================ -
## Shutdown parallel backend  ####
## ============================ -

parallel::stopCluster(cl)

## ---------------------------- -
## Postprocessing of the MCMC results  ####
## ---------------------------- -

# record start of postprocessing
run_meta <- log_step(run_meta, "Postprocessing MCMC started")
run_meta$postprocessing_start <- Sys.time()

# ----------------------------------- -
# Burn-in  ####
# ----------------------------------- -
mcmc.out.list.burn <- postBURNin(
  mcmc.out.list = modelMCMC$mcmc.out.list,
  burnin.length = burnin.length
)

# Thinning
mcmc.out.list.thin <- postTHINN(
  mcmc.out.list = modelMCMC$mcmc.out.list,
  thinningsteps = thinningsteps
)

# Burned + Thinned
mcmc.out.list.burn.thin <- postTHINN(
  mcmc.out.list = mcmc.out.list.burn,
  thinningsteps = thinningsteps
)

# ----------------------------------- -
# Convert scores to mcmc.list objects  ####
# ----------------------------------- -

# Original scores
list.mc.out.score <- mcmc.list(
  lapply(modelMCMC$mcmc.out.list, function(x) mcmc(x$scores))
)

# Thinned only
list.mc.out.thin.score <- mcmc.list(
  lapply(mcmc.out.list.thin, function(x) mcmc(x$scores))
)

# Burned + thinned
list.mc.out.burn.thin.score <- mcmc.list(
  lapply(mcmc.out.list.burn.thin, function(x) mcmc(x$scores))
)

# ----------------------------------- -
# Combine DAG arrays for burned + thinned  ####
# ----------------------------------- -
list.mc.out.burn.thin.dag <- abind::abind(
  lapply(mcmc.out.list.burn.thin, `[[`, "dags"),
  along = 3
)

# ----------------------------------- -
# Record postprocessing end time ####
# ----------------------------------- -
run_meta$postprocessing_end <- Sys.time()
run_meta <- log_step(run_meta, "Postprocessing MCMC finished")

# calculate postprocessing runtime in hours
run_meta$postprocessing_hours <- as.numeric(
  difftime(run_meta$postprocessing_end, run_meta$postprocessing_start, units = "hours")
)

run_meta <- log_step(run_meta, "MCMC plotting and QA start")
run_meta$plotting_start <- Sys.time()

#   ============================ -
#   MCMC Plots ####
#   ============================ -

# ----------------------------------- -
## Quality Assessment ####
# ----------------------------------- -
MCMCQualityAssessment(
  mcmcoutput = list.mc.out.score,
  mcmcoutputburnedthinned = list.mc.out.burn.thin.score,
  FILENAMEbase = paste0(validation_dir, "/"),
  filename = "GelmanPlot"
)
run_meta <- log_step(run_meta, "Quality Assessment: GelmanPlot completed")

# ----------------------------------- -
## Trace Plots (thinned + burned) ####
# ----------------------------------- -
TracePlot(
  mcouthin1 = mcmc.out.list.burn.thin[[1]],
  mcoutthin2 = mcmc.out.list.burn.thin[[2]],
  mcoutthin3 = mcmc.out.list.burn.thin[[3]],
  mcoutthin4 = mcmc.out.list.burn.thin[[4]],
  MCMC.SEEDS = MCMC.SEEDS,
  dat = datBN_red,
  dists = dists,
  METHOD = "mle",
  modeloutput = modelMCMC,
  filename = paste0("TracePlot_thinned", thinningsteps, "_burned", burnin.length),
  FILENAMEbase = paste0(validation_dir, "/")
)
run_meta <- log_step(run_meta, "TracePlot completed")

# ----------------------------------- -
# Consensus DAG  ####
# ----------------------------------- -
dag <- apply(list.mc.out.burn.thin.dag, 1:2, mean)

plotDAGthresholds(
  dag = dag,
  dists = dists,
  dat = datBN,
  threshold = c(0.5, 0.6),
  FILENAMEbase = paste0(validation_dir, "/"),
  filename = "DAG_different_thresholds"
)
run_meta <- log_step(run_meta, "Consensus DAG plots completed")

# ----------------------------------- -
## Arc strength threshold  ####
# ----------------------------------- -
arc.stren.sign.threshold <- arc.stren.threshold(dag, method = "l1")

ecdf <- ecdf_plot_ggplot(
  strengthdag = dag,
  arc.stren.sign.threshold = arc.stren.sign.threshold,
  filename = "ECDF_arc_sten_threshold2",
  FILENAMEbase = paste0(validation_dir, "/")
)
run_meta <- log_step(run_meta, "ECDF plot completed")

# DAG at significance threshold
plotDAGthresholds(
  threshold = arc.stren.sign.threshold,
  dag = dag,
  dists = dists,
  dat = datBN_red,
  FILENAMEbase = paste0(validation_dir, "/"),
  filename = "DAG_arcstreng_threshold"
)
run_meta <- log_step(run_meta, "DAG plot at arc strength threshold completed")

# Heatmap and adjacency matrix
ThresholdHeatmap(
  dag = dag,
  banmat = banmat,
  threshold = arc.stren.sign.threshold,
  FILENAMEbase = validation_dir,
  filename = "ThresholdHeatmap"  
)

ThresholdHeatmapSingle(
  dag = dag,
  banmat = banmat,
  threshold = arc.stren.sign.threshold,
  FILENAMEbase = validation_dir,
  filename = "ThresholdHeatmapSingle")


run_meta <- log_step(run_meta, "ThresholdHeatmaps completed")

AdjacencyMatrix(
  banmat = banmat,
  dag = dag,
  threshold = arc.stren.sign.threshold,
  filename = "AdjacencyMatrix",  
  FILENAMEbase = validation_dir
)
run_meta <- log_step(run_meta, "AdjacencyMatrix completed")

# DAG comparison
abnMCMCabnComparision(
  threshold = arc.stren.sign.threshold ,
  dag = dag,
  dists = dists,
  dat = datBN,
  prior = 1,
  banmat = banmat,
  retainmat = retainmat,
  max.parents = modelmaxparents$max.par,
  filename = "DAG_Comparison" ,
  FILENAMEbase = paste0(validation_dir, "/"),
  group.var   = "Source"
)
run_meta <- log_step(run_meta, "DAG comparison completed")

## ---------------------------- -
## MCMC Consensus DAG -> bnlearn DAG + CPDAG  ####
## ---------------------------- -

## 1. Threshold posterior arc strengths
threshold <- arc.stren.sign.threshold 

dag_threshold <- dag

dag_threshold[dag_threshold > threshold] <- 1
dag_threshold[dag_threshold <= threshold] <- 0

# Remove self-loops
diag(dag_threshold) <- 0

run_meta <- log_step(
  run_meta,
  paste0(
    "Thresholded MCMC consensus DAG at arc strength = ",
    threshold
  )
)


## 2. Ensure node names are available
nodes <- colnames(dag_threshold)

# Fallback if colnames are missing
if (is.null(nodes)) {
  
  nodes <- colnames(datBN_red)
  
  colnames(dag_threshold) <- nodes
  rownames(dag_threshold) <- nodes
}

nodes <- as.character(nodes)

# Make sure both dimensions have names
rownames(dag_threshold) <- nodes
colnames(dag_threshold) <- nodes



## 3. Convert thresholded consensus DAG to bnlearn object

dagMCMC_bn <- bnlearn::empty.graph(nodes)

# ABN -> bnlearn orientation
bnlearn::amat(dagMCMC_bn) <- t(dag_threshold)

run_meta <- log_step(
  run_meta,
  "Converted thresholded MCMC consensus DAG to bnlearn object"
)



## 4. Check whether graph is acyclic
if (!bnlearn::acyclic(dagMCMC_bn)) {
  
  stop(
    "The thresholded MCMC consensus graph contains cycles. ",
    "A CPDAG cannot be computed from this graph."
  )
}

run_meta <- log_step(
  run_meta,
  "Verified that thresholded MCMC consensus DAG is acyclic"
)



## 5. Construct CPDAG
cons.cpdag <- bnlearn::cpdag(dagMCMC_bn)

cons.amat <- bnlearn::amat(cons.cpdag)

run_meta <- log_step(
  run_meta,
  "Constructed CPDAG from thresholded MCMC consensus DAG"
)



## 6. Apply banmatrix constraints
# Ensure banmatrix has node names
if (is.null(rownames(banmat))) {
  rownames(banmat) <- colnames(banmat) <- seq_len(ncol(banmat))
}

# Match banmatrix to consensus DAG node order
banmat_sub <- banmat[nodes, nodes, drop = FALSE]

# Apply forbidden-direction constraints
#
# The transpose is retained from the previous implementation:
# banmat[i, j] = 1 corresponds to a forbidden j -> i arc
#
# Thus:
# cons.amat[i, j] == 1
# AND
# banmat_sub[j, i] == 1
# => remove i -> j

cons.amat[
  cons.amat == 1 &
    t(banmat_sub) == 1
] <- 0



## 7. Convert constrained CPDAG to bnlearn object

cons.mpdag <- bnlearn::empty.graph(nodes)

bnlearn::amat(cons.mpdag) <- cons.amat

run_meta <- log_step(
  run_meta,
  "Applied banmatrix constraints to consensus CPDAG"
)



## 8. Plot DAG and constrained CPDAG / MPDAG
pdf(
  file = file.path(
    validation_dir,
    paste0(run_name, "_cpdag_MCMC.pdf")
  ),
  width = 14,
  height = 7
)

par(mfrow = c(1, 2))

# Thresholded MCMC consensus DAG
plot(
  dagMCMC_bn,
  main = "MCMC consensus DAG",
  cex = 0.8
)

# Consensus CPDAG constrained by banmatrix
plot(
  cons.mpdag,
  main = "MCMC consensus MPDAG",
  cex = 0.8)

par(mfrow = c(1, 1))

dev.off()

run_meta <- log_step(
  run_meta,
  "Consensus DAG and banmatrix-constrained MPDAG plotted"
)


## Plot: Consenus DAG with edge strengths
# PDF
pdf(
  file = file.path(
    validation_dir,
    paste0(run_name, "_dagMP_MCMC.pdf")
  ),
  width = 8,
  height = 8
)

dag_threshold <- dag
dag_threshold[dag_threshold > arc.stren.sign.threshold ] <- 1
dag_threshold[dag_threshold <= arc.stren.sign.threshold ] <- 0

edgestren <- dag
edgestren[edgestren <= arc.stren.sign.threshold ] <- 0

cons.dag.plt.edgestrength <- plotAbn(
  dag = dag_threshold,
  data.df = dat,
  data.dists = dists,
  digits = 2,
  edge.strength = edgestren,
  plot = TRUE,
  main = arc.stren.sign.threshold
)

dev.off()


# SVG
svg(
  file = file.path(
    validation_dir,
    paste0(run_name, "_dagMP_MCMC.svg")
  ),
  width = 8,
  height = 8
)

cons.dag.plt.edgestrength <- plotAbn(
  dag = dag_threshold,
  data.df = dat,
  data.dists = dists,
  digits = 2,
  edge.strength = edgestren,
  plot = TRUE,
  main = arc.stren.sign.threshold
)

dev.off()

run_meta <- log_step(
  run_meta,
  "Most probable DAG of MCMC plotted"
)


## ============================ -
## Finalize and save metadata  ####
## ============================ -

run_meta$end_time <- Sys.time()

runtime_secs <- as.numeric(difftime(
  run_meta$end_time,
  run_meta$start_time,
  units = "secs"
))

run_meta$runtime_hours <- runtime_secs / 3600

meta_file <- file.path(out_dir, paste0(run_name, "_meta.txt"))

sink(meta_file)

cat("=====================================\n")
cat("Bayesian Network Run Metadata\n")
cat("=====================================\n\n")

cat("Run name:\n  ", run_meta$run_name, "\n\n")
cat("Output directory:\n  ", out_dir, "\n\n")

cat("Start time:\n  ", format(run_meta$start_time), "\n\n")
cat("End time:\n  ", format(run_meta$end_time), "\n\n")

cat("Runtime (hours):\n  ",
    round(run_meta$runtime_hours, 3), "\n\n")

cat("--------------------------------------------------\n")
cat("Steps executed:\n")
cat("--------------------------------------------------\n")

for (s in run_meta$steps) {
  cat(" - ", s, "\n")
}

cat("\n--------------------------------------------------\n")
cat("Settings:\n")
cat("--------------------------------------------------\n")

print_setting <- function(x, indent = 2) {
  pad <- paste(rep(" ", indent), collapse = "")
  for (nm in names(x)) {
    if (is.list(x[[nm]])) {
      cat(pad, nm, ":\n", sep = "")
      print_setting(x[[nm]], indent + 2)
    } else {
      cat(pad, nm, ": ", x[[nm]], "\n", sep = "")
    }
  }
}

print_setting(run_meta$settings)


cat("\n--------------------------------------------------\n")
cat("Data information:\n")
cat("--------------------------------------------------\n")

cat(" n_obs  : ", run_meta$data_info$n_obs, "\n", sep = "")
cat(" n_vars : ", run_meta$data_info$n_vars, "\n\n", sep = "")

cat(" Variables:\n\n")
for (v in run_meta$data_info$vars) {
  cat("", v, ", ", sep = "/ ")
}

cat("\n--------------------------------------------------\n")
cat("Validation (MCMC):\n")
cat("--------------------------------------------------\n")

cat(" Output directory:\n  ", run_meta$validation$output_dir, "\n\n")

cat(" MCMC seeds:\n  ",
    paste(run_meta$settings$validation$mcmc_seeds, collapse = ", "),
    "\n\n", sep = ", ")

cat(" MCMC scheme:\n  ",
    paste(run_meta$settings$validation$mcmc_iterations, collapse = ", "),
    "\n\n", sep = "")

cat(" Burn-in length:\n  ",
    run_meta$settings$validation$burnin_length, "\n\n", sep = "")

cat(" Thinning steps:\n  ",
    run_meta$settings$validation$thinning_steps, "\n", sep = "")


cat("\n--------------------------------------------------\n")
cat("R session information:\n")
cat("--------------------------------------------------\n")
print(sessionInfo())

sink()
