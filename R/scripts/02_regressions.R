############################################################################### -
#### File Content ####
# Script: 02_regressions.R
#
# Purpose:
#   This script evaluates associations between genetic risk measures and
#   aneurysm characteristics using mixed-effects regression models.
#
# Genetic Risk Measures:
#   - Polygenic Risk Score (PRS)
#   - Allele Count (AC)
#
# Outcomes:
#   - IA size at rupture (log, scaled)
#   - Age at rupture (scaled)
#   - IA multiplicity (binary)
#   - IA location (binary)
#
# Statistical Models:
#   - Linear mixed-effects models (lmer) for continuous outcomes
#   - Generalized linear mixed-effects models (glmer) for binary outcomes
#
# Random Effects:
#   - Source included as a random intercept to account for clustering
#     across contributing cohorts.
#
# Multiple Testing:
#   - Holm-adjusted p-values are calculated separately for PRS and AC
#     analyses across the four primary outcomes.
#
# Output:
#   - Forest plots for all model formulations
#   - Model diagnostics
#   - Combined publication figures
#   - Supplementary figures
#
# Notes:
#   - Source is treated as a clustering variable.
#   - Variables ending in "_Z" represent z-standardised variables.
#   - All analyses are performed using the synthetic dataset provided
#     in this repository.

###############################################################################-

###############################################################################-
#### PART I: POLYGENIC RISK SCORE (PRS) ANALYSES ####
# Forest Plots
#
# The pipeline fits six nested model formulations:
#
#   1. Base model
#   2. Base + ancestry PCs
#   3. Clinical model
#   4. Family history model
#   5. Risk-factor model
#   6. Personal-history model
#
# Model Diagnostics
#
# Diagnostic plots are generated to assess:
#
#   - Residual behaviour
#   - Homoscedasticity
#   - Model fit
#   - Deviations from normality
#
# Continuous outcomes:
#   - Tukey-Anscombe plots
#   - QQ-plots
#
###############################################################################-
## Exposure Variables
exp_variables <- c("PRStotal_Z", "Sex") # Main PRS analysis

##### For IA Size  ####
size_results_prs  <- run_pipeline(
  data = dat_s,
  outcome = "IA_size.at.time.of.SAH.log_Z",
  model_type = "lmer",
  random = "1|Source",
  clinical = "Multiple.IAs")


a <- ForestPlot(size_results_prs$base, term1 = "decreased IA size",
                term2 = "increased IA size")$plot
b <- ForestPlot(size_results_prs$pcs, term1 = "decreased IA size",
                term2 = "increased IA size")$plot
c <- ForestPlot(size_results_prs$clinical, term1 = "decreased IA size",
                term2 = "increased IA size")$plot
d <- ForestPlot(size_results_prs$famhist, term1 = "decreased IA size",
                term2 = "increased IA size")$plot
e <- ForestPlot(size_results_prs$riskfactors, term1 = "decreased IA size",
                term2 = "increased IA size")$plot
f <- ForestPlot(size_results_prs$personal, term1 = "decreased IA size",
                term2 = "increased IA size")$plot

## have at a look at the individual models by just calling the model object, e.g. b,f
b
f

# ## saving
# ggsave(
#   "figures/Appendix_size_prs_base.pdf",
#   a,
#   width = 150,
#   height = 80,
#   units = "mm",
#   device = cairo_pdf)
# 
# ggsave(
#   "figures/Appendix_size_prs_pcs.pdf",
#   b,
#   width = 150,
#   height = 80,
#   units = "mm",
#   device = cairo_pdf)
# 
# ggsave(
#   "figures/Appendix_size_prs_clinical.pdf",
#   c,
#   width = 150,
#   height = 80,
#   units = "mm",
#   device = cairo_pdf)
# 
# ggsave(
#   "figures/Appendix_size_prs_famhist.pdf",
#   d,
#   width = 150,
#   height = 80,
#   units = "mm",
#   device = cairo_pdf)
# 
# ggsave(
#   "figures/Appendix_size_prs_riskfactors.pdf",
#   e,
#   width = 150,
#   height = 80,
#   units = "mm",
#   device = cairo_pdf)
# 
# ggsave(
#   "figures/Appendix_size_prs_personal.pdf",
#   f,
#   width = 150,
#   height = 80,
#   units = "mm",
#   device = cairo_pdf)

size_diagnostics <- diagnose_pipeline(size_results_prs)
size_diagnostics$table

TA_size <- TukeyAnscombe_pipeline(size_results_prs, ncol = 3)
TA_size

# ggsave(
#   "figures/Appendix_TA_size_prs.pdf",
#   TA_size,
#   width = 250,
#   height = 150,
#   units = "mm",
#   device = cairo_pdf)

QQ_size <- QQ_pipeline(
  size_results_prs,
  ncol = 3)

QQ_size

# ggsave(
#   "figures/Appendix_QQ_size_prs.pdf",
#   QQ_size,
#   width = 250,
#   height = 150,
#   units = "mm",
#   device = cairo_pdf)


##### For Age at Rupture Size  ####
age_results_prs <-
  run_pipeline(
    data=dat_s,
    outcome="Age.at.time.of.SAH_Z",
    model_type="lmer",
    random = "1|Source",
    clinical = "IAlocation_cow")

a <- ForestPlot(age_results_prs$base, term1 = "decreased age",
                term2 = "increased age")$plot
b <- ForestPlot(age_results_prs$pcs, term1 = "decreased age",
                term2 = "increased age")$plot
c <- ForestPlot(age_results_prs$clinical, term1 = "decreased age",
                term2 = "increased age")$plot
d <- ForestPlot(age_results_prs$famhist, term1 = "decreased age",
                term2 = "increased age")$plot
e <- ForestPlot(age_results_prs$riskfactors, term1 = "decreased age",
                term2 = "increased age")$plot
f <- ForestPlot(age_results_prs$personal, term1 = "decreased age",
                term2 = "increased age")$plot


# model diagnostics
age_diagnostics <- diagnose_pipeline(age_results_prs)
age_diagnostics$table

TA_age <- TukeyAnscombe_pipeline(age_results_prs,ncol = 3)
TA_age

QQ_age <- QQ_pipeline(
  age_results_prs,
  ncol = 3)

QQ_age


##### For IA Multiplicity  ####
multiple_results_prs <-
  run_pipeline(
    data=dat_s,
    outcome="Multiple.IAs",
    model_type="glmer",
    random = "1|Source",
    clinical = "IAlocation_cow")

a <- ForestPlot(multiple_results_prs$base, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot
b <- ForestPlot(multiple_results_prs$pcs, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot
c <- ForestPlot(multiple_results_prs$clinical, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot
d <- ForestPlot(multiple_results_prs$famhist, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot
e <- ForestPlot(multiple_results_prs$riskfactors, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot
f <- ForestPlot(multiple_results_prs$personal, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot


# model diagnostics
multiple_diagnostics <- diagnose_pipeline(multiple_results_prs)
multiple_diagnostics$table

TA_multiple <- TukeyAnscombe_pipeline(multiple_results_prs,ncol = 3)
TA_multiple


##### For IA Location  ####
location_results_prs <-
  run_pipeline(
    data=dat_s,
    outcome="IAlocation_binary",
    model_type="glmer",
    random = "1|Source",
    clinical = "Multiple.IAs")

a <- ForestPlot(location_results_prs$base, term1 = "Low-Medium locations",
                term2 = "High locations")$plot
b <- ForestPlot(location_results_prs$pcs, term1 = "Low-Medium locations",
                term2 = "High locations")$plot
c <- ForestPlot(location_results_prs$clinical, term1 = "Low-Medium locations",
                term2 = "High locations")$plot
d <- ForestPlot(location_results_prs$famhist, term1 = "Low-Medium locations",
                term2 = "High locations")$plot
e <- ForestPlot(location_results_prs$riskfactors, term1 = "Low-Medium locations",
                term2 = "High locations")$plot
f <- ForestPlot(location_results_prs$personal, term1 = "Low-Medium locations",
                term2 = "High locations")$plot

# model diagnostics
location_diagnostics <- diagnose_pipeline(location_results_prs)
location_diagnostics$table

TA_location <- TukeyAnscombe_pipeline(location_results_prs,ncol = 3)
TA_location


###############################################################################.-
#### PART II: ALLELE COUNT (AC) ANALYSES ####
###############################################################################-
exp_variables <- c("allele_count_Z","Sex")

##### For IA Size ####
size_results_ac  <- run_pipeline(
  data = dat_s,
  outcome = "IA_size.at.time.of.SAH.log_Z",
  model_type = "lmer",
  random = "1|Source",
  clinical = "Multiple.IAs")

a <- ForestPlot(size_results_ac $base, term1 = "decreased IA size",
                term2 = "increased IA size")$plot
b <- ForestPlot(size_results_ac $pcs, term1 = "decreased IA size",
                term2 = "increased IA size")$plot
c <- ForestPlot(size_results_ac $clinical, term1 = "decreased IA size",
                term2 = "increased IA size")$plot
d <- ForestPlot(size_results_ac $famhist, term1 = "decreased IA size",
                term2 = "increased IA size")$plot
e <- ForestPlot(size_results_ac $riskfactors, term1 = "decreased IA size",
                term2 = "increased IA size")$plot
f <- ForestPlot(size_results_ac $personal, term1 = "decreased IA size",
                term2 = "increased IA size")$plot


## model diagnostics
size_diagnostics <- diagnose_pipeline(size_results_ac)
size_diagnostics$table

TA_size <- TukeyAnscombe_pipeline(size_results_ac, ncol = 3)
TA_size


QQ_size <- QQ_pipeline(
  size_results_ac,
  ncol = 3)

QQ_size


##### For Age at Rupture Size ####
age_results_ac <-
  run_pipeline(
    data=dat_s,
    outcome="Age.at.time.of.SAH_Z",
    model_type="lmer",
    random = "1|Source",
    clinical = "IAlocation_cow")

a <- ForestPlot(age_results_ac$base, term1 = "decreased age",
                term2 = "increased age")$plot
b <- ForestPlot(age_results_ac$pcs, term1 = "decreased age",
                term2 = "increased age")$plot
c <- ForestPlot(age_results_ac$clinical, term1 = "decreased age",
                term2 = "increased age")$plot
d <- ForestPlot(age_results_ac$famhist, term1 = "decreased age",
                term2 = "increased age")$plot
e <- ForestPlot(age_results_ac$riskfactors, term1 = "decreased age",
                term2 = "increased age")$plot
f <- ForestPlot(age_results_ac$personal, term1 = "decreased age",
                term2 = "increased age")$plot


# model diagnostics
age_diagnostics <- diagnose_pipeline(age_results_ac)
age_diagnostics$table

TA_age <- TukeyAnscombe_pipeline(age_results_ac,ncol = 3)
TA_age

QQ_age <- QQ_pipeline(
  age_results_ac,
  ncol = 3)

QQ_age


##### For IA Multiplicity  ####
multiple_results_ac <-
  run_pipeline(
    data=dat_s,
    outcome="Multiple.IAs",
    model_type="glmer",
    random = "1|Source",
    clinical = "IAlocation_cow")

a <- ForestPlot(multiple_results_ac$base, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot
b <- ForestPlot(multiple_results_ac$pcs, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot
c <- ForestPlot(multiple_results_ac$clinical, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot
d <- ForestPlot(multiple_results_ac$famhist, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot
e <- ForestPlot(multiple_results_ac$riskfactors, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot
f <- ForestPlot(multiple_results_ac$personal, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs")$plot


# model diagnostics
multiple_diagnostics <- diagnose_pipeline(multiple_results_ac)
multiple_diagnostics$table

TA_multiple <- TukeyAnscombe_pipeline(multiple_results_ac,ncol = 3)
TA_multiple


##### For IA Location  ####
location_results_ac <-
  run_pipeline(
    data=dat_s,
    outcome="IAlocation_binary",
    model_type="glmer",
    random = "1|Source",
    clinical = "Multiple.IAs")

a <- ForestPlot(location_results_ac$base, term1 = "Low-Medium locations",
                term2 = "High locations")$plot
b <- ForestPlot(location_results_ac$pcs, term1 = "Low-Medium locations",
                term2 = "High locations")$plot
c <- ForestPlot(location_results_ac$clinical, term1 = "Low-Medium locations",
                term2 = "High locations")$plot
d <- ForestPlot(location_results_ac$famhist, term1 = "Low-Medium locations",
                term2 = "High locations")$plot
e <- ForestPlot(location_results_ac$riskfactors, term1 = "Low-Medium locations",
                term2 = "High locations")$plot
f <- ForestPlot(location_results_ac$personal, term1 = "Low-Medium locations",
                term2 = "High locations")$plot

# model diagnostics
location_diagnostics <- diagnose_pipeline(location_results_ac)
location_diagnostics$table

TA_location <- TukeyAnscombe_pipeline(location_results_ac,ncol = 3)
TA_location

###############################################################################-
#### PART III: COMBINED PUBLICATION FIGURES ####
# Combined Forest Plots
#
# Publication figures compare estimates obtained from:
#
#   - Polygenic risk scores (PRS)
#   - Allele counts (AC)
#
# Holm-adjusted p-values are displayed to account for multiple testing
# across the four primary outcomes.
###############################################################################-

## Adjust for multiple testing:
prs_holm <- create_holm_table(
  size_models = size_results_prs,
  age_models = age_results_prs,
  multiplicity_models = multiple_results_prs,
  location_models = location_results_prs,
  term = "PRStotal_Z")

ac_holm <- create_holm_table(
  size_models = size_results_ac,
  age_models = age_results_ac,
  multiplicity_models = multiple_results_ac,
  location_models = location_results_ac,
  term = "allele_count_Z")


## Plot for Size with adjusted p-values:
size_forest <- CombinedPRSForestPlot(
  prs_models = size_results_prs,
  ac_models = size_results_ac,
  p_adjust_table = prs_holm,
  outcome_name = "size",
  show = c("PRS"),
  term1 = "decreased IA size",
  term2 = "increased IA size")

size_forest$plot

ggsave(
  filename = "figures/02_IA_size_PRS.pdf",
  plot = size_forest$plot,
  width = 5,
  height = 3.75,
  units = "in",
  device = cairo_pdf,
  bg = "white")


size_forest_AC <- CombinedPRSForestPlot(
  prs_models = size_results_prs,
  ac_models = size_results_ac,
  p_adjust_table = ac_holm,
  outcome_name = "size",
  show = "AC",
  term1 = "decreased IA size",
  term2 = "increased IA size")

size_forest_AC$plot


ggsave(
  filename = "figures/Appendix_IA_size_AC.pdf",
  plot = size_forest_AC$plot,
  width = 5,
  height = 3.75,
  units = "in",
  device = cairo_pdf,
  bg = "white")


## Plot for Age with adjusted p-values:
age_forest <- CombinedPRSForestPlot(
  prs_models = age_results_prs,
  ac_models = age_results_ac,
  p_adjust_table = prs_holm,
  outcome_name = "age",
  show = c("PRS"),
  term1 = "decreased age",
  term2 = "increased age")

age_forest$plot

ggsave(
  filename = "figures/02_age_PRS.pdf",
  plot = age_forest$plot,
  width = 5,
  height = 3.75,
  units = "in",
  device = cairo_pdf,
  bg = "white")

age_forest_AC <- CombinedPRSForestPlot(
  prs_models = age_results_prs,
  ac_models = age_results_ac,
  p_adjust_table = ac_holm,
  outcome_name = "age",
  show = "AC",
  term1 = "decreased age",
  term2 = "increased age")

age_forest_AC$plot

ggsave(
  filename = "figures/Appendix_age_AC.pdf",
  plot = age_forest_AC$plot,
  width = 5,
  height = 3.75,
  units = "in",
  device = cairo_pdf,
  bg = "white")


## Plot for Multiplicity with adjusted p-values:
multiplicity_forest <- CombinedPRSForestPlot(
  prs_models = multiple_results_prs,
  ac_models = multiple_results_ac,
  show = c("PRS"),
  p_adjust_table = prs_holm,
  outcome_name = "multiplicity",
  term1 = "favors solitary IAs",
  term2 = "favors multiple IAs")

multiplicity_forest$plot


ggsave(
  filename = "figures/02_IA_multiplicity_PRS.pdf",
  plot = multiplicity_forest$plot,
  width = 5,
  height = 3.75,
  units = "in",
  device = cairo_pdf,
  bg = "white")


multiplicity_forest_AC <- CombinedPRSForestPlot(
  prs_models = multiple_results_prs,
  ac_models = multiple_results_ac,
  show = c("AC"),
  p_adjust_table = ac_holm,
  outcome_name = "multiplicity",
  term1 = "favors solitary IAs",
  term2 = "favors multiple IAs")

multiplicity_forest_AC$plot


ggsave(
  filename = "figures/Appendix_IA_multiplicity_AC.pdf",
  plot = multiplicity_forest_AC$plot,
  width = 5,
  height = 3.75,
  units = "in",
  device = cairo_pdf,
  bg = "white")


## Plot for Location with adjusted p-values:
location_forest <- CombinedPRSForestPlot(
  prs_models = location_results_prs,
  ac_models = location_results_ac,
  p_adjust_table = prs_holm,
  outcome_name = "location",
  show = c("PRS"),
  term1 = "Low-Medium",
  term2 = "High")

location_forest$plot

ggsave(
  filename = "figures/02_IA_location_PRS.pdf",
  plot = location_forest$plot,
  width = 5,
  height = 3.75,
  units = "in",
  device = cairo_pdf,
  bg = "white")


location_forest_AC <- CombinedPRSForestPlot(
  prs_models = location_results_prs,
  ac_models = location_results_ac,
  p_adjust_table = ac_holm,
  outcome_name = "location",
  show = c("AC"),
  term1 = "Low-Medium",
  term2 = "High")

location_forest_AC$plot

ggsave(
  filename = "figures/Appendix_IA_location_AC.pdf",
  plot = location_forest_AC$plot,
  width = 5,
  height = 3.75,
  units = "in",
  device = cairo_pdf,
  bg = "white")


###############################################################################-
#### PART IV: MODELS FROM THE ABN STRUCTURE ####
# Models Derived from the ABN Network Structure
#
# To facilitate interpretation of selected dependencies identified in the
# additive Bayesian network (ABN), mixed-effects regression models are
# fitted using the corresponding variable combinations.
###############################################################################-
##### PRS total ####
fit_RIA_multiple_me <- glmer(Multiple.IAs~Sex+
                               IA_size.at.time.of.SAH.log_Z+
                               Smoking_Current_Former_No+
                               PRStotal_Z+
                               (1|Source),
                             family = "binomial",
                             data=dat_s)

ForestPlot(fit_RIA_multiple_me, term1 = "favors solitary IAs",
           term2 = "favors multiple IAs", exponentiate = TRUE)$plot



fit_RIA_multiple_me_AC <- glmer(Multiple.IAs~Sex+
                                  IA_size.at.time.of.SAH.log_Z+
                                  Smoking_Current_Former_No+
                                  allele_count_Z+
                                  (1|Source),
                                family = "binomial",
                                data=dat_s)

ForestPlot(fit_RIA_multiple_me_AC, term1 = "favors solitary IAs",
           term2 = "favors multiple IAs", exponentiate = TRUE)$plot



##### Cluster 2 ####
fit_RIA_multiple_me_cluster <- glmer(Multiple.IAs~Sex+
                                       Smoking_Current_Former_No+
                                       PRS2_Z+
                                       (1|Source),
                                     family = "binomial",
                                     data=dat_s)

ForestPlot(fit_RIA_multiple_me_cluster, term1 = "favors solitary IAs",
           term2 = "favors multiple IAs", exponentiate = TRUE)$plot

