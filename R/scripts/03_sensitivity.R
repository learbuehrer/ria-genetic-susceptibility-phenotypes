###############################################################################-
#### File Content ####
# Script: 03_sensitivity.R
#
# Purpose:
#   This script evaluates the robustness of the primary findings using a
#   series of sensitivity analyses.
#
# Sensitivity Analyses:
#   1. Finnish vs. non-Finnish cohorts
#   2. Cohort-specific effect estimates
#   3. Random-slope mixed-effects models
#
# Genetic Predictors:
#   - Polygenic Risk Score (PRS)
#   - Allele Count (AC)
#
# Outcomes:
#   - IA size at rupture
#   - Age at rupture
#   - IA multiplicity
#   - IA location
#
# Notes:
#   - Source is treated as a clustering variable.
#   - Variables ending in "_Z" represent z-standardised variables.
#   - All analyses are performed using the synthetic dataset provided
#     in this repository.
#
###############################################################################-

###############################################################################-
# PART I: FINNISH VS NON-FINNISH COHORTS ####
# Finnish cohorts:
#   - Helsinki
#   - Kuopio
#
# Non-Finnish cohorts:
#   - All remaining contributing sources
#
# This analysis evaluates whether the primary associations are
# driven by the Finnish genetic background.
###############################################################################-

### Descriptives ####
dat_s$fins <- ifelse(dat_s$Source %in% c("Helsinki", "Kuopio"), "Fin", "Non-Fin")

dat_s_fins <- dat_s %>% 
  filter(fins == "Fin") %>% 
  droplevels()

dat_s_nonfins <- dat_s %>% 
  filter(!fins == "Fin") %>% 
  droplevels()


###### Table 1 ####
vars = c("Sex",
         "Age.at.time.of.SAH","Hypertension",  
         "Smoking_Current_Former_No", "Multiple.IAs", 
         "IAlocation_cow", 
         "IAlocation_binary",
         "IA_size.at.time.of.SAH",
         "allele_count_Z", 
         "PC1_AVG_Z", "PC2_AVG_Z", "PC3_AVG_Z",
         "PC4_AVG_Z", "PC5_AVG_Z")

vars_nonnorm = c("Age.at.time.of.SAH", 
                 "IA_size.at.time.of.SAH",
                 "allele_count_Z", 
                 "PC1_AVG_Z", "PC2_AVG_Z", "PC3_AVG_Z",
                 "PC4_AVG_Z")

invisible(tab <- tab_df_html(vars = vars,
                             vars_cat = NULL ,
                             vars_nonnorm = vars_nonnorm,
                             dat=dat_s, strat = "fins",
                             caption="Descriptive statistics at baseline,
       stratified by fins.",
                             Levels = T,
                             Missing = F))
tab

###############################################################################-
# Principal Component Correlations
#
# Principal component correlations are examined after excluding
# Finnish cohorts because the Finnish population represents a well-known
# genetic isolate within Europe.
###############################################################################-
## Look at correlations between scores:
# Select only the relevant columns
cluster_data <- dat_s_nonfins[,c("PC1_AVG_Z", "PC2_AVG_Z", "PC3_AVG_Z", 
                                 "PC4_AVG_Z", "PC5_AVG_Z")]

cor_matrix <- cor(cluster_data, use = "pairwise.complete.obs", method = "pearson")


colnames(cor_matrix) <- c("PC1","PC2", "PC3", "PC4", "PC5")
rownames(cor_matrix) <- c("PC1","PC2", "PC3", "PC4", "PC5")


## publication:
p <- ggcorrplot(
  cor_matrix,
  type = "lower",
  lab = TRUE,
  lab_size = 3,
  colors = c("#005A8C", "white", "#A48700"),
  outline.color = "white"
) +
  theme_classic() +
  xlab("") +
  ylab("") +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")

p

###############################################################################-
#### Regression Models ####
# Models are fit separately in:
#
#   - Complete cohort
#   - Finnish cohorts
#   - Non-Finnish cohorts
#
# Forest plots compare effect estimates across these populations.
###############################################################################-
#### Original PRS
exp_variables <- c("PRStotal_Z", "Sex")

#### define the cohorts:
cohort_data <- list(
  "Complete cohort" = dat_s,
  "Fins" = dat_s_fins,
  "Non-Fins" = dat_s_nonfins)


##### For IA Size  ####
size_prs <- run_cohort_comparison(
  cohort_data = cohort_data,
  outcome = "IA_size.at.time.of.SAH.log_Z",
  model_type = "lmer",
  clinical = "Multiple.IAs",
  term1 = "decreased IA size",
  term2 = "increased IA size")

size_prs$figure


##### For Age at Rupture Size ####
age_prs <- run_cohort_comparison(
  cohort_data = cohort_data,
  outcome = "Age.at.time.of.SAH_Z",
  model_type = "lmer",
  clinical = "IAlocation_cow",
  term1 = "decreased age",
  term2 = "increased age")

age_prs$figure

##### For IA Multiplicity ####
mult_prs <- run_cohort_comparison(
  cohort_data = cohort_data,
  outcome = "Multiple.IAs",
  model_type = "glmer",
  clinical = "IAlocation_cow",
  effect_label = "OR",
  null = 1,
  term1 = "favors solitary IAs",
  term2 = "favors multiple IAs")

mult_prs$figure


##### For IA Location  ####
loc_prs <- run_cohort_comparison(
  cohort_data = cohort_data,
  outcome = "IAlocation_binary",
  model_type = "glmer",
  clinical = "Multiple.IAs",
  effect_label = "OR",
  null = 1,
  term1 = "favors Low-Medium",
  term2 = "favors High")

loc_prs$figure


###############################################################################-
#### Regression: Combined  ####
###############################################################################-
# Forest Plot Preparation
#
# Estimates from the complete cohort, Finnish cohorts and
# non-Finnish cohorts are combined to visually assess consistency
# of effect sizes.
###############################################################################-
exp_variables <- c("PRStotal_Z","Sex")

##### Size ####
size_results_prs_nonfins  <- run_pipeline(
  data = dat_s_nonfins,
  outcome = "IA_size.at.time.of.SAH.log_Z",
  model_type = "lmer",
  random = "1|Source",
  clinical = "Multiple.IAs")

size_results_prs_fins <- run_pipeline(
  data = dat_s_fins,
  outcome = "IA_size.at.time.of.SAH.log_Z",
  model_type = "lmer",
  random = "1|Source",
  clinical = "Multiple.IAs")

size_sensitivity <- SensitivityPRSForestPlot(
  main_models = size_results_prs,
  nonfins_models = size_results_prs_nonfins,
  fins_models = size_results_prs_fins,
  prs_term = "PRStotal_Z",
  clinical_label = "IA multiplicity",
  title = "",
  xlab = "β",
  x_limits = c(-0.3, 0.3),
  term1 = "decreased IA size",
  term2 = "increased IA size")

size_sensitivity$plot

##### Age #### 
age_results_prs_nonfins  <- run_pipeline(
  data = dat_s_nonfins,
  outcome = "Age.at.time.of.SAH_Z",
  model_type = "lmer",
  random = "1|Source",
  clinical = "IAlocation_cow")

age_results_prs_fins  <- run_pipeline(
  data = dat_s_fins,
  outcome = "Age.at.time.of.SAH_Z",
  model_type = "lmer",
  random = "1|Source",
  clinical = "IAlocation_cow")


age_sensitivity <- SensitivityPRSForestPlot(
  main_models = age_results_prs,
  nonfins_models = age_results_prs_nonfins,
  fins_models = age_results_prs_fins,
  prs_term = "PRStotal_Z",
  clinical_label = "IA Location",
  title = "",
  xlab = "β",
  x_limits = c(-0.15, 0.25),
  term1 = "decreased age",
  term2 = "increased age")

age_sensitivity$plot

##### Multiplicity ####
multiple_results_prs_nonfins  <- run_pipeline(
  data = dat_s_nonfins,
  outcome = "Multiple.IAs",
  model_type = "glmer",
  random = "1|Source",
  clinical = "IAlocation_cow")

multiple_results_prs_fins  <- run_pipeline(
  data = dat_s_fins,
  outcome = "Multiple.IAs",
  model_type = "glmer",
  random = "1|Source",
  clinical = "IAlocation_cow")


multiplicity_sensitivity <- SensitivityPRSForestPlot(
  main_models = multiple_results_prs,
  nonfins_models = multiple_results_prs_nonfins,
  fins_models = multiple_results_prs_fins,
  prs_term = "PRStotal_Z",
  clinical_label = "IA Location",
  title = "",
  xlab = "OR",
  x_limits = c(0.5,1.5),
  term1 = "favors solitary IAs",
  term2 = "favors multiple IAs")

multiplicity_sensitivity$plot


#####   Location ####
location_results_prs_nonfins  <- run_pipeline(
  data = dat_s_nonfins,
  outcome = "IAlocation_binary",
  model_type = "glmer",
  random = "1|Source",
  clinical = "Multiple.IAs")

location_results_prs_fins  <- run_pipeline(
  data = dat_s_fins,
  outcome = "IAlocation_binary",
  model_type = "glmer",
  random = "1|Source",
  clinical = "Multiple.IAs")


location_sensitivity <- SensitivityPRSForestPlot(
  main_models = location_results_prs,
  nonfins_models = location_results_prs_nonfins,
  fins_models = location_results_prs_fins,
  prs_term = "PRStotal_Z",
  clinical_label = "IA multiplicity",
  title = "",
  xlab = "OR",
  x_limits = c(0.5, 1.2),
  term1 = "Low-Medium location",
  term2 = "High location")

location_sensitivity$plot


###############################################################################-
# PART II: RANDOM-SLOPE MODELS ####
#
# Primary analyses assume that the PRS effect is identical across
# sources while allowing source-specific intercepts.
#
# Sensitivity analyses additionally allow the PRS effect to vary
# between sources:
#
#     (1 + PRStotal_Z | Source)
#
# This evaluates potential heterogeneity in genetic effects across
# contributing cohorts.
#
#
# Note:
# Random-slope models frequently approached singular fits due to the
# limited number of contributing cohorts and should therefore be
# interpreted as exploratory sensitivity analyses.
# Here only the results of the main models are shown, all other models
# can be printed identically using $model.
###############################################################################-

exp_variables <- c("Sex", "PRStotal_Z")

##### For IA Size ####
size_results_me <-
  run_pipeline(
    data=dat_s,
    outcome="IA_size.at.time.of.SAH.log_Z",
    model_type="lmer",
    random = "1|Source",
    clinical = "Multiple.IAs")


size_results_me_SA <-
  run_pipeline(
    data=dat_s,
    outcome="IA_size.at.time.of.SAH.log_Z",
    model_type="lmer",
    random = "1 + PRStotal_Z|Source",
    clinical = "Multiple.IAs")


## PCS
a <- ForestPlot(size_results_me$pcs, term1 = "decreased IA size",
                term2 = "increased IA size", title = "~(1|Source)")$plot
b <- ForestPlot(size_results_me_SA$pcs, term1 = "decreased IA size",
                term2 = "increased IA size", title = "~(1+PRStotal_Z|Source)")$plot
ggarrange(a,b, ncol=2)


random_slopes <- ranef(size_results_me_SA$pcs)$Source %>%
  tibble::rownames_to_column("Source") %>%
  dplyr::select(Source, PRStotal_Z) %>%
  rename(random_slope = PRStotal_Z)

ggplot(random_slopes,
       aes(x = reorder(Source, random_slope),
           y = random_slope)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    x = "Source",
    y = "Random slope deviation",
    title = "Source-specific random slopes for PRS"
  ) +
  theme_minimal()


##### For Age at Rupture Size ####
age_results_me <-
  run_pipeline(
    data=dat_s,
    outcome="Age.at.time.of.SAH_Z",
    model_type="lmer",
    random="1|Source",
    clinical = "IAlocation_cow")

age_results_me_SA <-
  run_pipeline(
    data=dat_s,
    outcome="Age.at.time.of.SAH_Z",
    model_type="lmer",
    random="1 + PRStotal_Z|Source",
    clinical = "IAlocation_cow")



## PCS
a <- ForestPlot(age_results_me$pcs, term1 = "decreased IA age",
                term2 = "increased IA age", title = "~(1|Source)")$plot
b <- ForestPlot(age_results_me_SA$pcs, term1 = "decreased IA age",
                term2 = "increased IA age", title = "~(1+PRStotal_Z|Source)")$plot
ggarrange(a,b, ncol=2)



### random slope:
random_slopes <- ranef(age_results_me_SA$pcs)$Source %>%
  tibble::rownames_to_column("Source") %>%
  dplyr::select(Source, PRStotal_Z) %>%
  rename(random_slope = PRStotal_Z)

ggplot(random_slopes,
       aes(x = reorder(Source, random_slope),
           y = random_slope)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    x = "Source",
    y = "Random slope deviation",
    title = "Source-specific random slopes for PRS"
  ) +
  theme_minimal()


##### For IA Multiplicity #####
multiple_results_me <-
  run_pipeline(
    data=dat_s,
    outcome="Multiple.IAs",
    model_type="glmer",
    random="1|Source",
    clinical = "IAlocation_cow")

multiple_results_me_SA <-
  run_pipeline(
    data=dat_s,
    outcome="Multiple.IAs",
    model_type="glmer",
    random="1+PRStotal_Z|Source",
    clinical = "IAlocation_cow")


## PCS
a <- ForestPlot(multiple_results_me$pcs, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs", title = "~(1|Source)")$plot
b <- ForestPlot(multiple_results_me_SA$pcs, term1 = "favors solitary IAs",
                term2 = "favors multiple IAs", title = "~(1+PRStotal_Z|Source)")$plot
ggarrange(a,b, ncol=2)



### random slope:
random_slopes <- ranef(multiple_results_me_SA$pcs)$Source %>%
  tibble::rownames_to_column("Source") %>%
  select(Source, PRStotal_Z) %>%
  rename(random_slope = PRStotal_Z)

ggplot(random_slopes,
       aes(x = reorder(Source, random_slope),
           y = random_slope)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    x = "Source",
    y = "Random slope deviation",
    title = "Source-specific random slopes for PRS"
  ) +
  theme_minimal()


##### For IA Location ####
location_results_me <-
  run_pipeline(
    data=dat_s,
    outcome="IAlocation_binary",
    model_type="glmer",
    random="1|Source",
    clinical = "Multiple.IAs")

location_results_me_SA <-
  run_pipeline(
    data=dat_s,
    outcome="IAlocation_binary",
    model_type="glmer",
    random="1+PRStotal_Z|Source",
    clinical = "Multiple.IAs")


## PCS
a <- ForestPlot(location_results_me$pcs, term1 = "Low-Medium locations",
                term2 = "High locations", title = "~(1|Source)")$plot
b <- ForestPlot(location_results_me_SA$pcs, term1 = "Low-Medium locations",
                term2 = "High locations", title = "~(1+PRStotal_Z|Source)")$plot
ggarrange(a,b, ncol=2)




### random slope:
random_slopes <- ranef(location_results_me_SA$pcs)$Source %>%
  tibble::rownames_to_column("Source") %>%
  dplyr::select(Source, PRStotal_Z) %>%
  rename(random_slope = PRStotal_Z)

ggplot(random_slopes,
       aes(x = reorder(Source, random_slope),
           y = random_slope)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    x = "Source",
    y = "Random slope deviation",
    title = "Source-specific random slopes for PRS"
  ) +
  theme_minimal()


###############################################################################-
# PART III: SOURCE-SPECIFIC ANALYSES ####
#
# To assess heterogeneity between contributing cohorts,
# separate regression models are fitted within each source.
#
# These analyses are exploratory and intended to identify
# potential cohort-specific patterns.
###############################################################################-
exp_variables <- c("PRStotal_Z", "Sex")


## remove Barcelona from the data set due to a too small sample size
dat_sourceplots <- dat_s %>%
  filter(!Source == "Barcelona") %>%
  droplevels()


##### IA Size ####
size_source <- plot_source_specific(
  data = dat_sourceplots,
  outcome = "IA_size.at.time.of.SAH.log_Z",
  model_type = "lm",
  clinical = "Multiple.IAs",
  title = "Sensitivity analysis: source-specific models (IA size)")

size_source$plot

##### Age at Rupture ####
age_source <- plot_source_specific(
  data = dat_sourceplots,
  outcome = "Age.at.time.of.SAH_Z",
  model_type = "lm",
  clinical = "IAlocation_cow",
  title = "Sensitivity analysis: source-specific models (Age at Rupture)")

age_source$plot


##### IA Multiplicity ####
multi_source <- plot_source_specific(
  data = dat_sourceplots,
  outcome = "Multiple.IAs",
  model_type = "glm",
  clinical = "IAlocation_binary",
  exponentiate = TRUE,
  measure = "OR",
  null_value = 1,
  title = "Sensitivity analysis: source-specific models (IA multiplicity)")

multi_source$plot

# Note:
# Source-specific logistic models showed convergence issues for GERFHS
# and Oxford in the synthetic data set, likely due to limited event counts within these cohorts.

##### IA location ####
location_source <- plot_source_specific(
  data = dat_sourceplots,
  outcome = "IAlocation_binary",
  model_type = "glm",
  clinical = "Multiple.IAs",
  exponentiate = TRUE,
  measure = "OR",
  null_value = 1,
  title = "Sensitivity analysis: source-specific models (IA location)")

location_source$plot
# Note:
# Source-specific logistic models showed convergence issues for GERFHS in the
# synthetic data set, likely due to limited event counts within these cohorts.