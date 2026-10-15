###############################################################################-
#### File Content ####
# Script: 01_descriptives.R
#
# Purpose:
#   This script generates the descriptive statistics and exploratory figures
#   used to characterise the study cohort.
#
# Input:
#   - Synthetic study dataset (dat_s)
#   - Project-specific helper functions sourced in the setup script
#
# Output:
#   Tables:
#     - Overall descriptive statistics
#     - Descriptive statistics stratified by source
#     - Computation of intraclass correlation coefficients (ICC) for selected variables
#
#   Figures:
#     - Cohort-specific baseline characteristics
#     - Distribution of polygenic risk scores (PRS)
#     - Distribution of allele counts (AC)
#     - Principal component (PC) visualisations
#     - Correlation analyses
#
# Notes:
#   - The repository contains synthetic data generated from the original
#     study cohort. These data are provided exclusively to reproduce the
#     analytical workflow and do not represent real study participants.
#
#   - Variables ending in "_Z" denote standardised variables
#     (z-transformed variables) on the synthetic dataset:
#
#         Z = (x - mean(x)) / sd(x)
#
#   - Make sure 00_setup.R has been sourced before running this script 
#     to ensure that the required packages/ data anf functions are loaded.
#
###############################################################################-


###############################################################################-
#### Tab 1: DESCRIPTIVE STATISTICS ####
###############################################################################-
vars = c("Sex",
         "Hypertension",  
         "Smoking_Current_Former_No", 
         "Positive.famillial.history",
         "Multiple.IAs",
         "Age.at.time.of.SAH",
         "IA_size.at.time.of.SAH",
         "IAlocation_cow")

vars_nonnorm = c(
  "Age.at.time.of.SAH","Multiple.IAs",
  "IA_size.at.time.of.SAH")


## Overall (in HTML format)
invisible(tab <- tab_df_html(vars = vars,
                             vars_cat = NULL ,
                             vars_nonnorm = vars_nonnorm,
                             dat=dat_s, overall = TRUE,
                             caption="Descriptive statistics at baseline.",
                             Levels = FALSE,
                             Missing = FALSE))

tab


## By Source (in Latex format)
invisible(tab <- tab_df(vars = vars,
                        vars_cat = NULL ,
                        vars_nonnorm = vars_nonnorm,
                        dat=dat_s, strat = "Source",
                        caption="Descriptive statistics at baseline,
       stratified by source.",
                        Levels = FALSE,
                        Missing = FALSE,
                        SMD=FALSE,
                        p=FALSE,
                        label = "tab:tab1PerSource",
                        filename="tables/TableOnePerSource.tex"))


###############################################################################-
#### Fig. 1a: Baseline characteristics by source ####
# Binary variables are displayed as cohort-specific percentages.
# Continuous variables are displayed as cohort-specific boxplots.
# Dashed lines indicate overall cohort means.
###############################################################################-

# Variable definitions
binary_labels <- c(
  Female           = "Female Sex",
  CurrentSmoking   = "Current Smoker",
  Hypertension     = "Hypertension Aware",
  FamilyHistory    = "Family History Aware",
  Multiplicity     = "Multiple IAs",
  HighRiskLocation = "High-risk IA Location")

continuous_labels <- c(
  AgeAtRupture = "Age at Rupture",
  IASize       = "IA Size at Rupture")

source_levels <- rev(levels(dat_s$Source))


# Binary variables
binary_summary <-
  dat_s %>%
  group_by(Source) %>%
  summarise(
    Female =
      mean(Sex == "Female", na.rm = TRUE) * 100,
    CurrentSmoking =
      mean(
        Smoking_Current_Former_No == "Current",
        na.rm = TRUE) * 100,
    Hypertension =
      mean(
        Hypertension == "AnyType",
        na.rm = TRUE) * 100,
    FamilyHistory =
      mean(
        Positive.famillial.history == "Yes",
        na.rm = TRUE
      ) * 100,
    Multiplicity =
      mean(
        Multiple.IAs == "Yes",
        na.rm = TRUE) * 100,
    HighRiskLocation =
      mean(
        IAlocation_binary == "High",
        na.rm = TRUE) * 100,
    .groups = "drop")

binary_plot_data <-
  binary_summary %>%
  pivot_longer(
    -Source,
    names_to = "Variable",
    values_to = "Value")

binary_overall <-
  dat_s %>%
  summarise(
    Female =
      mean(Sex == "Female", na.rm = TRUE) * 100,
    CurrentSmoking =
      mean(
        Smoking_Current_Former_No == "Current",
        na.rm = TRUE) * 100,
    Hypertension =
      mean(
        Hypertension == "AnyType",
        na.rm = TRUE) * 100,
    FamilyHistory =
      mean(
        Positive.famillial.history == "Yes",
        na.rm = TRUE) * 100,
    Multiplicity =
      mean(
        Multiple.IAs == "Yes",
        na.rm = TRUE) * 100,
    HighRiskLocation =
      mean(
        IAlocation_binary == "High",
        na.rm = TRUE) * 100
  ) %>%
  pivot_longer(
    everything(),
    names_to = "Variable",
    values_to = "Overall")


# Continuous variables
continuous_plot_data <-
  dat_s %>%
  select(
    Source,
    AgeAtRupture = Age.at.time.of.SAH,
    IASize = IA_size.at.time.of.SAH.log
  ) %>%
  pivot_longer(
    -Source,
    names_to = "Variable",
    values_to = "Value") %>%
  filter(is.finite(Value))

continuous_overall <-
  continuous_plot_data %>%
  group_by(Variable) %>%
  summarise(
    Overall = mean(Value, na.rm = TRUE),
    .groups = "drop")



# Labels and factor ordering
label_variable <- function(x, labels) {
  factor(
    x,
    levels = names(labels),
    labels = labels)}

binary_plot_data <-
  binary_plot_data %>%
  mutate(
    Variable = label_variable(
      Variable,
      binary_labels),
    Source = factor(
      Source,
      levels = source_levels))

binary_overall <-
  binary_overall %>%
  mutate(
    Variable = label_variable(
      Variable,
      binary_labels))

continuous_plot_data <-
  continuous_plot_data %>%
  mutate(
    Variable = label_variable(
      Variable,
      continuous_labels),
    Source = factor(
      Source,
      levels = source_levels))

continuous_overall <-
  continuous_overall %>%
  mutate(
    Variable = label_variable(
      Variable,
      continuous_labels))



# Common plotting theme
theme_sources <-
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 10.5),
    axis.text.y = element_text(
      size = 8),
    axis.text.x = element_text(
      size = 8),
    panel.spacing = unit(
      1.1,
      "lines"))

# Binary panel
p_binary <-
  ggplot(
    binary_plot_data,
    aes(
      x = Value,
      y = Source,
      colour = Source)) +
  geom_vline(
    data = binary_overall,
    aes(xintercept = Overall),
    colour = "grey50",
    linewidth = 0.7,
    linetype = "dashed") +
  geom_point(
    size = 2.7,
    show.legend = FALSE) +
  facet_wrap(
    ~Variable,
    ncol = 3,
    axes = "all",
    axis.labels = "all") +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%")) +
  scale_colour_manual(
    values = sourceColorPalette) +
  labs( x = "Percentage", y = NULL) +
  theme_sources +
  theme(
    axis.title.x = element_text(
      face = "bold",size = 10))


# Continuous panel
facet_x_labels <- data.frame(
  Variable = factor(
    c(
      "Age at Rupture",
      "IA Size at Rupture"),
    levels = c(
      "Age at Rupture",
      "IA Size at Rupture")),
  x = c(55, 2.5),
  y = -Inf,
  label = c(
    "Years",
    "log mm"))

p_continuous <-
  ggplot(
    continuous_plot_data,
    aes(
      x = Value,
      y = Source,
      fill = Source)) +
  geom_vline(
    data = continuous_overall,
    aes(xintercept = Overall),
    colour = "#4D4D4D",
    linewidth = 0.6,
    linetype = "dashed") +
  geom_boxplot(
    width = 0.55,
    linewidth = 0.6,
    outlier.size = 1.2,
    outlier.alpha = 0.7) +
  facet_grid(
    ~Variable,
    scales = "free_x",
    axes = "all") +
  scale_fill_manual(
    values = sourceColorPalette,
    guide = "none") +
  scale_colour_manual(
    values = sourceColorPalette,
    guide = "none") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(
        limits = c(10, 95),
        breaks = seq(10, 95, by = 20),
        expand = c(0, 0)),
      scale_x_continuous(
        limits = c(-1, 5),
        breaks = seq(-1, 5, by = 1),
        expand = c(0, 0)))) +
  geom_text(
    data = facet_x_labels,
    aes(
      x = x,
      y = y,
      label = label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.5,
    vjust = 4) +
  coord_cartesian(
    clip = "off") +
  labs(
    x = NULL,
    y = NULL) +
  theme_sources +
  theme(
    axis.ticks.y = element_blank())

# Combine panels
p_combined <- patchwork::wrap_plots(
  p_binary,
  p_continuous,
  ncol = 1,
  heights = c(1, 1)) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_combined

# Save figure
ggsave(
  p_combined,
  filename = "01_Baseline_characteristics_by_source.pdf",
  path = "figures",
  width = 10,
  height = 8,
  units = "in",
  dpi = 1200)

###############################################################################-
# Fig. 1b: ancestry PCs 1 + 2 ####
###############################################################################-
p <- ggplot()+
  geom_point(data=dat_s, aes(x=PC1_AVG_Z, y=PC2_AVG_Z, color=Source), alpha=0.6,
             show.legend = F)+
  labs(title="", x="PC1", y="PC2") +
  theme_classic(base_size = 14) +
  scale_color_manual(values = sourceColorPalette)+
  theme(legend.position = "bottom")

p


ggsave(
  p,
  filename = "01_PC12.pdf",
  path = "figures",
  width = 8,
  height = 5,
  units = "in",
  dpi = 1200)

###############################################################################-
#### Fig. 1c: PRS distributions per source ####
###############################################################################-
plots <- map("PRStotal_Z", function(var) {
  facet_labels <- dat_s %>%
    count(Source) %>%
    mutate(
      Source_label = paste0(Source, "\n(n = ", n, ")"),
      Source_label = factor(
        Source_label,
        levels = paste0(
          levels(dat_s$Source),
          "\n(n = ",
          n[match(levels(dat_s$Source), Source)],")")))
  dat_plot <- dat_s %>%
    left_join(facet_labels, by = "Source")
  
  ggplot(
    dat_plot,
    aes(
      x = .data[[var]],
      fill = Source,
      colour = Source)) +
    geom_density(
      alpha = 1,
      linewidth = 0.6) +
    facet_wrap(
      ~ Source_label,
      scales = "fixed") +
    scale_fill_manual(
      values = sourceColorPalette) +
    scale_colour_manual(
      values = sourceColorPalette) +
    labs(x = "Polygenic risk score", y = "Density") +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      strip.text = element_text(
        size = 14,
        face = "bold"),
      strip.background = element_blank(),
      panel.spacing = grid::unit(1, "lines"),
      legend.position = "none")})

plots

ggsave(
  filename = "01_PRS_distribution_by_source.pdf",
  plot = plots[[1]],
  path = "figures",
  width = 8,
  height = 5,
  units = "in",
  device = cairo_pdf)


###############################################################################-
#### Fig. A1a: AC distribution ####
###############################################################################-
plots_AC <- map("allele_count", function(var) {
  facet_labels <- dat_s %>%
    count(Source) %>%
    mutate(
      Source_label = paste0(Source, "\n(n = ", n, ")"),
      Source_label = factor(
        Source_label,
        levels = paste0(
          levels(dat_s$Source),
          "\n(n = ",
          n[match(levels(dat_s$Source), Source)], ")")))
  
  dat_plot <- dat_s %>%
    left_join(facet_labels, by = "Source")
  
  ggplot(
    dat_plot,
    aes(
      x = .data[[var]],
      fill = Source,
      colour = Source)) +
    geom_density(
      alpha = 1,
      linewidth = 0.6) +
    facet_wrap(
      ~ Source_label,
      nrow = 4,
      ncol = 3,
      scales = "fixed")+
    scale_fill_manual(
      values = sourceColorPalette) +
    scale_colour_manual(
      values = sourceColorPalette) +
    labs(x = "Allele Count", y = "Density") +
    theme_classic(base_size = 9) +
    theme(
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 9),
      strip.text = element_text(
        size = 8,
        face = "bold"),
      strip.background = element_blank(),
      panel.spacing = grid::unit(1, "lines"),
      legend.position = "none")})

plots_AC

ggsave(
  filename = "Appendix_AC_distribution_by_source.pdf",
  plot = plots_AC,
  path = "figures",
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf)

###############################################################################-
#### Fig. A1b: PRS vs. AC correlation ####
###############################################################################-
cor_test <- cor.test(
  dat_s$PRStotal_Z,
  dat_s$allele_count_Z,
  method = "pearson",
  use = "complete.obs")

r <- unname(cor_test$estimate)
p <- cor_test$p.value

ggplot(dat_s) +
  geom_point(
    aes(x = PRStotal_Z, y = allele_count),
    alpha = 0.6) +
  geom_smooth(
    aes(x = PRStotal_Z, y = allele_count),
    method = "lm",
    color = "black",
    se = TRUE) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = sprintf("r = %.2f", r),
    hjust = -0.2,
    vjust = 1.5,
    size = 4) +
  labs(
    x = "PRS",
    y = "AC") +
  theme_minimal()

ggsave(
  filename = "Appendix_PRS_AC_correlation.pdf",
  path = "figures",
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf,
  dpi=1200)

###############################################################################- 
#### Fig. A2: ancestry PC correlation matrix ####
###############################################################################-
# Select only the relevant columns
cluster_data <- dat_s[,c("PC1_AVG_Z", "PC2_AVG_Z", "PC3_AVG_Z", 
                         "PC4_AVG_Z", "PC5_AVG_Z")]

# Calculate the correlation matrix
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
  outline.color = "white") +
  theme_classic() +
  xlab("") +
  ylab("") +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")

p

ggsave(p,
       filename = "Appendix_PCs_correlation.pdf",
       path = "figures",
       width = 5,
       height = 5,
       units = "in",
       dpi = 1200)


###############################################################################-
#### Fig. A3: ancestry PC 1-5 ####
###############################################################################-
plots <- lapply(2:5, function(i) {
  pc_y <- paste0("PC", i, "_AVG_Z")
  ggplot(dat_s, aes(
    x = PC1_AVG_Z,
    y = .data[[pc_y]],
    color = Source)) +
    geom_point(size = 2, alpha = 0.8) +
    
    stat_ellipse(
      aes(group = Source),
      linewidth = 0.6,
      alpha = 0.5) +
    scale_color_manual(values = sourceColorPalette, name = "") +
    theme_minimal(base_size = 12) +
    labs(
      x = "PC1",
      y = paste0("PC", i),
      color = "Source") +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 20))})

# Arrange in 2 × 2 grid with one shared legend
p_combined <- patchwork::wrap_plots(
  plots,
  ncol = 2) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_combined

ggsave(
  p_combined,
  filename = "Appendix_PCs.pdf",
  path = "figures",
  width = 10,
  height = 8,
  units = "in",
  dpi = 1200)

###############################################################################-
#### Fig: 3-D representation of PC1-3 ####
###############################################################################-
plot_ly(dat_s, x = ~PC1_AVG_Z, y = ~PC2_AVG_Z, z= ~PC3_AVG_Z,
        color= ~Source,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 4),
        colors = sourceColorPalette) %>%
  layout(title = "3D Illustration of the first three PCs for simulated RIA patients",
         scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3")))


################################################################################-
#### ICC computations ####
################################################################################-
keep <- c(
  "Source",
  "Sex",
  "Positive.famillial.history",
  "Hypertension",
  "Smoking_Current_Former_No",
  "IAlocation_cow",
  "Multiple.IAs",
  "PRStotal_Z",
  "allele_count_Z",
  "PC1_AVG_Z", 
  "PC2_AVG_Z", 
  "PC3_AVG_Z", 
  "PC4_AVG_Z", 
  "PC5_AVG_Z")

dat <- dat_s %>%
  select(all_of(keep)) %>%
  na.omit() %>%
  mutate(Source = droplevels(as.factor(Source)))

vars <- keep[keep != "Source"]

# Convert factors to numeric values to enable variance-component
# estimation from linear mixed models used for ICC calculation.
for (v in vars) {
  if (is.factor(dat[[v]])) {
    if (nlevels(dat[[v]]) == 2) {
      dat[[v]] <- as.numeric(dat[[v]] == levels(dat[[v]])[2])
    } else {
      dat[[v]] <- as.numeric(dat[[v]])
    }
  }
}

results <- lapply(vars, function(v) {
  f <- as.formula(paste(v, "~ 1 + (1 | Source)"))
  m <- lmer(f, data = dat)
  vc <- as.data.frame(VarCorr(m))
  
  between <- vc$vcov[vc$grp == "Source"]
  within  <- vc$vcov[vc$grp == "Residual"]
  icc <- between / (between + within)
  
  data.frame(
    variable = v,
    between_var = round(between,3),
    within_var = round(within,3),
    icc = round(icc,3))})

results <- dplyr::bind_rows(results)

results %>%
  arrange(desc(icc)) %>%
  gt() %>%
  cols_label(
    variable = "Variable",
    between_var = "Between Var",
    within_var = "Within Var",
    icc = "ICC" ) %>%
  data_color(
    columns = icc,
    fn = scales::col_numeric(
      palette = c("white", "steelblue"),
      domain = NULL)) %>%
  tab_header(
    title = "Intraclass Correlation Coefficients (ICC)",
    subtitle = "Random intercept models by Source")
