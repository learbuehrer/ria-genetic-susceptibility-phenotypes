###############################################################################-
# File Content ####
# File: A_functions.R
#
# Purpose:
#   Collection of reusable helper functions used throughout the project.
#
# Main Components:
#   1. Table generation
#   2. Model specification and fitting
#   3. Analysis pipelines
#   4. Model diagnostics
#   5. Forest plots
#   6. Combine Forest plots
#   7. Sensitivity analyses fins
#   8. Sensitivity analyses Forest plots
#   9. Source-specific analyses
#   10. Random-effects visualisation
#   11. Multiple testing correction
#
# Notes:
#   Variables ending in "_Z" are z-standardised variables.
#
# Dependencies:
#   This file is sourced by 00_setup.R and subsequently used by all
#   analysis scripts.
#
###############################################################################-


###############################################################################-
# 1. TABLE GENERATION ####
###############################################################################-
##  Formatting helper for LaTeX bold text
bold <- function(x) {
  paste0('{\\textbf{', x, '}}')}

## Table 1 for R Markdown / Quarto
## Creates a TableOne table and automatically returns HTML or LaTeX
## HTML Table 1
tab_df_html <- function(vars, vars_cat, vars_nonnorm, strat = NULL, dat,
                        caption = NULL, SMD = TRUE, overall = TRUE,
                        p = FALSE, Levels = TRUE, Missing = TRUE,
                        digits = 1) {
  
  ## Create TableOne object
  if (is.null(strat)) {
    tab <- CreateTableOne(
      vars = vars,
      factorVars = vars_cat,
      data = dat,
      includeNA = TRUE
    )
  } else {
    tab <- CreateTableOne(
      vars = vars,
      factorVars = vars_cat,
      data = dat,
      strata = strat,
      addOverall = overall,
      includeNA = TRUE
    )
  }
  
  ## Print table
  tab_mat <- print(
    tab,
    nonnormal = vars_nonnorm,
    smd = SMD && !is.null(strat),
    contDigits = digits,
    showAllLevels = Levels,
    missing = Missing,
    printToggle = FALSE
  )
  
  ## Convert to data.frame while preserving row names
  tab_df <- cbind(
    Variable = rownames(tab_mat),
    as.data.frame(tab_mat, check.names = FALSE)
  )
  
  rownames(tab_df) <- NULL
  
  ## Clean up row labels
  tab_df[[1]] <- stringr::str_replace_all(tab_df[[1]], "   ", "  ")
  tab_df[[1]] <- stringr::str_replace_all(tab_df[[1]], "NA", "Missing")
  
  ## Replace missing values only in character columns
  tab_df[] <- lapply(tab_df, function(x) {
    if (is.character(x)) tidyr::replace_na(x, "Missing") else x
  })
  
  ## Rename columns
  names(tab_df)[names(tab_df) == "Missing"] <- "Missing (%)"
  names(tab_df)[names(tab_df) == "p"] <- "p-val"
  names(tab_df)[names(tab_df) == "level"] <- "Levels"
  
  ## Remove columns if requested
  remove_cols <- if (!p) {
    intersect(c("p", "test"), names(tab_df))
  } else {
    intersect("test", names(tab_df))
  }
  
  if (length(remove_cols) > 0) {
    tab_df <- tab_df[, !names(tab_df) %in% remove_cols, drop = FALSE]
  }
  
  ## Create HTML table
  kableExtra::kbl(
    tab_df,
    caption = caption,
    escape = FALSE,
    align = rep("l", ncol(tab_df))
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c(
        "striped",
        "hover",
        "condensed",
        "responsive"
      ),
      full_width = FALSE)}

## Table 1 as LaTeX data frame
tab_df <- function(vars, vars_cat, vars_nonnorm, strat, dat,
                   caption, label, filename, sideways = "table",
                   fontsize = "footnotesize", SMD = TRUE,
                   overall = TRUE, p = FALSE, Levels = TRUE,
                   Missing = TRUE) {
  
  tab <- CreateTableOne(
    vars = vars,
    factorVars = vars_cat,
    data = dat,
    strata = strat,
    addOverall = overall,
    includeNA = TRUE)
  
  tab <- print(
    tab,
    nonnormal = vars_nonnorm,
    smd = SMD,
    contDigits = 1,
    showAllLevels = Levels,
    missing = Missing)
  
  ## Sanitize row names for LaTeX
  rownames(tab) <- str_replace_all(
    rownames(tab), "   ", "\\\\quad ")
  rownames(tab) <- str_replace_all(
    rownames(tab), "%", "\\\\%")
  rownames(tab) <- str_replace_all(
    rownames(tab), "_", "\\\\_")
  rownames(tab) <- str_replace_all(
    rownames(tab), "NA", "Missing")
  
  tab[, 1] <- str_replace_na(tab[, 1], replacement = "Missing")
  
  colnames(tab)[colnames(tab) == "Missing"] <- "Missing (\\%)"
  colnames(tab)[colnames(tab) == "p"] <- "p-val"
  colnames(tab)[colnames(tab) == "level"] <- "Levels"
  
  if (!p) {
    ind <- which(colnames(tab) %in% c("p", "test"))
  } else {
    ind <- which(colnames(tab) == "test")}
  
  tab.df <- cbind(
    Variable = rownames(tab),
    tab[, -ind])
  
  ## Create LaTeX table
  xtab.tab <- xtable(
    tab.df,
    caption = caption,
    label = label,
    align = paste0("ll|", strrep("r", ncol(tab.df) - 1)))
  
  hlines <- c(-1, 0, 1, nrow(xtab.tab))
  
  print(
    xtab.tab,
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.colnames.function = bold,
    hline.after = hlines,
    caption.placement = "top",
    floating = TRUE,
    type = "latex",
    floating.environment = sideways,
    size = paste0("\\", fontsize),
    file = filename)}

###############################################################################-
# 2. MODEL SPECIFICATION ####
###############################################################################-
## Define regression formulas according to analysis level
build_formula <- function(
    outcome,
    level = c(
      "base",
      "pcs",
      "clinical",
      "famhist",
      "riskfactors",
      "personal"),
    random = NULL,
    clinical = NULL) {
  
  level <- match.arg(level)
  
  ## Base model: genetic score + Sex
  rhs <- exp_variables
  
  ## Add PCs
  if (level %in% c(
    "pcs",
    "clinical",
    "famhist",
    "riskfactors",
    "personal")) {
    rhs <- c(rhs, paste0("PC", 1:5, "_AVG_Z"))} # use the first 5 PCs
  
  ## Add outcome-specific clinical variable
  if (level == "clinical") {
    
    if (is.null(clinical)) {
      stop("For level = 'clinical', please supply 'clinical'.")
    }
    
    rhs <- c(rhs, clinical)
  }
  
  ## Add family history
  if (level %in% c("famhist", "personal")) {
    
    rhs <- c(
      rhs,
      "Positive.famillial.history"
    )
  }
  
  
  ## Add smoking + hypertension
  if (level %in% c("riskfactors", "personal")) {
    
    rhs <- c(
      rhs,
      "Smoking_Current_Former_No",
      "Hypertension") }
  
  ## Remove duplicates
  rhs <- unique(rhs)
  
  ## Create RHS
  rhs_formula <- paste(rhs, collapse = " + ")
  
  ## Add random effect
  if (!is.null(random)) {
    rhs_formula <- paste0(
      rhs_formula,
      " + (", random, ")"
    )
  }
  
  ## Return formula
  as.formula(
    paste(outcome, "~", rhs_formula)
  )
}

## Fit a regression model
fit_model <- function(
    data,
    outcome,
    model_type,
    level = "base",
    random = NULL,
    clinical = NULL) {
  form <- build_formula(
    outcome = outcome,
    level = level,
    random = random,
    clinical = clinical
  )
  
  fit <- switch(
    model_type,
    
    lm = lm(
      formula = form,
      data = data),
    
    lmer = lmerTest::lmer(
      formula = form,
      data = data),
    
    glm = glm(
      formula = form,
      data = data,
      family = binomial),
    
    glmer = lme4::glmer(
      formula = form,
      data = data,
      family = binomial),
    
    multinom = nnet::multinom(
      formula = form,
      data = data,
      trace = FALSE,
      model = TRUE))
  
  return(fit)}

## Fit model using an explicitly supplied formula
fit_model_s <- function(data, outcome = NULL, model_type,
                        level = "base", random = NULL,
                        formula = NULL) {
  
  if (is.null(formula)) {
    form <- build_formula(
      outcome = outcome,
      level = level,
      random = random
    )
  } else {
    form <- formula
  }
  
  switch(
    model_type,
    
    lm = lm(formula = form, data = data),
    
    lmer = lmerTest::lmer(
      formula = form,
      data = data
    ),
    
    glm = glm(
      formula = form,
      data = data,
      family = binomial
    ),
    
    glmer = lme4::glmer(
      formula = form,
      data = data,
      family = binomial
    ),
    
    multinom = nnet::multinom(
      formula = form,
      data = data,
      trace = FALSE,
      model = TRUE))}


###############################################################################-
# 3. ANALYSIS PIPELINES  ####
###############################################################################-
## Model order
MODEL_NAMES <- c(
  "base",
  "pcs",
  "clinical",
  "famhist",
  "riskfactors",
  "personal"
)

## Run all hierarchical model specifications
run_pipeline <- function(
    data,
    outcome,
    model_type,
    random = NULL,
    clinical = NULL
) {
  
  models <- list(
    
    ## a = Genetic Score + Sex
    base = fit_model(
      data = data,
      outcome = outcome,
      model_type = model_type,
      level = "base",
      random = random
    ),
    
    ## b = base + PCs
    pcs = fit_model(
      data = data,
      outcome = outcome,
      model_type = model_type,
      level = "pcs",
      random = random
    ),
    
    ## c = PCs + outcome-specific clinical variable
    clinical = fit_model(
      data = data,
      outcome = outcome,
      model_type = model_type,
      level = "clinical",
      random = random,
      clinical = clinical
    ),
    
    ## d = PCs + family history
    famhist = fit_model(
      data = data,
      outcome = outcome,
      model_type = model_type,
      level = "famhist",
      random = random
    ),
    
    ## e = PCs + smoking + hypertension
    riskfactors = fit_model(
      data = data,
      outcome = outcome,
      model_type = model_type,
      level = "riskfactors",
      random = random
    ),
    
    ## f = PCs + family history + smoking + hypertension
    personal = fit_model(
      data = data,
      outcome = outcome,
      model_type = model_type,
      level = "personal",
      random = random
    )
  )
  
  return(models)
}


## Model diagnostics for a single model
model_diagnostics <- function(model) {
  
  ## AIC / BIC
  aic_val <- AIC(model)
  bic_val <- BIC(model)
  
  ## Marginal R2
  r2_marginal <- tryCatch(
    {
      performance::r2_nakagawa(model)$R2_marginal
    },
    error = function(e) NA,
    warning = function(w) NA
  )
  
  
  ## Return one-row data frame
  data.frame(
    AIC = aic_val,
    BIC = bic_val,
    R2_marginal = r2_marginal,
    stringsAsFactors = FALSE
  )
}

###############################################################################-
# 4. MODEL DIAGNOSTICS ####
###############################################################################-
## Tukey-Anscombe diagnostic plot
TukeyAnscombePlot <- function(model,
                              title = NULL) {
  
  ## Determine model type
  is_glmm <- inherits(model, "glmerMod")
  is_glm  <- inherits(model, "glm")
  is_lmer <- inherits(model, "lmerMod")
  
  
  ## Fitted values
  fitted_values <- fitted(model)
  
  
  ## Residuals
  ## Pearson residuals work for all three model types
  residuals_values <- residuals(
    model,
    type = "pearson"
  )
  
  
  ## Data
  plot_df <- data.frame(
    fitted = fitted_values,
    residuals = residuals_values
  )
  
  ## Axis labels
  if (is_lmer) {
    
    y_label <- "Pearson residuals"
    
  } else if (is_glm || is_glmm) {
    
    y_label <- "Pearson residuals"
    
  } else {
    
    stop(
      "Unsupported model class: ",
      paste(class(model), collapse = ", ")
    )
  }
  
  ## Plot
  ggplot(
    plot_df,
    aes(
      x = fitted,
      y = residuals
    )
  ) +
    
    geom_point(
      alpha = 0.6
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    
    geom_smooth(
      method = "loess",
      se = TRUE
    ) +
    
    theme_classic() +
    
    labs(
      title = title,
      x = "Fitted values",
      y = y_label
    )
}


TukeyAnscombe_pipeline <- function(models,
                                   ncol = 3) {
  
  labels <- letters[seq_along(models)]
  
  plots <- lapply(
    seq_along(models),
    function(i) {
      
      TukeyAnscombePlot(
        model = models[[i]],
        title = labels[i]
      )
      
    }
  )
  
  ggpubr::ggarrange(
    plotlist = plots,
    ncol = ncol,
    nrow = ceiling(length(plots) / ncol)
  )
}

diagnose_pipeline <- function(models) {
  ## Basic diagnostics
  diagnostics <- lapply(
    names(models),
    function(x) {
      out <- model_diagnostics(models[[x]])
      out$Model <- x
      out
    }
  )
  
  diagnostics <- dplyr::bind_rows(diagnostics) %>%
    dplyr::select(Model, dplyr::everything())
  
  ## Explicit likelihood-ratio comparisons
  ## Model numbers to compare:
  ## 2 vs 1  = PCs vs Base
  ## 3 vs 2  = Clinical vs Base
  ## 4 vs 2  = FamHist vs PCs
  ## 5 vs 2  = RiskFactors vs PCs
  ## 6 vs 5  = Personal vs RiskFactors
  ## 6 vs 4  = Personal vs FamHist
  comparisons <- data.frame(
    model_1 = c(2, 3, 4, 5, 6, 6),
    model_0 = c(1, 2, 2, 2, 5, 4))
  
  
  ## Run LR tests
  lr_results <- lapply(
    seq_len(nrow(comparisons)),
    function(i) {
      
      m1_id <- comparisons$model_1[i]
      m0_id <- comparisons$model_0[i]
      
      m1_name <- names(models)[m1_id]
      m0_name <- names(models)[m0_id]
      
      m0 <- models[[m0_id]]
      m1 <- models[[m1_id]]
      
      ## Likelihood-ratio test
      lr <- anova(m0, m1)
      
      ## Identify columns
      lr_stat_col <- which(
        names(lr) %in% c("Chisq", "L.Ratio")
      )
      
      lr_p_col <- which(
        names(lr) %in% c("Pr(>Chisq)", "p.value")
      )
      
      ## Extract second row
      LR_value <- if (length(lr_stat_col) > 0) {
        lr[[lr_stat_col[1]]][2]
      } else {
        NA_real_
      }
      
      LR_p_value <- if (length(lr_p_col) > 0) {
        lr[[lr_p_col[1]]][2]
      } else {
        NA_real_
      }
      
      data.frame(
        Model = m1_name,
        Compared_to = m0_name,
        LR = LR_value,
        LR_p = LR_p_value,
        stringsAsFactors = FALSE
      )
    }
  )
  
  lr_results <- dplyr::bind_rows(lr_results)
  
  
  ## Combine diagnostics + LR results
  diagnostics <- diagnostics %>%
    dplyr::left_join(
      lr_results,
      by = "Model"
    ) %>%
    dplyr::mutate(
      LR_p_formatted = biostatUZH::formatPval(LR_p)
    )
  
  return(
    list(
      table = diagnostics,
      models = models))}

## Normal Q-Q diagnostic plot for lmer models
QQPlot <- function(model,
                   title = NULL) {
  
  ## Check model type
  is_lmer <- inherits(model, "lmerMod")
  
  if (!is_lmer) {
    stop(
      "QQPlot() only supports lmer models. Model class: ",
      paste(class(model), collapse = ", ")
    )
  }
  
  
  ## Residuals
  residuals_values <- residuals(
    model,
    type = "pearson"
  )
  
  
  ## Data
  plot_df <- data.frame(
    residuals = residuals_values
  )
  
  
  ## Plot
  ggplot(
    plot_df,
    aes(
      sample = residuals
    )
  ) +
    
    stat_qq() +
    
    stat_qq_line() +
    
    theme_classic() +
    
    labs(
      title = title,
      x = "Theoretical quantiles",
      y = "Sample quantiles"
    )
}


QQ_pipeline <- function(models,
                        ncol = 3) {
  
  ## Keep only lmer models
  lmer_models <- models[
    vapply(
      models,
      function(x) inherits(x, "lmerMod"),
      logical(1)
    )
  ]
  
  if (length(lmer_models) == 0) {
    stop("No lmer models found in the supplied model list.")
  }
  
  labels <- letters[seq_along(lmer_models)]
  
  plots <- lapply(
    seq_along(lmer_models),
    function(i) {
      
      QQPlot(
        model = lmer_models[[i]],
        title = labels[i]
      )
      
    }
  )
  
  ggpubr::ggarrange(
    plotlist = plots,
    ncol = ncol,
    nrow = ceiling(length(plots) / ncol)
  )
}

###############################################################################-
# 5. FOREST PLOTS ####
###############################################################################-

## Dictionary for readable variable labels
term_dictionary <- c(
  "SexMale" = "Male",
  "allele_count_Z" = "Allele Count",
  "Positive.famillial.historyYes" = "Aware of Family History",
  "Positive.famillial.historyNo" = "No Family History",
  "Smoking_Current_Former_NoFormer" = "Former Smoker",
  "Smoking_Current_Former_NoCurrent" = "Current Smoker",
  "HypertensionAnyType" = "Aware of Hypertension",
  "IA_size.at.time.of.SAH.log_Z" = "IA Size (log Z)",
  "Age.at.time.of.SAH.log_Z" = "Age at Rupture (Z)",
  "IAlocation_cowCommunicating CoW segment" = "Communicating CoW Segment",
  "IAlocation_cowDistal peripheral branches" = "Distal Peripheral Branches",
  "IAlocation_cowICA" = "ICA",
  "IAlocation_cowMCA" = "MCA",
  "IAlocation_cowNon CoW ICA segments" = "Non CoW ICA Segments",
  "IAlocation_cowPosterior circulation" = "Posterior Circulation",
  "PRStotal_Z" = "PRS (Z)",
  "Multiple.IAsYes" = "Multiple IAs",
  "PC1_AVG_Z" = "PC1 (Z)",
  "PC2_AVG_Z" = "PC2 (Z)",
  "PC3_AVG_Z" = "PC3 (Z)",
  "PC4_AVG_Z" = "PC4 (Z)",
  "PC5_AVG_Z" = "PC5 (Z)")


## Standard forest plot for lm, lmer, glm and glmer models
ForestPlot <- function(model, title = NULL, exponentiate = NULL,
                       measure = NULL, col_def = "black",
                       sig_col = "firebrick", hide_intercept = TRUE,
                       hide_random = TRUE, term1 = "Lower",
                       term2 = "Higher", cohort = NULL,
                       digits_effect = 2, digits_p = 2) {
  
  ## Determine whether to exponentiate automatically
  if (is.null(exponentiate)) {
    exponentiate <- inherits(model, c("glm", "glmerMod"))
  }
  
  ## Axis label
  if (is.null(measure)) {
    measure <- ifelse(exponentiate, "OR", "β")
  }
  
  ## Extract model estimates and confidence intervals
  sum_df <- broom.mixed::tidy(
    model,
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = exponentiate
  )
  
  ## Add p-values for glmer models if necessary
  if (!"p.value" %in% names(sum_df)) {
    if (inherits(model, "glmerMod")) {
      sum_df <- sum_df %>%
        mutate(
          p.value = 2 * pnorm(
            abs(statistic),
            lower.tail = FALSE
          )
        )
    } else {
      stop("No p.value was returned for this model.")
    }
  }
  
  ## Significance based on unrounded p-value
  sum_df <- sum_df %>%
    mutate(
      sig = ifelse(
        !is.na(p.value) & p.value < 0.05,
        "Significant",
        "Not significant"
      )
    )
  
  ## Remove intercept
  if (hide_intercept) {
    sum_df <- filter(sum_df, term != "(Intercept)")
  }
  
  ## Remove random-effect parameters
  if (hide_random) {
    sum_df <- filter(
      sum_df,
      !grepl("sd__|cor__|\\.sigma", term)
    )
  }
  
  ## Readable variable labels
  sum_df$term_label <- dplyr::recode(
    sum_df$term,
    !!!term_dictionary,
    .default = sum_df$term
  )
  
  ## Create effect and p-value annotations
  sum_df <- sum_df %>%
    mutate(
      effect_text = sprintf(
        paste0("%.", digits_effect, "f"),
        estimate
      ),
      
      p_text = case_when(
        is.na(p.value) ~ "p = —",
        p.value < 0.0001 ~ "p < 0.0001",
        p.value < 0.05 ~ paste0(
          "p = ",
          formatC(
            p.value,
            format = "f",
            digits = 4
          )
        ),
        TRUE ~ paste0(
          "p = ",
          formatC(
            p.value,
            format = "f",
            digits = digits_p
          )
        )
      ),
      
      Label = paste0(
        term_label,
        " (",
        measure,
        " = ",
        effect_text,
        ", ",
        p_text,
        ")"
      )
    )
  
  sum_df$cohort <- cohort
  
  ## Plot limits
  xmin <- min(sum_df$conf.low, na.rm = TRUE)
  xmax <- max(sum_df$conf.high, na.rm = TRUE)
  xrange <- xmax - xmin
  padding <- 0.15 * xrange
  
  ## Factor labels
  sum_df$Label <- factor(
    sum_df$Label,
    levels = rev(sum_df$Label)
  )
  
  ## Forest plot
  p <- ggplot(
    sum_df,
    aes(
      x = estimate,
      y = Label,
      colour = sig
    )
  ) +
    geom_point(size = 3) +
    scale_linewidth_identity() +
    geom_errorbar(
      aes(
        ymin = Label,
        ymax = Label,
        xmin = conf.low,
        xmax = conf.high
      ),
      orientation = "y",
      width = 0.25
    )+
    
    geom_vline(
      xintercept = ifelse(
        exponentiate,
        1,
        0
      ),
      linetype = 2,
      colour = "grey50"
    ) +
    
    scale_colour_manual(
      values = c(
        "Significant" = sig_col,
        "Not significant" = col_def
      )
    ) +
    
    coord_cartesian(
      xlim = c(
        xmin - padding,
        xmax + padding
      ),
      clip = "off"
    ) +
    
    theme_classic(base_size = 10) +
    
    theme(
      axis.text.x = element_text(
        size = 10,
        margin = margin(t = 4)),
      
      axis.text.y = element_text(
        size = 10 ),
      
      axis.title.x = element_text(
        size = 10,
        margin = margin(t = 8)),
      
      axis.title.y = element_blank(),
      
      legend.position = "none",
      
      plot.margin = margin(
        t = 10,
        r = 20,
        b = 30,
        l = 20
      )
    ) +
    
    labs(
      x = measure,
      y = NULL,
      title = title
    )
  
  ## Directional labels
  p <- p +
    annotate(
      "text",
      x = xmin - padding / 2,
      y = 0,
      label = term1,
      hjust = 0.5,
      vjust = 3
    ) +
    
    annotate(
      "text",
      x = xmax + padding / 2,
      y = 0,
      label = term2,
      hjust = 0.5,
      vjust = 3
    )
  
  ## Return data and plot
  list(
    data = sum_df,
    plot = p
  )
}

###############################################################################-
# 6. COMBINED FOREST PLOTS ####
###############################################################################-
## Combine forest-plot data across cohorts
combine_forest_data <- function(model_complete, model_fins, model_nonfins) {
  complete <- ForestPlot(
    model_complete,
    cohort = "Complete cohort")$data
  
  fins <- ForestPlot(
    model_fins,
    cohort = "Fins")$data
  
  nonfins <- ForestPlot(
    model_nonfins,
    cohort = "Non-Fins")$data
  
  bind_rows(
    complete %>%
      mutate(
        Label = term_label,
        cohort = "Complete cohort"),
    
    fins %>%
      mutate(
        Label = term_label,
        cohort = "Fins"),
    
    nonfins %>%
      mutate(
        Label = term_label,
        cohort = "Non-Fins"))}


## Plot estimates across cohorts
plot_combined_forest <- function(
    dat,
    title = NULL,
    effect_label = NULL,
    exponentiate = FALSE,
    null = NULL,
    x_limits = NULL,
    term1 = "Lower",
    term2 = "Higher"
) {
  if (is.null(effect_label)) {
    effect_label <- ifelse(exponentiate, "OR", "β")
  }
  
  if (is.null(null)) {
    null <- ifelse(exponentiate, 1, 0)
  }
  
  variable_levels <- dat %>%
    filter(cohort == "Complete cohort") %>%
    pull(Label) %>%
    unique()
  
  variable_levels <- rev(variable_levels)
  
  dat <- dat %>%
    mutate(
      Label = factor(Label, levels = variable_levels),
      y = as.numeric(Label),
      y_position = case_when(
        cohort == "Complete cohort" ~ y + 0.20,
        cohort == "Non-Fins" ~ y,
        cohort == "Fins" ~ y - 0.20,
        TRUE ~ y
      ),
      significance = ifelse(
        p.value < 0.05,
        "Significant",
        "Not significant"
      )
    )
  
  if (is.null(x_limits)) {
    x_range <- range(
      c(dat$conf.low, dat$conf.high),
      na.rm = TRUE
    )
    
    padding <- diff(x_range) * 0.05
    
    if (padding == 0 || !is.finite(padding)) {
      padding <- 0.1
    }
    
    x_limits <- c(
      x_range[1] - padding,
      x_range[2] + padding
    )
  }
  
  p <- ggplot(
    dat,
    aes(
      x = estimate,
      colour = cohort,
      shape = significance
    )
  ) +
    geom_errorbar(
      aes(
        y = y_position,
        xmin = conf.low,
        xmax = conf.high
      ),
      orientation = "y",
      width = 0.25,
      linewidth = 0.6
    )+
    
    geom_point(
      aes(y = y_position),
      size = 3,
      stroke = 0.9
    ) +
    geom_vline(
      xintercept = null,
      linetype = 2,
      colour = "grey50"
    ) +
    scale_colour_manual(
      values = c(
        "Complete cohort" = "black",
        "Non-Fins" = "#4D6A7A",
        "Fins" = "#003580"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Significant" = 16,
        "Not significant" = 1
      )
    ) +
    scale_y_continuous(
      breaks = seq_along(variable_levels),
      labels = variable_levels,
      expand = expansion(
        mult = c(0.04, 0.04)
      )
    ) +
    coord_cartesian(
      xlim = x_limits,
      clip = "off"
    ) +
    labs(
      x = effect_label,
      y = NULL,
      colour = "Cohort",
      shape = "Significance",
      title = title
    ) +
    theme_classic(base_size = 10) +
    theme(
      text = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      plot.title = element_text(
        size = 11,
        face = "plain"
      ),
      legend.position = "bottom",
      plot.margin = margin(
        t = 8,
        r = 10,
        b = 42,
        l = 10
      )
    )
  
  ## Directional labels below the x-axis
  p <- p +
    annotation_custom(
      grid::textGrob(
        term1,
        x = grid::unit(0, "npc"),
        y = grid::unit(-0.13, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 10)
      )
    ) +
    annotation_custom(
      grid::textGrob(
        term2,
        x = grid::unit(1, "npc"),
        y = grid::unit(-0.13, "npc"),
        just = c("right", "top"),
        gp = grid::gpar(fontsize = 10)
      )
    )
  
  return(p)
}

CombinedPRSForestPlot <- function(
    prs_models,
    ac_models,
    prs_term = "PRStotal_Z",
    ac_term = "allele_count_Z",
    title = NULL,
    xlab = NULL,
    x_limits = NULL,
    clinical_label = NULL,
    prs_label = "PRS",
    ac_label = "Allele count",
    digits_effect = 2,
    digits_p = 4,
    term1 = "Lower",
    term2 = "Higher",
    show = c("PRS", "AC"),
    p_adjust_table = NULL,
    outcome_name = NULL,
    model_names = MODEL_NAMES
) {
  
  show <- match.arg(show)
  
  ## Check model lists
  if (!all(model_names %in% names(prs_models))) {
    stop(
      "prs_models must contain: ",
      paste(model_names, collapse = ", ")
    )
  }
  
  if (!all(model_names %in% names(ac_models))) {
    stop(
      "ac_models must contain: ",
      paste(model_names, collapse = ", ")
    )
  }
  
  
  ## Detect model type
  is_binomial <- inherits(
    prs_models$base,
    c("glm", "glmerMod")
  )
  
  ac_is_binomial <- inherits(
    ac_models$base,
    c("glm", "glmerMod")
  )
  
  if (is_binomial != ac_is_binomial) {
    stop(
      "PRS and allele-count models must use the same model family."
    )
  }
  
  
  ## Effect scale
  if (is_binomial) {
    
    exponentiate <- TRUE
    null_value <- 1
    
    if (is.null(xlab)) {
      xlab <- "OR"
    }
    
  } else {
    
    exponentiate <- FALSE
    null_value <- 0
    
    if (is.null(xlab)) {
      xlab <- "β"
    }
  }
  
  
  ## Determine clinical variable and display label
  if (is.null(clinical_label)) {
    
    clinical_formula <- formula(
      prs_models$clinical
    )
    
    clinical_terms <- attr(
      terms(clinical_formula),
      "term.labels"
    )
    
    excluded_terms <- c(
      prs_term,
      "Sex",
      paste0("PC", 1:5, "_AVG_Z")
    )
    
    ## Remove random-effect terms such as "1 | Source"
    clinical_terms <- clinical_terms[
      !grepl("\\|", clinical_terms)
    ]
    
    ## Identify the clinical variable
    clinical_var <- setdiff(
      clinical_terms,
      excluded_terms
    )
    
    if (length(clinical_var) == 0) {
      clinical_label <- "Clinical"
    } else {
      clinical_var <- clinical_var[1]
      
      ## Convert internal variable names to publication labels
      clinical_label <- dplyr::recode(
        clinical_var,
        "Multiple.IAs" = "IA Multiplicity",
        "IAlocation_cow" = "IA Location",
        .default = clinical_var
      )
    }
  }
  
  
  ## Y-axis labels
  y_labels <- c(
    "Genetics\n + Sex",
    "+ PCs",
    paste0("+ ", clinical_label),
    "+ Family History",
    "+ Smoking\n+ Hypertension",
    "+ Family History\n + Smoking\n+ Hypertension"
  )
  
  
  ## Extract one coefficient
  extract_term <- function(
    model,
    term_name,
    outcome_label
  ) {
    
    dat <- broom.mixed::tidy(
      model,
      effects = "fixed",
      conf.int = TRUE,
      exponentiate = exponentiate
    )
    
    dat <- dat %>%
      dplyr::filter(
        .data$term == term_name
      )
    
    if (nrow(dat) == 0) {
      stop(
        "Term '", term_name,
        "' was not found in one of the models."
      )
    }
    
    dat %>%
      dplyr::slice(1) %>%
      dplyr::transmute(
        estimate = estimate,
        conf.low = conf.low,
        conf.high = conf.high,
        p.value = p.value,
        outcome = outcome_label
      )
  }
  
  
  ## Extract PRS
  prs_dat <- purrr::map_dfr(
    model_names,
    function(m) {
      
      out <- extract_term(
        model = prs_models[[m]],
        term_name = prs_term,
        outcome_label = prs_label
      )
      
      out$model <- m
      out
    }
  )
  
  if (!is.null(p_adjust_table)) {
    
    prs_adj <- p_adjust_table %>%
      dplyr::filter(
        outcome == outcome_name
      ) %>%
      dplyr::select(
        model,
        p.adj
      )
    
    prs_dat <- prs_dat %>%
      dplyr::left_join(
        prs_adj,
        by = "model"
      )
    
  } else {
    
    prs_dat$p.adj <- prs_dat$p.value
  }
  
  
  ## Extract allele count
  ac_dat <- purrr::map_dfr(
    model_names,
    function(m) {
      
      out <- extract_term(
        model = ac_models[[m]],
        term_name = ac_term,
        outcome_label = ac_label
      )
      
      out$model <- m
      out
    }
  )
  
  if (!is.null(p_adjust_table)) {
    
    ac_adj <- p_adjust_table %>%
      dplyr::filter(
        outcome == outcome_name
      ) %>%
      dplyr::select(
        model,
        p.adj
      )
    
    ac_dat <- ac_dat %>%
      dplyr::left_join(
        ac_adj,
        by = "model"
      )
    
  } else {
    
    ac_dat$p.adj <- ac_dat$p.value
  }
  
  
  ## Select PRS or allele count
  if (show == "PRS") {
    
    plot_dat <- prs_dat
    
  } else if (show == "AC") {
    
    plot_dat <- ac_dat
    
  } else {
    
    plot_dat <- dplyr::bind_rows(
      prs_dat,
      ac_dat
    )
  }
  
  
  ## Prepare plotting data
  plot_dat <- plot_dat %>%
    dplyr::mutate(
      main_model = model == "pcs",
      
      model = factor(
        model,
        levels = model_names
      ),
      point_size = ifelse(main_model, 4.5, 3),
      ci_width   = ifelse(main_model, 1.1, 0.6),
      
      outcome = factor(
        outcome,
        levels = if (show == "PRS") {
          prs_label
        } else {
          ac_label
        }
      ),
      
      ## Base = 6 (top)
      ## Personal = 1 (bottom)
      model_y = 7 - as.numeric(model),
      
      ## Effect text
      effect_text = sprintf(
        paste0(
          ifelse(
            is_binomial,
            "OR = %.",
            "β = %."
          ),
          digits_effect,
          "f"
        ),
        estimate
      ),
      
      ## P-value text
      p_text = dplyr::case_when(
        
        is.na(p.adj) ~
          "pHolm = —",
        
        p.adj < 0.001 ~
          "pHolm < 0.001",
        
        TRUE ~ paste0(
          "pHolm = ",
          formatC(
            p.adj,
            format = "f",
            digits = digits_p
          )
        )
      ),
      
      ## Combined annotation
      result_text = paste(
        effect_text,
        p_text,
        sep = ", "
      ),
      
      ## Significance
      significance = dplyr::if_else(
        p.adj < 0.05,
        "Significant",
        "Not significant"
      )
    )
  
  
  ## Determine x-axis limits
  if (is.null(x_limits)) {
    
    x_range <- range(
      c(
        plot_dat$conf.low,
        plot_dat$conf.high,
        null_value
      ),
      na.rm = TRUE
    )
    
    xrange <- diff(x_range)
    
    if (xrange == 0) {
      xrange <- ifelse(
        is_binomial,
        0.5,
        1
      )
    }
    
    x_limits <- c(
      x_range[1] - 0.05 * xrange,
      x_range[2] + 0.05 * xrange
    )
    
    x_limits[1] <- min(
      x_limits[1],
      null_value
    )
    
    x_limits[2] <- max(
      x_limits[2],
      null_value
    )
  }
  
  
  ## Position effect/p text outside plotting region
  xrange <- diff(x_limits)
  text_x <- x_limits[2] + 0.02 * xrange
  
  
  ## Plot
  p <- ggplot(
    plot_dat,
    aes(
      x = estimate,
      y = model_y
    )
  ) +
    
    ## Separator between + PCs and additional adjustments
    annotate(
      "segment",
      x = x_limits[1],
      xend = x_limits[2],
      y = 4.5,
      yend = 4.5,
      linewidth = 0.5,
      colour = "grey50",
      linetype = 3
    )+
    
    ## Confidence intervals
    geom_segment(
      aes(
        x = conf.low,
        xend = conf.high,
        y = model_y,
        yend = model_y,
        linewidth = ci_width
      )
    ) +
    scale_linewidth_identity() +
    
    ## Point estimates
    geom_point(
      aes(
        colour = outcome,
        shape = significance
      ),
      size = 3,
      stroke = 0.9
    ) +
    
    ## β/OR + p-value
    geom_text(
      aes(
        x = text_x,
        label = result_text,
        fontface = ifelse(
          model == "pcs",
          "bold",
          "plain"
        )
      ),
      hjust = 0,
      size = 2.8,
      show.legend = FALSE
    ) +
    
    ## Null/reference line
    geom_vline(
      xintercept = null_value,
      linetype = 2,
      colour = "grey50"
    ) +
    
    ## PRS vs allele count
    scale_colour_manual(
      values = stats::setNames(
        "black",
        if (show == "PRS") {
          prs_label
        } else {
          ac_label
        }
      ),
      breaks = if (show == "PRS") {
        prs_label
      } else {
        ac_label
      },
      labels = if (show == "PRS") {
        "PRS"
      } else {
        "Allele count"
      },
      name = NULL
    ) +
    
    ## Significant vs non-significant
    scale_shape_manual(
      values = c(
        "Significant" = 16,
        "Not significant" = 1
      ),
      breaks = c(
        "Significant",
        "Not significant"
      ),
      labels = c(
        "p < 0.05",
        "p ≥ 0.05"
      ),
      name = NULL
    ) +
    
    ## Y-axis
    scale_y_continuous(
      breaks = 1:6,
      labels = rev(y_labels),
      expand = expansion(
        mult = c(
          0.04,
          0.04
        )
      )
    ) +
    
    ## X-axis
    coord_cartesian(
      xlim = x_limits,
      clip = "off"
    ) +
    
    labs(
      x = xlab,
      y = NULL,
      title = title,
      colour = NULL,
      shape = NULL
    ) +
    
    ## Directional labels
    annotation_custom(
      grid::textGrob(
        term1,
        x = grid::unit(0, "npc"),
        y = grid::unit(-0.08, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 10)
      )
    ) +
    
    annotation_custom(
      grid::textGrob(
        term2,
        x = grid::unit(1, "npc"),
        y = grid::unit(-0.08, "npc"),
        just = c("right", "top"),
        gp = grid::gpar(fontsize = 10)
      )
    ) +
    
    ## Theme
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      legend.position = "bottom",
      
      plot.title = element_text(
        size = 14,
        face = "plain"
      ),
      
      plot.margin = margin(
        t = 4,
        r = 100,
        b = 15,
        l = 4
      )
    ) +
    
    guides(
      colour = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    )
  
  
  ## Return
  return(
    list(
      data = plot_dat,
      plot = p
    )
  )
}

###############################################################################-
# 7. SENSITIVITY ANALYSIS FINS ####
###############################################################################-

run_cohort_comparison <- function(
    cohort_data,
    outcome,
    model_type,
    clinical,
    term1,
    term2,
    effect_label = NULL,
    null = NULL,
    model_names = MODEL_NAMES,
    plot_titles = c(
      base = "Base",
      pcs = "PCs",
      clinical = "Clinical",
      famhist = "Family History Awareness",
      riskfactors = "Risk Factors",
      personal = "Personal"
    )
) {
  
  ## Fit models
  results <- lapply(
    cohort_data,
    run_pipeline,
    outcome = outcome,
    model_type = model_type,
    random = "1|Source",
    clinical = clinical
  )
  
  ## Combine forest data
  forest_data <- lapply(
    model_names,
    function(x) {
      combine_forest_data(
        results[["Complete cohort"]][[x]],
        results[["Fins"]][[x]],
        results[["Non-Fins"]][[x]]
      )
    }
  )
  
  names(forest_data) <- model_names
  
  ## Common x-axis limits
  all_forest_data <- bind_rows(forest_data)
  
  x_range <- range(
    c(
      all_forest_data$conf.low,
      all_forest_data$conf.high
    ),
    na.rm = TRUE
  )
  
  padding <- diff(x_range) * 0.05
  
  x_limits <- c(
    x_range[1] - padding,
    x_range[2] + padding
  )
  
  ## Create plots
  plots <- Map(
    function(dat, title) {
      
      plot_combined_forest(
        dat,
        title = title,
        x_limits = x_limits,
        effect_label = effect_label,
        null = null,
        term1 = term1,
        term2 = term2
      )
      
    },
    forest_data,
    plot_titles
  )
  
  ## Final figure
  figure <- ggpubr::ggarrange(
    plotlist = plots,
    ncol = 2,
    nrow = 3,
    heights = c(1.5, 2, 2),
    common.legend = TRUE,
    legend = "bottom"
  )
  
  list(
    results = results,
    forest_data = forest_data,
    figure = figure
  )
}

###############################################################################-
# 8. SENSITIVITY FOREST PLOTS ####
###############################################################################-

SensitivityPRSForestPlot <- function(
    main_models,
    nonfins_models,
    fins_models,
    prs_term = "PRStotal_Z",
    title = NULL,
    xlab = NULL,
    x_limits = NULL,
    clinical_label = NULL,
    main_label = "Complete cohort",
    nonfins_label = "Non-Fins",
    fins_label = "Fins",
    digits_effect = 2,
    digits_p = 3,
    term1 = "Lower",
    term2 = "Higher",
    model_names = MODEL_NAMES) {
  
  ## Check model lists
  if (!all(model_names %in% names(main_models))) {
    stop(
      "main_models must contain: ",
      paste(model_names, collapse = ", ")
    )
  }
  
  if (!all(model_names %in% names(nonfins_models))) {
    stop(
      "nonfins_models must contain: ",
      paste(model_names, collapse = ", ")
    )
  }
  
  if (!all(model_names %in% names(fins_models))) {
    stop(
      "fins_models must contain: ",
      paste(model_names, collapse = ", ")
    )
  }
  
  
  ## Detect model type
  is_binomial <- inherits(
    main_models$base,
    c("glm", "glmerMod")
  )
  
  nonfins_is_binomial <- inherits(
    nonfins_models$base,
    c("glm", "glmerMod")
  )
  
  fins_is_binomial <- inherits(
    fins_models$base,
    c("glm", "glmerMod")
  )
  
  if (!(
    is_binomial ==
    nonfins_is_binomial &&
    is_binomial ==
    fins_is_binomial
  )) {
    stop(
      "All model sets must use the same model family."
    )
  }
  
  
  ## Effect scale
  if (is_binomial) {
    
    exponentiate <- TRUE
    null_value <- 1
    
    if (is.null(xlab)) {
      xlab <- "OR"
    }
    
  } else {
    
    exponentiate <- FALSE
    null_value <- 0
    
    if (is.null(xlab)) {
      xlab <- "β"
    }
  }
  
  ## Determine clinical label
  if (is.null(clinical_label)) {
    
    clinical_formula <- formula(main_models$clinical)
    
    clinical_terms <- attr(
      terms(clinical_formula),
      "term.labels"
    )
    
    excluded_terms <- c(
      prs_term,
      "Sex",
      paste0("PC", 1:5, "_AVG_Z")
    )
    
    ## Remove random-effect terms
    clinical_terms <- clinical_terms[
      !grepl("\\|", clinical_terms)
    ]
    
    ## Identify the clinical variable
    clinical_var <- setdiff(
      clinical_terms,
      excluded_terms)
    
    if (length(clinical_var) == 0) {
      clinical_label <- "Clinical"
    } else {
      clinical_label <- clinical_var[1]
    }
  }
  
  ## Convert internal variable names to publication labels
  clinical_label <- dplyr::case_when(
    clinical_label == "Multiple.IAs" ~ "IA Multiplicity",
    clinical_label == "IAlocation" ~ "IA Location",
    TRUE ~ clinical_label )
  
  ## Y-axis labels
  y_labels <- c(
    "Genetics\n + Sex",
    "+ PCs",
    paste0("+ ", clinical_label),
    "+ Family History",
    "+ Smoking\n+ Hypertension",
    "+ Family History\n + Smoking\n+ Hypertension"
  )
  
  ## Extract one PRS coefficient
  extract_term <- function(
    model,
    term_name,
    cohort_label
  ) {
    
    dat <- broom.mixed::tidy(
      model,
      effects = "fixed",
      conf.int = TRUE,
      exponentiate = exponentiate
    )
    
    dat <- dat %>%
      dplyr::filter(
        .data$term == term_name
      )
    
    if (nrow(dat) == 0) {
      stop(
        "Term '", term_name,
        "' was not found in one of the models."
      )
    }
    
    dat %>%
      dplyr::slice(1) %>%
      dplyr::transmute(
        estimate = estimate,
        conf.low = conf.low,
        conf.high = conf.high,
        p.value = p.value,
        cohort = cohort_label
      )
  }
  
  
  ## Extract complete cohort
  main_dat <- purrr::map_dfr(
    model_names,
    function(m) {
      
      out <- extract_term(
        model = main_models[[m]],
        term_name = prs_term,
        cohort_label = main_label
      )
      
      out$model <- m
      out
    }
  )
  
  
  ## Extract Non-Fins
  nonfins_dat <- purrr::map_dfr(
    model_names,
    function(m) {
      
      out <- extract_term(
        model = nonfins_models[[m]],
        term_name = prs_term,
        cohort_label = nonfins_label
      )
      
      out$model <- m
      out
    }
  )
  
  
  ## Extract Fins
  fins_dat <- purrr::map_dfr(
    model_names,
    function(m) {
      
      out <- extract_term(
        model = fins_models[[m]],
        term_name = prs_term,
        cohort_label = fins_label
      )
      
      out$model <- m
      out
    }
  )
  
  
  ## Combine
  plot_dat <- dplyr::bind_rows(
    main_dat,
    nonfins_dat,
    fins_dat
  ) %>%
    dplyr::mutate(
      
      model = factor(
        model,
        levels = model_names
      ),
      
      cohort = factor(
        cohort,
        levels = c(
          main_label,
          nonfins_label,
          fins_label
        )
      ),
      
      ## Base = 6 (top)
      ## Personal = 1 (bottom)
      model_y = 7 - as.numeric(model),
      
      ## Separate the three cohorts vertically
      y = dplyr::case_when(
        cohort == main_label ~ model_y + 0.20,
        cohort == nonfins_label ~ model_y,
        cohort == fins_label ~ model_y - 0.20
      ),
      
      ## Effect text
      effect_text = sprintf(
        paste0(
          ifelse(
            is_binomial,
            "OR = %.",
            "β = %."
          ),
          digits_effect,
          "f"
        ),
        estimate
      ),
      
      ## P-value text
      p_text = dplyr::case_when(
        is.na(p.value) ~ "p = —",
        p.value < 0.001 ~ "p < 0.001",
        p.value < 0.05 ~ paste0(
          "p = ",
          formatC(
            p.value,
            format = "f",
            digits = 4
          )
        ),
        TRUE ~ paste0(
          "p = ",
          formatC(
            p.value,
            format = "f",
            digits = 4
          )
        )
      ),
      
      ## Combined annotation
      result_text = paste(
        effect_text,
        p_text,
        sep = ", "
      ),
      
      ## Significance
      significance = dplyr::if_else(
        p.value < 0.05,
        "Significant",
        "Not significant"
      )
    )
  
  
  ## Determine x-axis limits
  if (is.null(x_limits)) {
    
    x_range <- range(
      c(
        plot_dat$conf.low,
        plot_dat$conf.high,
        null_value
      ),
      na.rm = TRUE
    )
    
    xrange <- diff(x_range)
    
    if (xrange == 0) {
      xrange <- ifelse(
        is_binomial,
        0.5,
        1
      )
    }
    
    x_limits <- c(
      x_range[1] - 0.05 * xrange,
      x_range[2] + 0.05 * xrange
    )
    
    x_limits[1] <- min(
      x_limits[1],
      null_value
    )
    
    x_limits[2] <- max(
      x_limits[2],
      null_value
    )
  }
  
  ## Position effect/p-value text
  xrange <- diff(x_limits)
  text_x <- x_limits[2] + 0.02 * xrange
  
  
  ## Plot
  p <- ggplot(
    plot_dat,
    aes(
      x = estimate,
      y = y,
      colour = cohort
    )
  ) +
    
    ## Separator
    annotate(
      "segment",
      x = x_limits[1],
      xend = x_limits[2],
      y = 4.5,
      yend = 4.5,
      linewidth = 0.5,
      colour = "grey50",
      linetype = 3
    )+
    
    ## Confidence intervals
    geom_errorbar(
      aes(
        ymin = y,
        ymax = y,
        xmin = conf.low,
        xmax = conf.high
      ),
      orientation = "y",
      width = 0.25,
      linewidth = 0.6
    )+
    
    ## Point estimates
    geom_point(
      aes(
        shape = significance
      ),
      size = 3,
      stroke = 0.9
    ) +
    
    ## β/OR + p-value
    geom_text(
      aes(
        x = text_x,
        label = result_text,
        fontface = ifelse(
          model == "pcs",
          "bold",
          "plain"
        )
      ),
      hjust = 0,
      size = 3.2,
      show.legend = FALSE
    ) +
    
    ## Null/reference line
    geom_vline(
      xintercept = null_value,
      linetype = 2,
      colour = "grey50"
    ) +
    
    ## Cohort colours
    scale_colour_manual(
      values = stats::setNames(
        c(
          "black",
          "#4D6A7A",
          "#003580"
        ),
        c(
          main_label,
          nonfins_label,
          fins_label
        )
      ),
      breaks = c(
        main_label,
        nonfins_label,
        fins_label
      ),
      labels = c(
        "Complete cohort",
        "Non-Fins",
        "Fins"
      ),
      name = NULL
    ) +
    
    ## Significance symbols
    scale_shape_manual(
      values = c(
        "Significant" = 16,
        "Not significant" = 1
      ),
      breaks = c(
        "Significant",
        "Not significant"
      ),
      labels = c(
        "p < 0.05",
        "p ≥ 0.05"
      ),
      name = NULL
    ) +
    
    ## Y-axis
    scale_y_continuous(
      breaks = 1:6,
      labels = rev(y_labels),
      expand = expansion(
        mult = c(
          0.04,
          0.04
        )
      )
    ) +
    
    ## X-axis
    coord_cartesian(
      xlim = x_limits,
      clip = "off"
    ) +
    
    labs(
      x = xlab,
      y = NULL,
      title = title,
      colour = NULL,
      shape = NULL
    )+
    annotation_custom(
      grid::textGrob(
        term1,
        x = grid::unit(0, "npc"),
        y = grid::unit(-0.08, "npc"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 12)
      )
    ) +
    annotation_custom(
      grid::textGrob(
        term2,
        x = grid::unit(1, "npc"),
        y = grid::unit(-0.08, "npc"),
        just = c("right", "top"),
        gp = grid::gpar(fontsize = 12)
      )
    ) +
    
    ## Theme
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      legend.position = "bottom",
      
      plot.title = element_text(
        size = 14,
        face = "plain"
      ),
      
      plot.margin = margin(
        t = 10,
        r = 150,
        b = 10,
        l = 10
      )
    )+
    guides(
      colour = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    )
  
  
  ## Return
  return(
    list(
      data = plot_dat,
      plot = p
    )
  )
}

###############################################################################-
# 9. SOURCE-SPECIFIC ANALYSES ####
###############################################################################-
## Fit the selected model separately within each source
fit_models_by_source <- function(
    data,
    outcome,
    model_type = "lm",
    level = "base",
    clinical = NULL,
    source_var = "Source"
) {
  
  sources <- unique(data[[source_var]])
  
  models <- lapply(sources, function(src) {
    
    ## Subset to source
    data_src <- data %>%
      dplyr::filter(.data[[source_var]] == src)
    
    ## Build the original formula
    form <- build_formula(
      outcome = outcome,
      level = level,
      random = NULL,
      clinical = clinical
    )
    
    ## Variables included in formula
    vars <- all.vars(form)
    predictors <- setdiff(vars, outcome)
    
    ## Keep only variables actually present in the data
    predictors <- predictors[predictors %in% names(data_src)]
    
    ## Identify invariant categorical variables
    bad_vars <- predictors[
      sapply(
        data_src[predictors],
        function(x) {
          is.factor(x) &&
            dplyr::n_distinct(x[!is.na(x)]) < 2
        }
      )
    ]
    
    ## Remove invariant factors
    if (length(bad_vars) > 0) {
      
      message(
        "Source ", src,
        ": removing invariant factors: ",
        paste(bad_vars, collapse = ", ")
      )
      
      predictors <- setdiff(
        predictors,
        bad_vars
      )
      
      form <- reformulate(
        termlabels = predictors,
        response = outcome
      )
    }
    
    message(
      "Fitting Source ", src,
      " (n = ", nrow(data_src), ")"
    )
    
    ## Fit without random effect because Source is constant
    fit_model_s(
      data = data_src,
      model_type = model_type,
      formula = form
    )
  })
  
  names(models) <- sources
  
  models
}

plot_source_specific <- function(
    data,
    outcome,
    model_type,
    level = "clinical",
    clinical,
    title,
    exponentiate = FALSE,
    measure = "Estimate",
    null_value = 0,
    ncol = 3
) {
  
  ## Fit models
  results_source <- fit_models_by_source(
    data = data,
    outcome = outcome,
    model_type = model_type,
    level = level,
    clinical = clinical
  )
  
  ## Extract forest data
  source_forest_data <- purrr::map_dfr(
    names(results_source),
    function(src) {
      
      if (is.null(results_source[[src]])) {
        return(NULL)
      }
      
      fp <- ForestPlot(
        model = results_source[[src]],
        title = NULL,
        exponentiate = exponentiate,
        measure = measure
      )
      
      fp$data %>%
        mutate(Source = src)
    }
  )
  
  ## Reorder terms
  source_forest_data <- source_forest_data %>%
    mutate(
      term_label = as.character(term_label)
    )
  
  term_order <- source_forest_data %>%
    distinct(term, term_label) %>%
    pull(term_label)
  
  source_forest_data <- source_forest_data %>%
    mutate(
      term_label = factor(
        term_label,
        levels = rev(term_order)
      )
    )
  
  ## Plot
  p <- ggplot(
    source_forest_data,
    aes(
      x = estimate,
      y = term_label,
      colour = sig
    )
  ) +
    geom_point(size = 2.5) +
    
    geom_errorbar(
      aes(
        xmin = conf.low,
        xmax = conf.high
      ),
      orientation = "y",
      width = 0.25
    ) +
    
    geom_vline(
      xintercept = null_value,
      linetype = 2,
      colour = "grey50"
    ) +
    
    facet_wrap(
      ~ Source,
      ncol = ncol
    ) +
    
    scale_colour_manual(
      values = c(
        "Significant" = "firebrick",
        "Not significant" = "black"
      )
    ) +
    
    theme_classic(base_size = 11) +
    
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    ) +
    
    labs(
      x = measure,
      y = NULL,
      colour = NULL,
      title = title
    )
  
  list(
    models = results_source,
    data = source_forest_data,
    plot = p
  )
}



###############################################################################-
# 10. RANDOM-EFFECTS VISUALISATION ####
###############################################################################-
## Plot random slopes by stratum
plot_random_slope <- function(model, random_var, stratum_var = "Source",
                              title = NULL, level2 = NULL,
                              colorpalette = sourceColorPalette) {
  
  re <- ranef(model)[[stratum_var]]
  re_df <- as.data.frame(re)
  re_df[[stratum_var]] <- rownames(re)
  
  fix <- fixef(model)
  
  re_df <- re_df %>%
    mutate(
      abs_intercept = fix["(Intercept)"] + `(Intercept)`,
      abs_slope = .data[[paste0(random_var, level2)]]
    ) %>%
    filter(is.finite(abs_intercept), is.finite(abs_slope))
  
  line_df <- re_df %>%
    expand_grid(!!random_var := c(0, 1)) %>%
    mutate(
      x_numeric = ifelse(!!as.name(random_var) == 0, 0, 1),
      logit = abs_intercept + abs_slope * x_numeric,
      prob = plogis(logit))
  
  ggplot() +
    geom_line(
      data = line_df,
      aes(
        x = x_numeric,
        y = prob,
        group = !!as.name(stratum_var),
        color = !!as.name(stratum_var)
      ),alpha = 0.8) +
    geom_point(
      data = filter(line_df, !!as.name(random_var) == 0),
      aes(
        x = x_numeric,
        y = prob,
        color = !!as.name(stratum_var)),
      size = 2.5) +
    geom_text(
      data = filter(line_df, !!as.name(random_var) == 0),
      aes(
        x = x_numeric,
        y = prob,
        label = !!as.name(stratum_var)),
      nudge_x = 0.05,
      size = 2.5,
      check_overlap = TRUE) +
    scale_x_continuous(
      breaks = c(0, 1),
      labels = c("0", "1")) +
    scale_color_manual(values = colorpalette) +
    labs(
      x = random_var,
      y = "Probability",
      title = title,
      color = "STRATUM") +
    theme_minimal() +
    theme(legend.position = "bottom")}


## Plot random intercepts, slopes, and their relationship
plot_random_intercept_slope <- function(model, random_var,
                                        stratum_var = "Source",
                                        level2 = NULL,
                                        colorpalette = sourceColorPalette,
                                        title_intercept = NULL,
                                        title_slope = NULL) {
  
  re <- ranef(model, condVar = TRUE)
  stratum_re <- re[[stratum_var]]
  post_var <- attr(re[[stratum_var]], "postVar")
  
  re_df <- as.data.frame(stratum_re)
  re_df[[stratum_var]] <- rownames(stratum_re)
  
  re_df$se_intercept <- sqrt(post_var[1, 1, ])
  re_df$se_slope <- sqrt(post_var[2, 2, ])
  
  slope_name <- paste0(random_var, level2)
  
  ## Random intercept plot
  p_intercept <- ggplot(
    re_df,
    aes(
      x = reorder(!!as.name(stratum_var), `(Intercept)`),
      y = `(Intercept)`,
      color = !!as.name(stratum_var)
    )
  ) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(
        ymin = `(Intercept)` - 1.96 * se_intercept,
        ymax = `(Intercept)` + 1.96 * se_intercept
      ),
      width = 0
    ) +
    coord_flip() +
    scale_color_manual(values = colorpalette) +
    labs(
      x = "Stratum",
      y = "Random Intercept",
      title = title_intercept,
      color = "STRATUM"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ## Random slope plot
  p_slope <- ggplot(
    re_df,
    aes(
      x = reorder(
        !!as.name(stratum_var),
        !!as.name(slope_name)
      ),
      y = !!as.name(slope_name),
      color = !!as.name(stratum_var)
    )
  ) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(
        ymin = !!as.name(slope_name) - 1.96 * se_slope,
        ymax = !!as.name(slope_name) + 1.96 * se_slope
      ),
      width = 0
    ) +
    coord_flip() +
    scale_color_manual(values = colorpalette) +
    labs(
      x = "Stratum",
      y = "Random Slope",
      title = title_slope,
      color = "STRATUM"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ## Random intercept vs slope
  p_scatter <- ggplot(
    re_df,
    aes(
      x = `(Intercept)`,
      y = !!as.name(slope_name),
      color = !!as.name(stratum_var),
      label = !!as.name(stratum_var)
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(size = 2.5) +
    geom_text(
      nudge_x = 0.02,
      size = 2.5,
      check_overlap = TRUE
    ) +
    scale_color_manual(values = colorpalette) +
    labs(
      x = "Random Intercept",
      y = "Random Slope",
      title = "Intercept vs Slope",
      color = "STRATUM"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  list(
    p_intercept = p_intercept,
    p_slope = p_slope,
    p_scatter = p_scatter )}

###############################################################################-
# 11. MULTIPLE TESTING ADJUSTMENT ####
###############################################################################-
## Plot forest plots using adjusted p-values
extract_pvalue <- function(model, term,
                           model_names = MODEL_NAMES) {
  
  dat <- broom.mixed::tidy(model)
  
  idx <- which(dat$term == term)
  
  if (length(idx) == 0) {
    stop(
      paste0(
        "Term '",
        term,
        "' not found."
      )
    )
  }
  
  dat$p.value[idx[1]]
}

create_holm_table <- function(
    size_models,
    age_models,
    multiplicity_models,
    location_models,
    term = "PRStotal_Z",
    model_names = MODEL_NAMES
) {

  
  outcomes <- list(
    size = size_models,
    age = age_models,
    multiplicity = multiplicity_models,
    location = location_models
  )
  
  res <- list()
  
  for (m in model_names) {
    
    pvals <- sapply(
      outcomes,
      function(x) extract_pvalue(
        x[[m]],
        term
      )
    )
    
    res[[m]] <- data.frame(
      outcome = names(pvals),
      model = m,
      p.raw = pvals,
      p.adj = p.adjust(
        pvals,
        method = "holm"
      )
    )
  }
  
  dplyr::bind_rows(res)
}