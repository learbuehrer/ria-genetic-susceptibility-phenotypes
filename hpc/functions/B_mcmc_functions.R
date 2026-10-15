##############################################################################-
# Functions for the MCMC ABN analysis ####
# MCMC ABN helper functions
#
# Purpose:
#   Utility functions used for:
#   - max-parent selection
#   - MCMC network fitting
#   - convergence diagnostics
#   - network visualisation
#   - posterior arc-strength evaluation
#
##############################################################################-

##############################################################################-
# 1. Find optimal no. of parent nodes ####
# this function has been adapted for the HPC environment with parallel processing
##############################################################################-
findMaxParents <- function(
    METHOD      = "mle",
    SCORE       = "bic",
    dat,
    dists,
    banmat      = NULL,
    retainmat  = NULL,
    PRIOR       = 1,
    group.var  = NULL
) {
  
  ## Parallel safety
  if (!foreach::getDoParRegistered()) {
    warning("No parallel backend registered – running sequentially")
    foreach::registerDoSEQ()
  }
  
  cat(paste(
    "\nRun the exact search across incremental parent limits with method:",
    METHOD
  ))
  
  starttime <- Sys.time()

  ## Setup
  novars <- if (is.null(group.var)) ncol(dat) else ncol(dat) - 1
  
  
  startvalue <- if (!is.null(retainmat)) {
    max(1, unname(rowSums(retainmat)[rowSums(retainmat) > 1]))
  } else {
    1
  }
  

  ## Parallel loop over max parents
  net.scores <- foreach::foreach(
    i = startvalue:novars,
    .combine   = 'rbind',
    .packages  = c("abn", "mcmcabn"),
    .inorder   = TRUE,
    .export    = c()
  ) %dopar% {
    
    max.par <- i
    
    mycache <- buildScoreCache(
      data.df      = as.data.frame(dat),
      data.dists   = dists,
      dag.banned   = banmat,
      dag.retained = retainmat,
      max.parents  = max.par,
      method       = METHOD,
      group.var    = group.var
    )
    
    dag.mP <- mostProbable(
      score.cache  = mycache,
      score        = SCORE,
      prior.choice = PRIOR
    )
    
    fabn.mP <- fitAbn(
      object       = dag.mP,
      method       = METHOD,
      prior.choice = PRIOR,
      group.var    = group.var
    )
    
    list(i = list(dag.mP, c(i, SCORE, fabn.mP[[SCORE]])))
  }
  
  endtime <- Sys.time()
  
  cat(paste(
    "\n****************************\n",
    "End MCMC ABN max-parent search, Time used [h]:",
    round(difftime(endtime, starttime, units = "hours"), 4)
  ))
  
  ## Post-processing
  net.scores.dags.prelim <- unlist(net.scores, recursive = FALSE)
  net.scores.dags <- net.scores.dags.prelim[
    seq(1, length(net.scores.dags.prelim), 2)
  ]
  
  net.scores <- net.scores.dags.prelim[
    seq(2, length(net.scores.dags.prelim), 2)
  ]
  
  net.scores <- data.frame(
    npar       = as.integer(unlist(net.scores)[seq(1, length(unlist(net.scores)), 3)]),
    scoretype = unlist(net.scores)[seq(2, length(unlist(net.scores)), 3)],
    scorevalue= as.numeric(unlist(net.scores)[seq(3, length(unlist(net.scores)), 3)])
  )
  
  rownames(net.scores) <- NULL
  
  ## Select optimal max parents
  x <- net.scores[net.scores$scoretype == SCORE, ]
  max.par <- x$npar[which.min(x$scorevalue)][1]
  
  cat(paste("\n****************************\nmax. parents:", max.par))
  
  return(list(
    max.par         = max.par,
    net.scores      = net.scores,
    net.scores.dags = net.scores.dags
  ))
}



##############################################################################-
# 2. Looking at the maxparents: ####
##############################################################################-
# Plot of Score over the number of max parents:
ScorePlotMaxParent <- function(net.scores,FILENAMEbase){
  df <- net.scores %>%
    mutate(
      scoretype = as.factor(scoretype),
      npar = as.factor(npar),
      scorevalue = as.numeric(scorevalue)) %>%
    
    group_by(scoretype) %>%
    reframe(
      scorevalue.norm = (1 - abs(1 - scorevalue / max(scorevalue))) * 100,
      .groups = "keep",
      npar = npar,
      scorevalue = scorevalue)
  
  plt.rel <- ggplot(df, aes(x=npar, color = scoretype, y = scorevalue.norm)) +
    geom_point(aes(shape = scoretype, size = 1), show.legend = F)+
    # facet_wrap(.~dist_type)+
    labs(title =  "Relative network score per no. parent nodes",
         x="number of parent nodes",
         y="network score [%]") +
    scale_x_discrete(limits = unique(df$npar))+
    theme_bw()+
    labs(caption = paste0("Score type = ", net.scores$scoretype[1]))
  
  plt.abs <- ggplot(df, aes(color = scoretype, y = scorevalue, group = scoretype)) +
    geom_point(aes(x=npar, shape = scoretype, size = 1), show.legend = F)+
    # facet_wrap(.~dist_type)+
    
    labs(title = "Absolute network score per no. parent nodes",
         x="number of parent nodes",
         y="network score") +
    scale_x_discrete(limits = unique(df$npar))+
    theme_bw()+
    labs(caption = paste0("Score type = ", net.scores$scoretype[1]))
  
  plt.comb <- cowplot::plot_grid(plt.rel, plt.abs, labels = "AUTO", ncol = 1)
  plt.comb
  
  ggsave(plt.comb, filename = paste0(FILENAMEbase, "ScorePlotMaxParent.png"),
         width = 7,height = 9, dpi = 1200)}


# DAGs for the max parents:
DAGMaxParent <- function(net.scores.dags, net.scores, par1, par2, retainmat = NULL,
                         FILENAMEbase, filename, width = 12, height = 9, dpi = 1200) {
  width_px <- width * dpi
  height_px <- height * dpi
  
  png(filename = paste0(FILENAMEbase, filename, ".png"), width = width_px, height = height_px, res = dpi)
  par(mfrow = c(par1, par2))
  
  if (!is.null(retainmat)) {
    parents <- max(1, unname(max(rowSums(retainmat))))
  } else {
    parents <- 1
  }
  
  for (idx in seq_along(net.scores.dags)) {
    plot(net.scores.dags[[idx]])
    box("figure", col = "gray", lwd = 1)
    title(sub = paste0("DAG with max.par = ", parents, 
                       ", netscore = ", round(net.scores[idx, 3], 3)))
    parents <- parents + 1
  }
  
  dev.off()
}


##############################################################################-
# 3. Fit MCMCM abn() with max.par ####
##############################################################################-
mcmcmABN <- function(
    METHOD= "mle",  SCORE ="bic",              ## method. cache parameters
    dat, dists, banmat=NULL, retainmat=NULL,   ## data cache parameters
    FILENAMEbase,  FILENAME,                   ## computational parameters
    max.par = length(dat)-1,                   ## optimal number of max.par
    MCMC.SEEDS = c("0815", "31415", "12", "14"),## set seeds for evaluation
    MCMC.SCHEME = c(500,0,0),                  ## Meaning: c(RETURN.DAGS, THINNING, BURNIN.LEN)
    MCMC.PRIOR = 2,                            ## 2 = Koivisto Prior, 1 = uniformative prior
    PROB.REV = 0.05,                           ## prob of selecting a new edge reversal
    PROB.MBR = 0.05,                           ## prob of selecting a Markov blanket resampling move
    THRESHOLD = 0.5,                           ## threshold for arc strength
    fileaddendum="",                           ## file name addendum for cross-validation
    group.var=NULL
)

# REV and MBR are efficient in producing high scoring structure but 
# computationally costly compared to classical MCMC jumps.

{
  ## Parallel safety
  ## Do NOT create clusters here; assume the main script has registered
  ## a parallel backend. Fall back to sequential execution if not.
  if (!foreach::getDoParRegistered()) {
    warning("No parallel backend registered – running sequentially")
    foreach::registerDoSEQ()
  }
  
  ## Fit ABN with the given number of max.par
  cat(paste("\nfit ABN with max.par:"))
  starttime <- Sys.time()
  
  mycache.maxpar <- buildScoreCache(
    data.df = as.data.frame(dat),
    data.dists = dists,
    dag.banned = banmat,
    dag.retained = retainmat,
    max.parents = max.par,
    group.var = group.var,
    method = METHOD
  )
  
  dag.maxpar <- mostProbable(
    score.cache = mycache.maxpar,
    score = SCORE
  )
  
  fabn.maxpar <- fitAbn(
    object = dag.maxpar,
    method = METHOD,
    group.var = group.var,
  )
  
  endtime <- Sys.time()
  cat(paste(
    "\n***********************\nEnd ABN with max.par. Time used [h]:",
    round(difftime(endtime, starttime, units = "hours"), 4)
  ))
  

  ## Bootstrapping / MCMC
  starttime <- Sys.time()
  
  cat(paste("\n***********************\nnumber of seeds:", length(MCMC.SEEDS)))
  cat(paste("\nStart MCMC ABN with max.par:", max.par))
  cat(paste(
    "\nMCMC scheme (return dags, thinning, burnin length):",
    paste(MCMC.SCHEME, collapse = ", ")
  ))
  cat(paste("\nStarting Time:", starttime))
  
  mcmc.out.list <- foreach(
    SEED = MCMC.SEEDS,
    .packages = c("abn", "mcmcabn"),
    .inorder = TRUE
  ) %dopar% {
    
    mcmcabn(
      score.cache = mycache.maxpar,
      score = SCORE,
      data.dists = dists,
      max.parents = max.par,
      mcmc.scheme = MCMC.SCHEME,
      seed = SEED,
      verbose = FALSE,
      start.dag = dag.maxpar$dag,
      prob.rev = PROB.REV,
      # REV and MBR are efficient in producing high scoring structure but computationally costly compared to classical MCMC jumps.
      prob.mbr = PROB.MBR,
      prior.choice = MCMC.PRIOR
    )
  }
  
  endtime <- Sys.time()
  cat(paste(
    "\n***********************\n End MCMC ABN. Time used [h]:",
    round(difftime(endtime, starttime, units = "hours"), 4)
  ))
  
  return(list(
    "mcmc.out.list" = mcmc.out.list,
    "dag.maxpar"   = dag.maxpar,
    "fabn.maxpar"  = fabn.maxpar
  ))
}


##############################################################################-
# 4. Postpreprocessing ####
##############################################################################-
postBURNin <- function(mcmc.out.list, burnin.length){
  for (d in 1:length(mcmc.out.list)){
    mcmc.out.list[[d]]$burnin <- burnin.length
    mcmc.out.list[[d]]$dags <- mcmc.out.list[[d]]$dags[,,-(1:burnin.length)]
    mcmc.out.list[[d]]$scores <- mcmc.out.list[[d]]$scores[-(1:burnin.length)]
    mcmc.out.list[[d]]$alpha <- mcmc.out.list[[d]]$alpha[-(1:burnin.length)]
    mcmc.out.list[[d]]$method <- mcmc.out.list[[d]]$method[-(1:burnin.length)]
    mcmc.out.list[[d]]$rejection <- mcmc.out.list[[d]]$rejection[-(1:burnin.length)]
    mcmc.out.list[[d]]$heating <- mcmc.out.list[[d]]$heating[-(1:burnin.length)]
  }
  return(mcmc.out.list)
}


postTHINN <- function(mcmc.out.list, thinningsteps){
  thin <- c(TRUE, rep(FALSE, thinningsteps-1)) # see recycle rules
  mcmc.out.list.thinned <- mcmc.out.list
  
  for (d in 1:length(mcmc.out.list.thinned)){
    if(thinningsteps > mcmc.out.list.thinned[[d]]$iterations){
      stop(paste0("More thinning steps (", thinningsteps, ") than mcmc iterations (", mcmc.out.list.thinned[[d]]$iterations, ")."))
    }
    mcmc.out.list.thinned[[d]]$thinning <- thinningsteps
    mcmc.out.list.thinned[[d]]$dags <- mcmc.out.list.thinned[[d]]$dags[,,thin]
    mcmc.out.list.thinned[[d]]$scores <- mcmc.out.list.thinned[[d]]$scores[thin]
    mcmc.out.list.thinned[[d]]$alpha <- mcmc.out.list.thinned[[d]]$alpha[thin]
    mcmc.out.list.thinned[[d]]$method <- mcmc.out.list.thinned[[d]]$method[thin]
    mcmc.out.list.thinned[[d]]$rejection <- mcmc.out.list.thinned[[d]]$rejection[thin]
  }
  return(mcmc.out.list.thinned)}



##############################################################################-
# 5. MCMC Quality Assessment:  ####
##############################################################################-
## Gelman Diagnostics:
MCMCQualityAssessment <- function(mcmcoutputburnedthinned,mcmcoutput,
                                  FILENAMEbase,filename,  width = 6, height = 6, dpi = 1200,
                                  ymax= 1.5){
  cat("QUALITY ASSESSMENT:\n*********************\nGelman: of higher than 1.1 or 1.2, run chain longer to improve convergence\n")
  
  ## Gelman:
  print(gelman.diag(x = mcmcoutputburnedthinned,autoburnin = F))
  # Calculate the dimensions in pixels
  width_px <- width * dpi
  height_px <- height * dpi
  
  # Open a PNG device
  png(filename = paste0(FILENAMEbase, filename, ".png"), width = width_px, height = height_px, res = dpi)
  gelman.plot(mcmcoutput, autoburnin = F, main="Gelman Plot", ylim=c(1,ymax))
  dev.off()
  
  ## Raftery:
  cat("\n*********************\nRaftery: calculate no. of iterations and no. of burn-ins to satisfy specified condition\n*********************")
  print(raftery.diag((mcmcoutputburnedthinned)))
  
  ##Heidelberg and Welch:
  cat("\n*********************\nHeidelberg and Welch:  test H0: The Markov Chain is from a stationary distribution. If not passed, chain must run longer
      \n*********************")
  for (chain in 1:length(mcmcoutputburnedthinned)){
    cat("\n------------------------\n")
    print(paste("Chain no: ", chain))
    print(heidel.diag(mcmcoutputburnedthinned[[chain]]))}
}


## Trace Plot

TracePlot <- function(modeloutput, mcouthin1, mcoutthin2, mcoutthin3, mcoutthin4,
                      MCMC.SEEDS, dat, dists, METHOD, filename, FILENAMEbase, group.var=NULL){
  fabn <-fitAbn(dag = modeloutput$dag.maxpar$dag,
                data.df = dat,
                data.dists = dists, 
                method = METHOD,
                group.var = group.var) 
  
  max.score <- -fabn$bic
  
  
  dta <- data.frame(mcouthin1[2:4],
                    mcoutthin2[2:4],
                    mcoutthin3[2:4],
                    mcoutthin4[2:4]) 
  dta <- dta[,seq(1,3*length(MCMC.SEEDS),3)]
  names(dta) <- paste0("Run", 1:length(MCMC.SEEDS))
  dta$X <- (1:length(dta$Run1))
  
  dta_long <- dta %>%
    pivot_longer(
      cols = starts_with("Run"),
      names_to = "variable",
      values_to = "scores"
    )
  
  # Compute cumulative max
  dta_long <- dta_long %>%
    arrange(variable, X) %>%
    group_by(variable) %>%
    mutate(cummax = cummax(scores)) %>%
    ungroup()
  
  # Original MCMC scores plot
  original_plot <- ggplot(data = dta_long, aes(x = X, y = scores, color = variable)) +
    geom_line(alpha = 0.8, lwd = 1.1) +
    geom_hline(yintercept = max.score, linetype = "dashed", color = "red") +
    annotate(
      "text",
      x = 25,
      y = max.score,
      label = round(max.score, 2),
      vjust = -0.5,
      color = "red"
    ) +
    labs(x = "DAG index", y = "DAG scores", colour = "MCMC:") +
    ggpubr::theme_pubr() +
    ylim(min(dta_long$scores), max(dta_long$scores))
  
  # Density plot for y-axis
  y_density <- cowplot::axis_canvas(original_plot, axis = "y", coord_flip = TRUE) +
    geom_density(data = dta_long, aes(x = scores, fill = factor(variable)), alpha = 0.5) +
    coord_flip()
  
  # Cumulative max plot
  cummax_plt <- ggplot(data = dta_long, aes(x = X, y = cummax, color = variable)) +
    geom_line(alpha = 0.8, lwd = 1.1) +
    geom_point(aes(color = variable)) +
    geom_hline(yintercept = max.score, linetype = "dashed", color = "red") +
    annotate(
      "text",
      x = 25,
      y = max.score,
      label = round(max.score, 2),
      vjust = -0.5,
      color = "red"
    ) +
    labs(x = "DAG index", y = "Cumulative max score", colour = "MCMC:") +
    ggpubr::theme_pubr() +
    ylim(min(dta_long$scores), max(dta_long$scores))
  
  # Combined plots
  combined_plot <- cowplot::ggdraw(cowplot::insert_yaxis_grob(plot = original_plot, grob = y_density, position = "right"))
  combined_cummax_plot <- cowplot::ggdraw(cowplot::insert_yaxis_grob(plot = original_plot, grob = cummax_plt, position = "right"))
  
  # Save plot to file
  ggsave(
    filename = file.path(FILENAMEbase, paste0(filename, ".png")),
    plot = combined_plot,
    width = 9, height = 7, dpi = 1200
  )
  ggsave(
    filename = file.path(FILENAMEbase, paste0(filename, "_cummax.png")),
    plot = combined_cummax_plot,
    width = 9, height = 7, dpi = 1200
  )
}

##############################################################################-
# 6. Assessment of the DAG: ####
##############################################################################-
plotDAGthresholds <- function(threshold, dat, dists, dag,  
                              FILENAMEbase,filename, width = 12, height = 9, dpi = 1200){
  # Calculate the dimensions in pixels
  width_px <- width * dpi
  height_px <- height * dpi
  
  
  # Open a PNG device
  png(filename = paste0(FILENAMEbase, filename, ".png"), width = width_px, height = height_px, res = dpi)
  
  if((length(threshold))==1){
    par(mfrow=c(1,2))
    dag_threshold <- dag
    dag_threshold[dag_threshold > threshold] <- 1
    dag_threshold[dag_threshold <= threshold] <- 0
    
    edgestren <- dag
    edgestren[edgestren <= threshold] <- 0
    
    cons.dag.plt.edgestrength <- plotAbn(dag = dag_threshold,
                                         data.df = dat,
                                         data.dists = dists,
                                         digits = 2,
                                         edge.strength = edgestren,
                                         plot = T, 
                                         main=threshold)
    title(paste0("Threshold = ",threshold))
    
    # control with bnlearn
    dagBN = bnlearn::empty.graph(names(dat))
    bnlearn::amat(dagBN) = t(dag_threshold)
    bnlearn::graphviz.plot(dagBN,
                           shape = "rectangle",
                           main = "MCMC ABN DAG")
  }
  if((length(threshold)>1)){
    par(mfrow= c(2, ceiling(length(threshold)/2)))
    for(i in 1:length(threshold)){
      dag_threshold <- dag
      dag_threshold[dag_threshold > threshold[i]] <- 1
      dag_threshold[dag_threshold <= threshold[i]] <- 0
      
      edgestren <- round(dag, 2)
      edgestren[edgestren <= threshold[i]] <- 0
      
      cons.dag.plt.edgestrength <- plotAbn(dag = dag_threshold,
                                           data.df = dat,
                                           data.dists = dists,
                                           digits = 2,
                                           edge.strength = edgestren,
                                           plot = T, 
                                           main=threshold[i])
      
      title(sub = paste0("Threshold = ",threshold[i], ", n = ", nrow(dat)))}}
  # Close the device
  dev.off()}

# Arc Strength:
arc.stren.threshold <- function(strengthdag = dag, method = "l1",
                                width = 12, height = 8, dpi = 1200) {
  # Handle empty input
  if (length(strengthdag) == 0)
    return(0)
  
  # Compute ECDF and knots
  e <- ecdf(strengthdag)
  u <- knots(e)
  
  # Define norm function based on method
  norm <- switch(method,
                 "l1" = function(p) {
                   sum(diff(unique(c(0, u, 1))) * abs(e(unique(c(0, u[u < 1]))) - p))
                 },
                 "l2" = function(p) {
                   sum(diff(unique(c(0, u, 1))) * (e(unique(c(0, u[u < 1]))) - p)^2)
                 },
                 stop("Unsupported method. Use 'l1' or 'l2'."))
  
  # Optimize norm
  p0 <- optimize(f = norm, interval = c(0, 1))$minimum
  
  # Check boundaries
  if (norm(1) < norm(p0)) p0 <- 1
  if (norm(0) < norm(p0)) p0 <- 0
  
  # Compute threshold
  arc.stren.sign.threshold <- quantile(strengthdag, p0, type = 1, names = FALSE)
  
  return(arc.stren.sign.threshold)}


# EDCF-plot
ecdf_plot_ggplot <- function(strengthdag = dag,
                             arc.stren.sign.threshold = arc.stren.sign.threshold,
                             filename, FILENAMEbase,
                             width = 12, height = 8, dpi = 1200) {
  # Handle empty input
  if (length(strengthdag) == 0)
    return(0)
  
  # Compute ECDF
  e <- ecdf(strengthdag)
  x_vals <- sort(unique(strengthdag))
  y_vals <- e(x_vals)
  df <- data.frame(x = x_vals, y = y_vals)
  
  # Add (0, 0) manually
  df <- data.frame(x = c(0, x_vals), y = c(0, y_vals))
  
  
  # Create ggplot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_hline(yintercept = 1, color = "black", linewidth = 0.5, lty = 2) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5, lty = 2) +
    geom_step(color = "black") +
    geom_point()+
    geom_vline(xintercept = arc.stren.sign.threshold, linetype = "dashed", color = "darkred") +
    geom_vline(xintercept = 0.5, color = "gray", linewidth = 0.5) +
    labs(title = paste0("ECDF(DAG), threshold = ", arc.stren.sign.threshold),
         x = "Arc Strength", y = "Fn(x)") +
    theme_classic() +
    xlim(c(0, 1)) +
    ylim(c(0, 1))
  
  
  # Save plot
  ggsave(filename = paste0(FILENAMEbase, filename, ".png"),
         plot = p, width = width, height = height, dpi = dpi)
  
  return(invisible(p))}



# threshold map:
ThresholdHeatmap <- function(
    dag,
    threshold,
    banmat,
    filename,
    FILENAMEbase,
    width = 12,
    height = 7,
    dpi = 1200) {
  ## 1. Dictionary for readable variable labels
  
  
  term_dictionary <- c(
    "SexMale" = "Male",
    "Age.at.time.of.SAH_Z" = "Age at Rupture",
    "allele_count_Z" = "Allele Count",
    "PRStotal_Z" = "PRS",
    "PRS1_Z" = "PRS1",
    "PRS2_Z" = "PRS2",
    "PRS3_Z" = "PRS3",
    "PRS4_Z" = "PRS4",
    "PRS5_Z" = "PRS5",
    "Positive.famillial.history" = "Family History Awareness",
    "Smoking_Current_Former_No" = "Smoking Status",
    "Hypertension" = "Hypertension Awareness",
    "IA_size.at.time.of.SAH.log_Z" = "IA Size",
    "IAlocation_cow" =      "IA Location",
    "genotyping_methodIl" =      "Genotyping Method",
    "PRS_LOSO_Z" = "PRS LOSO",
    "PRStotal_Z" = "PRS",
    "Multiple.IAs" = "Multiple IAs",
    "PC1_AVG_Z" = "PC1",
    "PC2_AVG_Z" = "PC2",
    "PC3_AVG_Z" = "PC3",
    "PC4_AVG_Z" = "PC4",
    "PC5_AVG_Z" = "PC5" )
  
  
  
  ## 2. Ensure node names are consistent
  nodes <- rownames(banmat)
  
  if (is.null(nodes)) {
    stop("banmat must have row names.")
  }
  
  rownames(dag) <- nodes
  colnames(dag) <- nodes
  
  
  
  ## 3. Create readable labels
  readable_nodes <- nodes
  matched <- nodes %in% names(term_dictionary)
  readable_nodes[matched] <- term_dictionary[
    nodes[matched]]
  names(readable_nodes) <- nodes
  
  
  
  ## 4. Convert DAG matrix to long format
  df <- as.data.frame(dag)
  df$from <- rownames(df)
  
  df <- df %>%
    tidyr::pivot_longer(
      cols = -from,
      names_to = "to",
      values_to = "value"
    )
  
  df$from <- factor(
    df$from,
    levels = nodes)
  
  df$to <- factor(
    df$to,
    levels = nodes)
  
  df$diagonal <- df$from == df$to
  
  
  
  ## 5. Common fill scale
  fill_limits <- c(
    0,
    max(df$value, na.rm = TRUE))
  
  
  
  ## 6. Publication theme
  publication_theme <- ggplot2::theme_minimal(
    base_size = 11
  ) +
    ggplot2::theme(
      
      panel.grid = element_blank(),
      
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        colour = "black"
      ),
      
      axis.text.y = element_text(
        colour = "black"
      ),
      
      axis.title = element_text(
        colour = "black",
        size = 11
      ),
      
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        size = 10
      ),
      
      legend.text = element_text(
        size = 9
      ),
      
      legend.key.width = unit(
        1.5,
        "cm"
      ),
      
      legend.key.height = unit(
        0.4,
        "cm"
      ),
      
      plot.margin = margin(
        t = 8,
        r = 8,
        b = 4,
        l = 8
      )
    )
  
  
  
  ## 7. Plot A: raw posterior arc strengths
  a <- ggplot(
    df,
    aes(
      x = to,
      y = from,
      fill = value
    )
  ) +
    
    geom_tile(
      colour = "grey70",
      linewidth = 0.2
    ) +
    
    # Diagonal: white + black X
    geom_tile(
      data = df[df$diagonal, ],
      aes(
        x = to,
        y = from
      ),
      fill = "white",
      colour = "grey70",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    
    geom_point(
      data = df[df$diagonal, ],
      aes(
        x = to,
        y = from
      ),
      shape = 4,
      size = 2.5,
      stroke = 0.8,
      colour = "black",
      inherit.aes = FALSE
    ) +
    
    scale_fill_gradient(
      name = "Arc strength",
      low = "white",
      high = "black",
      limits = fill_limits,
      breaks = scales::pretty_breaks(n = 5),
      expand = c(0, 0)
    ) +
    
    scale_x_discrete(
      labels = readable_nodes,
      position = "top",
      expand = expansion(mult = c(0, 0))
    ) +
    
    scale_y_discrete(
      labels = readable_nodes,
      limits = rev(nodes),
      expand = expansion(mult = c(0, 0))
    ) +
    
    labs(
      title = "A  Posterior arc strengths",
      x = "Parent",
      y = "Child"
    ) +
    
    publication_theme
  
  
  
  ## 8. Thresholded data
  df_threshold <- df
  
  df_threshold$value[
    df_threshold$value < threshold
  ] <- 0
  
  
  
  ## 9. Plot B: thresholded posterior arc strengths
  b <- ggplot(
    df_threshold,
    aes(
      x = to,
      y = from,
      fill = value
    )
  ) +
    
    geom_tile(
      colour = "grey80",
      linewidth = 0.2
    ) +
    
    # Diagonal: white + black X
    geom_tile(
      data = df_threshold[df_threshold$diagonal, ],
      aes(
        x = to,
        y = from
      ),
      fill = "white",
      colour = "grey80",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    
    geom_point(
      data = df_threshold[df_threshold$diagonal, ],
      aes(
        x = to,
        y = from
      ),
      shape = 4,
      size = 2.5,
      stroke = 0.8,
      colour = "black",
      inherit.aes = FALSE
    ) +
    
    scale_fill_gradient(
      name = "Arc strength",
      low = "white",
      high = "black",
      limits = fill_limits,
      breaks = scales::pretty_breaks(n = 5),
      expand = c(0, 0)
    ) +
    
    scale_x_discrete(
      labels = readable_nodes,
      position = "top",
      expand = expansion(mult = c(0, 0))
    ) +
    
    scale_y_discrete(
      labels = readable_nodes,
      limits = rev(nodes),
      expand = expansion(mult = c(0, 0))
    ) +
    
    labs(
      title = paste0(
        "B  Thresholded arc strengths (≥ ",
        threshold,
        ")"
      ),
      x = "Parent",
      y = "Child"
    ) +
    
    publication_theme
  
  
  
  ## 10. Combine plots
  combined_plot <- ggpubr::ggarrange(
    a,
    b,
    nrow = 1,
    ncol = 2,
    common.legend = TRUE,
    legend = "bottom",
    widths = c(1, 1)
  )
  
  
  
  ## 11. Save figure
  ggsave(
    filename = file.path(
      FILENAMEbase,
      paste0(filename, ".png")
    ),
    plot = combined_plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

ThresholdHeatmapSingle <- function(
    dag,
    threshold,
    banmat,
    filename,
    FILENAMEbase,
    width = 8,
    height = 8,
    dpi = 1200
) {
  
  
  # Node names
  
  
  nodes <- rownames(banmat)
  
  if (is.null(nodes)) {
    stop("banmat must have row names defining the node names.")
  }
  
  if (nrow(dag) != length(nodes) || ncol(dag) != length(nodes)) {
    stop("Dimensions of dag and banmat do not match.")
  }
  
  rownames(dag) <- nodes
  colnames(dag) <- nodes
  
  
  
  # Readable variable labels
  
  term_dictionary <- c(
    "SexMale" = "Male",
    "Age.at.time.of.SAH_Z" = "Age at Rupture",
    "allele_count_Z" = "Allele Count",
    "PRStotal_Z" = "PRS",
    "PRS1_Z" = "PRS1",
    "PRS2_Z" = "PRS2",
    "PRS3_Z" = "PRS3",
    "PRS4_Z" = "PRS4",
    "PRS5_Z" = "PRS5",
    "Positive.famillial.history" = "Family History Awareness",
    "Smoking_Current_Former_No" = "Smoking Status",
    "Hypertension" = "Hypertension Awareness",
    "IA_size.at.time.of.SAH.log_Z" = "IA Size",
    "IAlocation_cow" =      "IA Location",
    "genotyping_methodIl" =      "Genotyping Method",
    "PRS_LOSO_Z" = "PRS LOSO",
    "PRStotal_Z" = "PRS",
    "Multiple.IAs" = "Multiple IAs",
    "PC1_AVG_Z" = "PC1",
    "PC2_AVG_Z" = "PC2",
    "PC3_AVG_Z" = "PC3",
    "PC4_AVG_Z" = "PC4",
    "PC5_AVG_Z" = "PC5" )
  
  readable_nodes <- nodes
  
  matched <- nodes %in% names(term_dictionary)
  
  readable_nodes[matched] <-
    unname(term_dictionary[nodes[matched]])
  
  names(readable_nodes) <- nodes
  
  
  
  # Convert matrix to long format
  
  
  df <- as.data.frame(dag)
  
  df$from <- rownames(df)
  
  df <- df |>
    tidyr::pivot_longer(
      cols = -from,
      names_to = "to",
      values_to = "value"
    )
  
  df$from <- factor(
    df$from,
    levels = nodes
  )
  
  df$to <- factor(
    df$to,
    levels = nodes
  )
  
  df$diagonal <- df$from == df$to
  
  
  
  # Keep diagonal separate
  
  
  df_plot <- df |>
    dplyr::filter(!diagonal)
  
  df_diag <- df |>
    dplyr::filter(diagonal)
  
  
  
  # Colour scale
  #
  # Below threshold: white -> grey
  # Above threshold: light blue -> dark blue
  
  
  # Colour scale
  #
  # <= threshold: white -> grey
  # > threshold: light teal -> dark teal
  
  max_value <- max(df_plot$value, na.rm = TRUE)
  
  if (max_value <= 0) {
    max_value <- 1
  }
  
  # Keep threshold within the plotted range
  threshold_plot <- min(
    max(threshold, 0),
    max_value
  )
  
  # Position of the threshold on the 0-1 scale
  threshold_position <- threshold_plot / max_value
  
  # Very small separation above the threshold so that
  # values > threshold enter the coloured range.
  epsilon <- min(
    0.001,
    (1 - threshold_position) / 2
  )
  
  if (epsilon <= 0) {
    epsilon <- 0.0001
  }
  
  colour_values <- c(
    0,
    threshold_position,
    min(1, threshold_position + epsilon),
    1
  )
  
  colour_palette <- c(
    "#FFFFFF",
    "#737373",
    "#9CC8BF",
    "#2F6F65"
  )
  
  colour_palette <- colour_palette[
    seq_along(colour_values)
  ]
  
  
  # Plot
  
  
  p <- ggplot(
    df_plot,
    aes(
      x = to,
      y = from,
      fill = value
    )
  ) +
    
    geom_tile(
      colour = "grey85",
      linewidth = 0.2
    ) +
    
    # Diagonal: white cells
    geom_tile(
      data = df_diag,
      aes(
        x = to,
        y = from
      ),
      fill = "white",
      colour = "grey85",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    
    # Diagonal: cross
    geom_text(
      data = df_diag,
      aes(
        x = to,
        y = from
      ),
      label = "×",
      size = 5,
      fontface = "plain",
      colour = "black",
      inherit.aes = FALSE
    ) +
    
    scale_fill_gradientn(
      name = "Posterior arc strength",
      colours = colour_palette,
      values = colour_values,
      limits = c(0, max_value),
      breaks = scales::pretty_breaks(n = 5),
      labels = scales::label_number(
        accuracy = 0.01
      ),
      expand = c(0, 0)
    ) +
    
    scale_x_discrete(
      labels = readable_nodes,
      position = "top",
      expand = expansion(mult = c(0, 0))
    ) +
    
    scale_y_discrete(
      labels = readable_nodes,
      limits = rev(nodes),
      expand = expansion(mult = c(0, 0))
    ) +
    
    labs(
      x = "Parent",
      y = "Child",
      title = "Posterior arc strengths",
      subtitle = paste0(
        "Threshold = ",
        formatC(
          threshold,
          format = "f",
          digits = 3
        )
      )
    ) +
    
    theme_minimal(
      base_size = 11
    ) +
    
    theme(
      panel.grid = element_blank(),
      
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 0,
        colour = "black",
        size = 9
      ),
      
      axis.text.y = element_text(
        colour = "black",
        size = 9
      ),
      
      axis.title.x = element_text(
        colour = "black",
        size = 11,
        margin = margin(t = 8)
      ),
      
      axis.title.y = element_text(
        colour = "black",
        size = 11,
        margin = margin(r = 8)
      ),
      
      plot.title = element_text(
        size = 13,
        face = "bold",
        colour = "black"
      ),
      
      plot.subtitle = element_text(
        size = 10,
        colour = "grey30"
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        size = 10
      ),
      
      legend.text = element_text(
        size = 9
      ),
      
      legend.key.width = unit(
        2.0,
        "cm"
      ),
      
      legend.key.height = unit(
        0.45,
        "cm"
      ),
      
      plot.margin = margin(
        10, 10, 10, 10
      )
    )
  
  
  
  # Save
  ggsave(
    filename = file.path(
      FILENAMEbase,
      paste0(filename, ".png")
    ),
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

# adjacency matrix
AdjacencyMatrix <- function(banmat, dag, threshold,
                            filename, FILENAMEbase, width=12, height=8, dpi = 1200){
  
  # Ban matrix long format
  df <- as.data.frame(banmat)
  df$from <- rownames(df)
  df <- df %>%
    tidyr::pivot_longer(
      cols = -from,
      names_to = "to",
      values_to = "value"
    )
  
  a <- ggplot(df, aes(x = to, y = from, fill = value)) +
    geom_tile(show.legend = F) +
    scale_fill_gradient(low = "palegreen3", high = "white") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(caption= "Banmatrix") +
    scale_y_discrete(limits = rev(unique(df$from))) +
    scale_x_discrete(expand = expansion(mult = c(0,0)), position = "top") +
    xlab("parent") +
    ylab("child")
  
  # Thresholded DAG long format
  dag_threshold <- dag
  dag_threshold[dag_threshold > threshold] <- 1
  dag_threshold[dag_threshold <= threshold] <- 0
  rownames(dag_threshold) <- colnames(dag_threshold) <- rownames(banmat)
  
  dfdag <- as.data.frame(dag_threshold)
  dfdag$from <- rownames(dfdag)
  dfdag <- dfdag %>%
    tidyr::pivot_longer(
      cols = -from,
      names_to = "to",
      values_to = "value"
    )
  
  b <- ggplot(dfdag, aes(x = to, y = from, fill = value)) +
    geom_tile(show.legend = F) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(caption= "DAG fitted") +
    scale_y_discrete(limits = rev(unique(df$from))) +
    scale_x_discrete(expand = expansion(mult = c(0,0)), position = "top") +
    xlab("parent") +
    ylab("child")
  
  combined <- ggpubr::ggarrange(a, b, nrow=1, ncol=2, labels = c("A", "B"))
  
  ggsave(plot = combined,
         filename = file.path(FILENAMEbase, paste0(filename, ".png")),
         width = width, height = height, dpi = dpi)
}

AdjacencyMatrixComparision <- function(banmat, dag_train, dag_test, 
                                       filename, FILENAMEbase, width=12, height=8, dpi = 1200){
  
  # Ban matrix long format
  df <- as.data.frame(banmat)
  df$from <- rownames(df)
  df <- df %>%
    tidyr::pivot_longer(
      cols = -from,
      names_to = "to",
      values_to = "value"
    )
  
  # DAG train long format
  df_train <- as.data.frame(dag_train)
  df_train$from <- rownames(df_train)
  df_train <- df_train %>%
    tidyr::pivot_longer(
      cols = -from,
      names_to = "to",
      values_to = "value"
    )
  
  a <- ggplot(df_train, aes(x = to, y = from, fill = value)) +
    geom_tile(show.legend = F) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(caption= "DAG train") +
    scale_y_discrete(limits = rev(unique(df$from))) +
    scale_x_discrete(expand = expansion(mult = c(0,0)), position = "top") +
    xlab("parent") +
    ylab("child")
  
  # DAG test long format
  df_test <- as.data.frame(dag_test)
  df_test$from <- rownames(df_test)
  df_test <- df_test %>%
    tidyr::pivot_longer(
      cols = -from,
      names_to = "to",
      values_to = "value"
    )
  
  b <- ggplot(df_test, aes(x = to, y = from, fill = value)) +
    geom_tile(show.legend = F) +
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(caption= "DAG test") +
    scale_y_discrete(limits = rev(unique(df$from))) +
    scale_x_discrete(expand = expansion(mult = c(0,0)), position = "top") +
    xlab("parent") +
    ylab("child")
  
  combined <- ggpubr::ggarrange(a, b, nrow=1, ncol=2, labels = c("A", "B"))
  
  ggsave(plot = combined,
         filename = file.path(FILENAMEbase, paste0(filename, ".png")),
         width = width, height = height, dpi = dpi)
}



DAGsabn <- function(seeds = c("0815", "1234"),dat, dists, banmat= NULL, retainmat= NULL,
                    max.parents=length(dat)-1, group.var = NULL,
                    filename, FILENAMEbase, width=12, height=8, dpi = 1200, prior=1, par1, par2, SCORE="bic"){
  # Calculate the dimensions in pixels
  width_px <- width * dpi
  height_px <- height * dpi
  
  # Open a PNG device
  png(filename = paste0(FILENAMEbase, filename, ".png"), width = width_px, height = height_px, res = dpi)
  par(mfrow= c(par1,par2))
  
  
  for (i in 1:length(seeds)){
    chacheMLE <- buildScoreCache(data.df = dat,
                                 data.dists = dists,
                                 method = "mle",
                                 dag.banned = banmat,
                                 dag.retained = retainmat,
                                 centre = TRUE,
                                 max.parents = max.parents,
                                 group.var = group.var)
    
    dag <- mostProbable(score.cache = chacheMLE, prior.choice = prior,
                        score=SCORE)
    plot(dag)
    box("figure", col="gray", lwd = 1)
    title(sub = paste0("seed = ", seeds[i]))}
  
  
  dev.off()}


abnMCMCabnComparision <- function(dat, dists, dag,  threshold,prior=1, banmat= NULL,
                                  retainmat = NULL,
                                  FILENAMEbase,filename, width = 12, height = 9, dpi = 1200,
                                  max.parents, seed= "0815",
                                  group.var=NULL){
  # Calculate the dimensions in pixels
  width_px <- width * dpi
  height_px <- height * dpi
  
  # Open a PNG device
  png(filename = paste0(FILENAMEbase, filename, ".png"), width = width_px, height = height_px, res = dpi)
  par(mfrow=c(1,2))
  dag_threshold <- dag
  dag_threshold[dag_threshold > threshold] <- 1
  dag_threshold[dag_threshold <= threshold] <- 0
  edgestren <- dag
  edgestren[edgestren <= threshold] <- 0
  cons.dag.plt.edgestrength <- plotAbn(dag = dag_threshold,
                                       data.df = dat,
                                       data.dists = dists,
                                       digits = 2,
                                       edge.strength = edgestren,
                                       plot = T, 
                                       main=threshold)
  title(sub = paste0("mcmcabn(), threshold = ",threshold, ", n = ", nrow(dat)))
  set.seed(seed)
  chacheMLE1 <- buildScoreCache(data.df = dat,
                                data.dists = dists,
                                method = "mle",
                                dag.banned = banmat,
                                dag.retained = retainmat,
                                centre = TRUE,
                                max.parents = max.parents,
                                group.var = group.var)
  
  dagabn <- mostProbable(score.cache = chacheMLE1, prior.choice = prior)
  plot(dagabn)
  title(sub = paste0("abn(), seed = ", seed, ", n = ", nrow(dat)))
  
  # Close the device
  dev.off()}