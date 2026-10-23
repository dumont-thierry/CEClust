CECsolver_bound_dispatch <- function(dots,optimization,penalty,optimizer_control,rsi_policy,
  rsi_saturation_max,save_results,output_dir,checkpoint_dir,resume,force_recompute,
  grid_strategy='fixed',adaptive_control=list(),reference_bootstrap=TRUE) {
  if(is.null(dots$Z)) {
    n <- if(is.null(dots$n))100L else CECinteger(dots$n,'n')
    seed <- if(is.null(dots$data_seed))13L else CECinteger(dots$data_seed,'data_seed',0L)
    sim <- CECwith_seed(seed,simulate_gaussian_mixture_1d(n=n,seed=seed))
    dots$Z <- sim$Z
    if(is.null(dots$density))dots$density <- sim$density
  }
  if(grid_strategy=='fixed') {
    if(is.null(dots$C_grid))dots$C_grid <- 1
    if(is.null(dots$lambda_grid))dots$lambda_grid <- seq(.5,2,by=.1)
  }
  family <- dots$familyType
  if(!is.null(family)&&!family%in%c('gaussVector','gaussUniv','gaussAndDiscreteVector'))
    stop('The selected grid engine currently requires Gaussian numeric data.',call.=FALSE)
  if(!is.null(dots$stability_approach)&&dots$stability_approach!='exhaustive')
    stop('The selected grid engine uses exhaustive original-bootstrap diagnostics.',call.=FALSE)
  if(!is.null(dots$Cquali)&&dots$Cquali!=Inf)stop('Qualitative floors are not supported by this grid engine.')
  display <- dots$display_data;labels <- dots$labels;label_name <- if(is.null(dots$label_name))'Label' else dots$label_name
  dataset_name <- dots$dataset_name;density <- dots$density
  # These historical options affect representation/progress only, not fitting.
  for(nm in c('familyType','Cquali','stability_approach','display_data','labels','label_name','dataset_name','density',
    'show_progress','global_progress','compact_results','keep_lambda_details','Nshots_warm','n','data_seed',
    'batch_size','dynamic_task_queue'))dots[[nm]] <- NULL
  if(!is.null(dots$final_descent_from_best)&&isTRUE(dots$final_descent_from_best))
    stop('final_descent_from_best is incompatible with frozen independent bank winners.',call.=FALSE)
  if(!is.null(dots$use_smoothed_keep)&&isTRUE(dots$use_smoothed_keep))
    stop('Smoothed filtering is not supported by the exhaustive grid engine.',call.=FALSE)
  for(nm in c('final_descent_from_best','use_smoothed_keep'))dots[[nm]] <- NULL
  args <- c(dots,list(optimization=optimization,penalty=penalty,optimizer_control=optimizer_control,
    rsi_policy=rsi_policy,rsi_saturation_max=rsi_saturation_max,save_results=FALSE,output_dir=output_dir,
    checkpoint_dir=checkpoint_dir,resume=resume,force_recompute=force_recompute,reference_bootstrap=reference_bootstrap))
  engine <- if(grid_strategy=='adaptive')CECfit_adaptive_grid else CECfit_solver_grid
  if(grid_strategy=='adaptive')args$adaptive_control <- adaptive_control
  unknown <- setdiff(names(args),names(formals(engine)))
  if(length(unknown))stop('Unsupported option(s) for the selected grid engine: ',paste(unknown,collapse=', '),call.=FALSE)
  out <- do.call(engine,args)
  out$meta$grid_strategy <- grid_strategy
  if(!is.null(display)) {
    if(NROW(display)!=NROW(out$Z))stop('display_data must have one row per observation.')
    out$display_data <- as.data.frame(display)
  }
  if(!is.null(labels)) {
    if(length(labels)!=NROW(out$Z))stop('labels must have one value per observation.')
    out$display_data[[label_name]] <- factor(labels)
  }
  out$density <- density;out$dataset_name <- dataset_name;out$meta$dataset <- dataset_name
  if(isTRUE(save_results)) {
    dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    file <- file.path(output_dir,paste0('CEC_',optimization,'_',penalty,'_',substr(out$meta$signature,1,12),'.rds'))
    out$result_file <- file;saveRDS(out,file)
  }
  out
}

CECexhaustive_summaries <- function(original,bootstrap,partition_metric) {
  z <- CECstability(original,bootstrap,partition_metric)
  summarise <- function(values,n_valid) {
    list(mean_sim=if(length(values))mean(values) else NA_real_,
      sd_sim=if(length(values)>1)stats::sd(values) else NA_real_,
      median_sim=if(length(values))stats::median(values) else NA_real_,
      min_sim=if(length(values))min(values) else NA_real_,max_sim=if(length(values))max(values) else NA_real_,
      n_valid=n_valid,medoid_mean_sim=NA_real_)
  }
  list(stab0=summarise(z$algorithm_agreements,length(original)),
    stabB=summarise(z$bootstrap_agreements,length(bootstrap)),ratio=z$RSI)
}

CECsolver_select_stable <- function(x,C=NULL,stab_algo_threshold=NULL,rsi_threshold=NULL,
  sat_threshold=NULL,use_smoothed=FALSE,min_consecutive=1L) {
  if(isTRUE(use_smoothed))stop('Solver diagnostics are filtered without smoothing.',call.=FALSE)
  min_consecutive <- CECinteger(min_consecutive,'min_consecutive')
  if(inherits(x,'CEC_solver_lambda')) {
    if(is.null(x$grid))stop('Use the parent grid for selection after compacting lambda diagnostics.')
    x <- x$grid
  }
  if(is.null(C))C <- cec_grid_best_row(x)$C[1L]
  cc <- x$C_grid[which.min(abs(x$C_grid-C))]
  s <- cec_grid_apply_thresholds(x,stab_algo_threshold,rsi_threshold,sat_threshold)
  s <- s[s$C==cc,,drop=FALSE]
  s <- s[order(s$lambda),,drop=FALSE]
  mask <- s$keep
  if(min_consecutive>1L) {
    runs <- rle(mask);runs$values <- runs$values & runs$lengths>=min_consecutive
    mask <- inverse.rle(runs)
  }
  s$keep <- mask
  list(summary=s,stable_lambda=s$lambda[mask],stable_index=which(mask),stable=mask,
    unstable=!mask,stable_stability=s$keep_stability,unstable_stability=!s$keep_stability,
    exclude_saturation=s$exclude_saturation,
    lambda_min=if(any(mask))min(s$lambda[mask]) else NA_real_,
    stab_algo_threshold=s$stab_algo_threshold[1L],rsi_threshold=s$rsi_threshold[1L],
    sat_threshold=s$sat_threshold[1L],min_consecutive=min_consecutive,use_smoothed=FALSE)
}

# Public help: man/CECrescoreGrid.Rd (maintained by hand).
CECrescoreGrid <- function(x,C_grid=if(inherits(x,'CEC_adaptive_grid'))x$adaptive$initial_C else x$C_grid,
  lambda_grid=if(inherits(x,'CEC_adaptive_grid'))x$adaptive$initial_lambda else x$lambda_grid,
  penalty=x$penalty,rsi_policy=x$meta$rsi_policy,
  grid_strategy=if(inherits(x,'CEC_adaptive_grid'))'adaptive' else 'fixed',
  adaptive_control=if(inherits(x,'CEC_adaptive_grid'))x$adaptive$control else list()) {
  if(!inherits(x,'CEC_solver_grid')||x$optimization!='em_bank')
    stop('CECrescoreGrid requires a retained EM-bank grid.',call.=FALSE)
  if(!identical(digest::digest(unname(as.matrix(x$Z))),x$meta$data_fingerprint))
    stop('The grid data have changed since the candidate banks were fitted.',call.=FALSE)
  grid_strategy <- match.arg(grid_strategy,c('fixed','adaptive'))
  engine <- if(grid_strategy=='adaptive')CECfit_adaptive_grid else CECfit_solver_grid
  args <- list(Z=x$Z,C_grid=C_grid,lambda_grid=lambda_grid,k0=x$meta$k0,B=x$meta$B,
    optimization='em_bank',penalty=penalty,optimizer_control=x$meta$control,
    algo_seed=x$meta$algo_seed,rsi_policy=rsi_policy,rsi_saturation_max=x$meta$rsi_saturation_max,
    stab_algo_threshold=x$meta$stab_algo_threshold,rsi_threshold=x$meta$rsi_threshold,
    sat_threshold=x$meta$sat_threshold,partition_metric=x$meta$partition_metric,
    bootstrap_indices=x$bootstrap_indices,banks=x$banks,
    reference_bootstrap=if(is.null(x$meta$reference_bootstrap))TRUE else x$meta$reference_bootstrap)
  if(grid_strategy=='adaptive')args$adaptive_control <- adaptive_control
  out <- do.call(engine,args);out$meta$grid_strategy <- grid_strategy
  out$display_data <- x$display_data;out$density <- x$density;out$dataset_name <- x$dataset_name
  out$meta$reused_banks <- TRUE;out
}
