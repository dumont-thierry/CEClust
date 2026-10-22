# Public help: man/CECfitBoundGrid.Rd (maintained by hand).
CECfitBoundGrid <- function(
  ...,
  save_results = FALSE,
  output_dir = file.path(tempdir(), "CEClust"),
  checkpoint_dir = FALSE,
  auto_checkpoint = FALSE,
  resume = FALSE,
  force_recompute = FALSE,
  optimization = c("cem", "em_bank"),
  penalty = c("none", "bic"),
  optimizer_control = list(),
  rsi_policy = c("all", "unsaturated", "nontrivial_unsaturated", "saturation_band"),
  rsi_saturation_max = .10,
  grid_strategy = c("fixed", "adaptive"),
  adaptive_control = list()
) {
  rsi_policy_missing <- missing(rsi_policy)
  selected <- !missing(optimization) || !missing(penalty) || length(optimizer_control) > 0 || !missing(rsi_policy) || !missing(grid_strategy) || length(adaptive_control)>0
  optimization <- match.arg(optimization); penalty <- match.arg(penalty); rsi_policy <- match.arg(rsi_policy)
  grid_strategy <- match.arg(grid_strategy)
  if(rsi_policy_missing && grid_strategy=='adaptive' && optimization=='cem')
    rsi_policy <- 'saturation_band'
  if(grid_strategy=='fixed' && length(adaptive_control))stop('Use grid_strategy="adaptive" with adaptive_control.')
  if (selected) {
    if (isTRUE(checkpoint_dir) || (isTRUE(auto_checkpoint) && (is.null(checkpoint_dir) || isFALSE(checkpoint_dir))))
      checkpoint_dir <- file.path(output_dir,"solver_checkpoints")
    return(CECsolver_bound_dispatch(list(...),optimization,penalty,optimizer_control,
      rsi_policy,rsi_saturation_max,save_results,output_dir,checkpoint_dir,resume,force_recompute,
      grid_strategy,adaptive_control))
  }
  out <- run_cec_grid(
    ...,
    save_results = save_results,
    output_dir = output_dir,
    checkpoint_dir = checkpoint_dir,
    auto_checkpoint = auto_checkpoint,
    resume = resume,
    force_recompute = force_recompute
  )
  class(out) <- unique(c("CEC_bound_grid", "CEC_grid", class(out)))
  out
}


# Public help: man/CECfitPreset.Rd (maintained by hand).
CECfitPreset <- function(
  preset = c("gaussian_1d", "iris", "breast_cancer"),
  ...,
  save_results = FALSE,
  output_dir = file.path(tempdir(), "CEClust"),
  checkpoint_dir = FALSE,
  auto_checkpoint = FALSE,
  resume = FALSE,
  force_recompute = FALSE
) {
  preset <- match.arg(preset)
  grid_args <- list(
    ...,
    save_results = save_results,
    output_dir = output_dir,
    checkpoint_dir = checkpoint_dir,
    auto_checkpoint = auto_checkpoint,
    resume = resume,
    force_recompute = force_recompute
  )
  out <- switch(
    preset,
    gaussian_1d = do.call(CECfitBoundGrid, grid_args),
    iris = do.call(run_cec_grid_iris, grid_args),
    breast_cancer = do.call(run_cec_grid_breast_cancer, grid_args)
  )
  class(out) <- unique(c("CEC_bound_grid", "CEC_grid", class(out)))
  out
}


# Public help: man/CECselectStableLambdas.Rd (maintained by hand).
CECselectStableLambdas <- function(x, best_parts = NULL, C = NULL, ...) {
  if (inherits(x, "CEC_solver_grid") || inherits(x, "CEC_solver_lambda"))
    return(CECsolver_select_stable(x, C = C, ...))
  dots <- list(...)
  requested_sat_threshold <- NULL
  if (!is.null(dots$sat_threshold) &&
      length(dots$sat_threshold) == 1L &&
      is.finite(dots$sat_threshold) &&
      dots$sat_threshold <= 0) {
    requested_sat_threshold <- dots$sat_threshold
    dots$sat_threshold <- .Machine$double.eps
  }

  select_keep <- function(lambda_diag, best_parts = NULL) {
    out <- do.call(
      identify_cec_grid_lambda_keep,
      c(
        list(lambda_diag = lambda_diag, best_parts = best_parts),
        dots
      )
    )
    if (!is.null(requested_sat_threshold)) {
      out$sat_threshold <- requested_sat_threshold
    }
    out
  }

  if (inherits(x, "CEC_bound_grid") || !is.null(x$runs_by_C)) {
    result <- if (is.null(C)) {
      best <- cec_grid_best_row(x)
      cec_grid_result_for_C(x, best$C[1L])
    } else {
      cec_grid_result_for_C(x, C)
    }
    return(select_keep(result$lambda_diag, result$best_parts))
  }

  if (!is.null(x$lambda_diag) && !is.null(x$best_parts)) {
    return(select_keep(x$lambda_diag, x$best_parts))
  }

  select_keep(x, best_parts)
}


# Public help: man/CECsummariseGrid.Rd (maintained by hand).
CECsummariseGrid <- function(x, ...) {
  cec_grid_summary_data(x)
}


# Public help: man/CECcompactGrid.Rd (maintained by hand).
CECcompactGrid <- function(x, keep_lambda_details = FALSE) {
  class_x <- class(x)
  out <- cec_grid_compact_grid_object(
    x,
    keep_lambda_details = keep_lambda_details
  )
  class(out) <- unique(c("CEC_bound_grid", "CEC_grid", class_x))
  out
}


# Public help: man/CECextractPartition.Rd (maintained by hand).
CECextractPartition <- function(x, C, lambda) {
  cec_grid_partition_for_cell(x, C = C, lambda = lambda)
}


# Public help: man/CECplotGrid.Rd (maintained by hand).
CECplotGrid <- function(x, ...) {
  plot_cec_grid_reo_map(x, ...)
}


# Public help: man/CECplotPartition.Rd (maintained by hand).
CECplotPartition <- function(
  x,
  C = NULL,
  lambda = NULL,
  x_var = NULL,
  y_var = NULL,
  shape_var = "None",
  ...
) {
  if(inherits(x,'CEC_adaptive_grid') && !is.null(C) && !is.null(lambda)) {
    at <- CECadaptive_cell(x,C,lambda);C <- x$summary$C[at];lambda <- x$summary$lambda[at]
  }
  if (is.null(C) || is.null(lambda)) {
    row <- cec_grid_best_row(x)
    C <- if (is.null(C)) row$C[1L] else C
    lambda <- if (is.null(lambda)) row$lambda[1L] else lambda
  }

  display_data <- cec_grid_display_data(x)
  qvars <- cec_grid_quant_vars(display_data)
  if (length(qvars) <= 1L) {
    return(plot_cec_grid_selected_partition_1d(x, C = C, lambda = lambda, ...))
  }

  choices <- cec_grid_display_choices(display_data)
  if (is.null(x_var)) x_var <- choices[1L]
  if (is.null(y_var)) y_var <- choices[min(2L, length(choices))]
  plot_cec_grid_selected_scatter(
    x,
    C = C,
    lambda = lambda,
    x_var = x_var,
    y_var = y_var,
    shape_var = shape_var,
    ...
  )
}


# Public help: man/CECplotPath.Rd (maintained by hand).
CECplotPath <- function(
  x,
  C = NULL,
  lambda = NULL,
  type = c("summary", "partitions"),
  lambda_stable_obj = NULL,
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.8,
  sat_threshold = 0.1,
  use_smoothed = FALSE,
  min_consecutive = 1,
  ...
) {
  type <- match.arg(type)

  if (is.null(C) || (type == "summary" && is.null(lambda))) {
    row <- cec_grid_best_row(x)
    C <- if (is.null(C)) row$C[1L] else C
    if (type == "summary") {
      lambda <- if (is.null(lambda)) row$lambda[1L] else lambda
    }
  }

  if (type == "partitions") {
    result <- cec_grid_result_for_C(x, C)
    Z <- result$Z
    if (is.data.frame(Z) || is.matrix(Z)) {
      if (NCOL(Z) != 1L) {
        stop(
          "type = 'partitions' is currently implemented for one-dimensional numeric data. ",
          "Use CECplotGrid(), CECplotPartition(), or CECexplore() for multivariate grids.",
          call. = FALSE
        )
      }
      Z <- as.numeric(as.data.frame(Z)[[1L]])
    }
    if (!is.numeric(Z)) {
      stop("type = 'partitions' requires one-dimensional numeric data.", call. = FALSE)
    }

    if (is.null(lambda_stable_obj)) {
      stable_sat_threshold <- sat_threshold
      if (length(stable_sat_threshold) == 1L &&
          is.finite(stable_sat_threshold) &&
          stable_sat_threshold <= 0) {
        stable_sat_threshold <- .Machine$double.eps
      }
      lambda_stable_obj <- CECselectStableLambdas(
        x,
        C = C,
        stab_algo_threshold = stab_algo_threshold,
        rsi_threshold = rsi_threshold,
        sat_threshold = stable_sat_threshold,
        use_smoothed = use_smoothed,
        min_consecutive = min_consecutive
      )
      lambda_stable_obj$sat_threshold <- sat_threshold
    }

    return(plotCECPartitionEvolution1D(
      Z = Z,
      best_partitions_obj = result$best_parts,
      lambda_stable_obj = lambda_stable_obj,
      ...
    ))
  }

  plot_cec_grid_selected_lambda_path_1d(
    x,
    C = C,
    selected_lambda = lambda,
    ...
  )
}


# Public help: man/CECexplore.Rd (maintained by hand).
CECexplore <- function(x, ...) {
  launch_cec_grid_shiny(x, ...)
}

