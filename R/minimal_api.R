#' Fit CEClust on a C by lambda grid.
#'
#' `CECfitBoundGrid()` is the high-level entry point for the extended CEClust
#' workflow. It evaluates CEClust across a grid of density bounds `C` and
#' regularisation values `lambda`, keeps the existing stability, saturation,
#' repair, and change-detection logic, and returns a compact object suitable for
#' diagnostics, plotting, or interactive exploration.
#'
#' @param ... Arguments forwarded to the internal grid engine. Common arguments
#'   include `Z`, `lambda_grid`, `C_grid`, `Cquali`, `familyType`, `k0`, `B`,
#'   `Nshots_fresh`, `Nshots_warm`, `Nloop`, `n_cores`, and
#'   `final_descent_from_best`.
#' @param save_results Logical. If `TRUE`, save the returned object as an `.rds`
#'   file. Defaults to `FALSE` in the public wrapper.
#' @param output_dir Directory used when `save_results = TRUE`.
#' @param checkpoint_dir,auto_checkpoint,resume,force_recompute Checkpoint
#'   controls. The public defaults do not write checkpoint files.
#' @return An object of class `"CEC_bound_grid"` containing the original data,
#'   per-`C` lambda-grid fits, a combined `summary` table, and metadata.
#' @export
CECfitBoundGrid <- function(
  ...,
  save_results = FALSE,
  output_dir = file.path(tempdir(), "CEClust"),
  checkpoint_dir = FALSE,
  auto_checkpoint = FALSE,
  resume = FALSE,
  force_recompute = FALSE
) {
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


#' Fit one of the built-in CEClust demonstration workflows.
#'
#' `CECfitPreset()` provides a small set of named entry points for common
#' examples used in documentation and article workflows.
#'
#' @param preset Preset data workflow. One of `"gaussian_1d"`, `"iris"`, or
#'   `"breast_cancer"`.
#' @param ... Additional arguments forwarded to the corresponding grid workflow.
#' @return A `"CEC_bound_grid"` object.
#' @examples
#' \dontrun{
#' iris_grid <- CECfitPreset("iris", k0 = 2, B = 2, Nshots_fresh = 1)
#' CECplotGrid(iris_grid)
#' CECexplore(iris_grid)
#' }
#' @export
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


#' Select stable lambda values from a CEClust grid result.
#'
#' @param x A `"CEC_bound_grid"` object, a single-`C` grid result, or a
#'   `lambda_diag` object returned by [CECfitLambdaGrid()].
#' @param best_parts Optional best-partitions object when `x` is a `lambda_diag`
#'   object.
#' @param C Optional bound value used when `x` contains several `C` values.
#' @param ... Threshold and smoothing arguments forwarded to the stability
#'   selector, such as `stab_algo_threshold`, `rsi_threshold`, `sat_threshold`,
#'   `use_smoothed`, and `min_consecutive`.
#' @return A list describing kept and rejected lambdas plus a diagnostic summary.
#' @export
CECselectStableLambdas <- function(x, best_parts = NULL, C = NULL, ...) {
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


#' Summarise a CEClust C by lambda grid.
#'
#' @param x A `"CEC_bound_grid"` object.
#' @param ... Reserved for future extensions.
#' @return A data frame with one row per `C` and `lambda` pair.
#' @export
CECsummariseGrid <- function(x, ...) {
  cec_grid_summary_data(x)
}


#' Compact a CEClust grid object.
#'
#' `CECcompactGrid()` removes dense intermediate matrices and per-task lambda
#' details that are not needed by the final Shiny explorer, phase diagrams, or
#' selected partitions.
#'
#' @param x A `"CEC_bound_grid"` object.
#' @param keep_lambda_details Logical. Keep the full `lambda_diag$details`
#'   objects. The default `FALSE` is the memory-saving mode.
#' @return A compact `"CEC_bound_grid"` object.
#' @export
CECcompactGrid <- function(x, keep_lambda_details = FALSE) {
  class_x <- class(x)
  out <- cec_grid_compact_grid_object(
    x,
    keep_lambda_details = keep_lambda_details
  )
  class(out) <- unique(c("CEC_bound_grid", "CEC_grid", class_x))
  out
}


#' Extract the selected partition from a CEClust grid.
#'
#' @param x A `"CEC_bound_grid"` object.
#' @param C,lambda Selected grid cell.
#' @return An integer vector of cluster labels.
#' @export
CECextractPartition <- function(x, C, lambda) {
  cec_grid_partition_for_cell(x, C = C, lambda = lambda)
}


#' Plot a CEClust C by lambda grid.
#'
#' @param x A `"CEC_bound_grid"` object.
#' @param ... Arguments forwarded to the grid plot engine. Use
#'   `show_legend = FALSE` to hide the default REO and rejection legend.
#' @return Invisibly returns the plotted summary data.
#' @export
CECplotGrid <- function(x, ...) {
  plot_cec_grid_reo_map(x, ...)
}


#' Plot a selected partition from a CEClust grid.
#'
#' @param x A `"CEC_bound_grid"` object.
#' @param C,lambda Optional selected cell. If omitted, the first retained cell
#'   in the summary is used.
#' @param x_var,y_var Optional variables for multivariate scatter plots. When
#'   omitted, PCA axes are used when available.
#' @param shape_var Optional qualitative variable used for point shapes.
#' @param ... Additional plotting arguments.
#' @return Invisibly returns the selected partition object.
#' @export
CECplotPartition <- function(
  x,
  C = NULL,
  lambda = NULL,
  x_var = NULL,
  y_var = NULL,
  shape_var = "None",
  ...
) {
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


#' Plot a lambda path from a CEClust grid.
#'
#' @param x A `"CEC_bound_grid"` object.
#' @param C Optional bound value. If omitted, the first retained cell is used.
#' @param lambda Optional lambda highlighted on the diagnostic summary path.
#' @param type Plot type. `"summary"` displays the REO path for one `C`;
#'   `"partitions"` displays the one-dimensional partition evolution over
#'   lambda and shades rejected regions.
#' @param lambda_stable_obj Optional output of [CECselectStableLambdas()]. When
#'   omitted and `type = "partitions"`, it is computed from the grid.
#' @param stab_algo_threshold,rsi_threshold,sat_threshold Thresholds used when
#'   computing rejected lambda regions for `type = "partitions"`.
#' @param use_smoothed,min_consecutive Stability-selection options forwarded to
#'   [CECselectStableLambdas()] for `type = "partitions"`.
#' @param ... Additional plotting arguments.
#' @return Invisibly returns the plotted summary table for `type = "summary"`
#'   or the partition-path plotting data for `type = "partitions"`.
#' @export
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


#' Launch the interactive CEClust explorer.
#'
#' `CECexplore()` keeps the Shiny workflow as the main interactive interface for
#' inspecting C by lambda grid results.
#'
#' @param x A `"CEC_bound_grid"` object.
#' @param ... Arguments forwarded to the Shiny explorer.
#' @return The value returned by `shiny::runApp()`.
#' @export
CECexplore <- function(x, ...) {
  launch_cec_grid_shiny(x, ...)
}
