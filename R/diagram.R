# Common entry point. Numerical solvers remain in their own source files.
CECinitialGrid <- function(X, powers = 1:10) {
  X <- CECnumeric_data(X)
  if (!is.numeric(powers) || length(powers) < 2L || any(!is.finite(powers)) ||
      any(powers != floor(powers)) || any(diff(powers) <= 0) ||
      any(diff(powers) != diff(powers)[1L]))
    stop("powers must be at least two increasing, equally spaced integers.", call. = FALSE)
  stats <- ca_stats(X, rep(1L, nrow(X)))
  if (stats$rank < ncol(X))
    stop("Cannot compute Cmin: the full-sample covariance is singular. Supply C_grid explicitly or remove redundant variables.", call. = FALSE)
  log_peak <- -.5 * (ncol(X) * log(2 * pi) + stats$logdet)
  Cmin <- exp(log_peak)
  C_grid <- exp(log_peak + powers * log(2))
  if (!is.finite(Cmin) || Cmin <= 0 || any(!is.finite(C_grid) | C_grid <= 0))
    stop("The density peak or initial C grid is outside numerical range.", call. = FALSE)
  list(Cmin = Cmin, C_grid = C_grid, powers = powers)
}

CECdiagramControl <- function(method = c("adaptive", "bands"), ...) {
  method <- match.arg(method)
  supplied <- list(...)
  if (length(supplied) && (is.null(names(supplied)) || any(!nzchar(names(supplied))) ||
      anyDuplicated(names(supplied))))
    stop("Control arguments must have unique names.", call. = FALSE)
  if (method == "bands") {
    defaults <- list(extra_starts = 2L, max_anchors = 10L, max_order = 40L,
      C_points = 11L, max_C_points = 41L, refine_rounds = 2L,
      projection_rounds = 2L, max_projection_queries = 200L)
    return(do.call(bands_control, utils::modifyList(defaults, supplied)))
  }
  defaults <- list(optimization = "cem", r0 = 20L, Nshots_fresh = 5L,
    Nloop = 1000L, rsi_policy = "saturation_band",
    adaptive_control = list(extend_left = FALSE))
  common <- c("Z", "C_grid", "lambda_grid", "k0", "B", "algo_seed", "n_cores",
    "penalty", "bootstrap_indices", "checkpoint_dir", "resume", "save_results", "output_dir")
  allowed <- setdiff(names(formals(CECfit_adaptive_grid)), common)
  unknown <- setdiff(names(supplied), allowed)
  if (length(unknown)) stop("Unknown adaptive control: ", paste(unknown, collapse = ", "), call. = FALSE)
  out <- utils::modifyList(defaults, supplied)
  if (!is.list(out$adaptive_control)) stop("adaptive_control must be a named list.", call. = FALSE)
  out$adaptive_control <- CECadaptive_control(out$adaptive_control,
    if (is.null(out$sat_threshold)) 0 else out$sat_threshold)
  out
}

CECdiagram <- function(X, method = c("adaptive", "bands"), L = 20L, B = 20L,
    workers = 1L, seed = 1L, penalty = c("none", "bic"), control = list(),
    C_grid = NULL, lambda_grid = NULL, bootstrap_indices = NULL,
    directory = tempfile("CEC-diagram-"), resume = FALSE) {
  method <- match.arg(method); penalty <- match.arg(penalty)
  if (!is.list(control)) stop("control must be a named list.", call. = FALSE)
  if (!is.character(directory) || length(directory) != 1L || is.na(directory) || !nzchar(directory))
    stop("directory must be a nonempty path.", call. = FALSE)
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume))
    stop("resume must be TRUE or FALSE.", call. = FALSE)
  L <- CECinteger(L, "L", 2L); B <- CECinteger(B, "B", 1L)
  workers <- CECinteger(workers, "workers"); seed <- CECinteger(seed, "seed", 0L)
  settings <- do.call(CECdiagramControl, c(list(method = method), control))
  if (method == "bands") {
    if (!is.null(C_grid) || !is.null(lambda_grid))
      stop("The bands engine constructs its own diagram. Use C_points/max_C_points in control, not C_grid or lambda_grid.", call. = FALSE)
    fit <- cec_bands(X, L = L, B = B, workers = workers, seed = seed,
      penalty = penalty, control = settings, bootstrap_indices = bootstrap_indices,
      directory = directory, resume = resume)
    initialization <- NULL
  } else {
    X <- CECnumeric_data(X)
    initialization <- if (is.null(C_grid)) CECinitialGrid(X) else NULL
    if (!is.null(initialization)) C_grid <- initialization$C_grid
    # Never silently overwrite a previous run unless resumption was requested.
    if (dir.exists(directory) && !resume && length(list.files(directory, all.files = TRUE, no.. = TRUE)))
      stop("Archive already exists; choose a new directory or resume=TRUE.", call. = FALSE)
    fit <- do.call(CECfit_adaptive_grid, c(list(Z = X, C_grid = C_grid,
      lambda_grid = lambda_grid, k0 = L, B = B, n_cores = workers,
      algo_seed = seed, penalty = penalty, bootstrap_indices = bootstrap_indices,
      checkpoint_dir = directory, resume = resume, save_results = FALSE), settings))
  }
  attr(fit, "CECdiagram") <- list(method = method, version = as.character(utils::packageVersion("CEClust")),
    initialization = initialization, directory = normalizePath(directory, winslash = "/", mustWork = TRUE))
  CECwrite_rds_atomic(fit, file.path(directory, "diagram.rds"))
  fit
}
