bands_integer <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != floor(x) ||
      x < minimum || x > .Machine$integer.max)
    stop(name, " must be an integer >= ", minimum, ".")
  as.integer(x)
}

bands_control <- function(orders = 1:20, starts = 5L, extra_starts = 5L,
    max_order = 320L, max_iter = 1000L, tolerance = 1e-8,
    max_anchors = 40L, C_points = 41L, max_C_points = 257L,
    refine_rounds = 2L, projection_rounds = 2L, max_projection_queries = 2000L,
    cache_mb = 64, diag_mode = c("components", "all"),
    algo_min = .85, rsi_min = .95, saturation_max = 0,
    measurement_saturation = .10, branches = c("free", "bounded")) {
  diag_mode <- match.arg(diag_mode)
  if (!length(orders) || any(!is.finite(orders) | orders < 1 | orders != floor(orders)))
    stop("orders must be positive integers.")
  branches <- unique(match.arg(branches, c("free", "bounded"), several.ok = TRUE))
  values <- list(starts = starts, extra_starts = extra_starts, max_order = max_order,
    max_iter = max_iter, max_anchors = max_anchors, C_points = C_points,
    max_C_points = max_C_points, refine_rounds = refine_rounds,
    projection_rounds = projection_rounds, max_projection_queries = max_projection_queries)
  for (name in names(values)) values[[name]] <- bands_integer(values[[name]], name,
    if (name %in% c("extra_starts", "refine_rounds", "projection_rounds")) 0L else 1L)
  if (values$C_points < 3L || values$max_C_points < values$C_points) stop("Invalid C grid budgets.")
  if (max(orders) > values$max_order) stop("orders exceed max_order.")
  thresholds <- c(algo_min, rsi_min, saturation_max, measurement_saturation)
  if (length(thresholds) != 4L || any(!is.finite(thresholds) | thresholds < 0 | thresholds > 1) ||
      saturation_max > measurement_saturation) stop("Invalid diagnostic thresholds.")
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
      length(cache_mb) != 1L || !is.finite(cache_mb) || cache_mb <= 0) stop("Invalid numerical control.")
  c(list(orders = sort(unique(as.integer(orders))), tolerance = tolerance,
    cache_mb = cache_mb, diag_mode = diag_mode, algo_min = algo_min,
    rsi_min = rsi_min, saturation_max = saturation_max,
    measurement_saturation = measurement_saturation, branches = branches), values)
}

bands_runtime_cache <- new.env(parent = emptyenv())

bands_runtime <- function() {
  if (exists("value", bands_runtime_cache, inherits = FALSE)) return(bands_runtime_cache$value)
  path <- system.file(package = "CEClust")
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  files <- files[grepl("(DESCRIPTION$|/R/|/libs/)", gsub("\\\\", "/", files))]
  value <- digest::digest(list(as.character(utils::packageVersion("CEClust")), tools::md5sum(files)))
  bands_runtime_cache$value <- value
  value
}

bands_rng_save <- function() {
  present <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  list(present = present, seed = if (present) get(".Random.seed", .GlobalEnv) else NULL)
}

bands_rng_restore <- function(state) {
  if (state$present) assign(".Random.seed", state$seed, .GlobalEnv)
  else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
}

bands_score <- function(X, labels, C, lambda = 1, penalty = c("none", "bic")) {
  X <- bands_data(X); penalty <- match.arg(penalty)
  bands_profile(bands_stats(X, bands_labels(labels)), C, lambda, ncol(X), penalty)
}
