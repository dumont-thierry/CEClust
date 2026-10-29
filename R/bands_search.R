bands_add_candidate <- function(bank, candidate) {
  if (!is.null(candidate) && is.null(bank[[candidate$key]])) bank[[candidate$key]] <- candidate
  bank
}

bands_baselines <- function(X, store, data_key) {
  one <- bands_candidate(X, rep(1L, nrow(X)), store, "trivial", data_key)
  endpoint <- bands_candidate(X, seq_len(nrow(X)), store, "saturated_endpoint", data_key)
  result <- list(one, endpoint); names(result) <- c(one$key, endpoint$key); result
}

bands_descent <- function(X, initial, C, lambda, branch, control, penalty, store,
    data_key, anchor, stage, shot, seed, requested_order = NA_integer_) {
  started <- proc.time()[3L]; bound <- if (branch == "free") Inf else C
  checkpoint_key <- digest::digest(list(data_key, initial, bound, lambda, branch,
    control, penalty, anchor, stage, shot, seed, bands_runtime()), algo = "xxhash64")
  checkpoint <- file.path(store$path, "descents", paste0(checkpoint_key, ".rds"))
  if (file.exists(checkpoint)) return(readRDS(checkpoint))
  if (is.null(initial)) fit <- list(valid = FALSE, stop = "unsupported_initial_projection")
  else fit <- tryCatch(bands_fit(X, initial, bound, lambda, penalty, control), error = identity)
  valid <- !inherits(fit, "error") && isTRUE(fit$valid)
  origin <- paste(stage, anchor, branch, shot, sep = ":")
  candidate <- if (valid) bands_candidate(X, fit$labels, store, origin, data_key) else NULL
  log <- data.frame(anchor = anchor, stage = stage, branch = branch, shot = shot,
    C = bound, lambda = lambda, seed = seed,
    initial_key = if (is.null(initial)) NA_character_ else digest::digest(bands_labels(initial)),
    candidate_key = if (valid) candidate$key else NA_character_, valid = valid,
    requested_order = requested_order,
    effective_order = if (is.null(initial)) NA_integer_ else length(unique(initial)),
    stop = if (inherits(fit, "error")) conditionMessage(fit) else fit$stop,
    iterations = if (is.null(fit$iterations)) 0L else fit$iterations,
    seconds = unname(proc.time()[3L] - started))
  result <- list(candidate = candidate, log = log,
    trace = if (inherits(fit, "error") || is.null(fit$trace)) numeric() else fit$trace)
  bands_save(result, checkpoint); result
}

bands_pilot <- function(X, indices, seed, control, penalty, store) {
  started <- proc.time()[3L]; Z <- X[indices, , drop = FALSE]
  data_key <- digest::digest(Z); C0 <- bands_C0(Z)
  candidates <- bands_baselines(Z, store, data_key); logs <- list()
  capacity <- min(nrow(unique(Z)), control$max_order)
  orders <- control$orders[control$orders <= capacity]
  for (order in orders) for (shot in seq_len(control$starts)) {
    initial_seed <- bands_seed(seed, "pilot", order, shot)
    set.seed(initial_seed); initial <- bands_initial(Z, order, shot)
    stats <- bands_stats(Z, initial)
    C <- if (all(stats$rank == ncol(Z)))
      max(exp(-ncol(Z) / 2 * log(2 * pi) - stats$logdet / 2)) * (1 + 1e-8) else C0 * 2^10
    for (branch in control$branches) {
      result <- bands_descent(Z, initial, C, 1, branch, control, penalty, store,
        data_key, 0L, "pilot", shot, initial_seed, order)
      candidates <- bands_add_candidate(candidates, result$candidate)
      logs[[length(logs) + 1L]] <- result$log
    }
  }
  list(candidates = candidates, indices = indices, seed = seed, data_key = data_key,
    n = nrow(Z), d = ncol(Z), C0 = C0, orders = orders,
    skipped_orders = setdiff(control$orders, orders), log = do.call(rbind, logs),
    diagram = bands_geometry(candidates, control, penalty),
    seconds = unname(proc.time()[3L] - started))
}

bands_training_bank <- function(S, job, pilot, anchors, control, penalty, store, runtime) {
  started <- proc.time()[3L]; data_key <- digest::digest(S)
  signature <- digest::digest(list(data_key, job, anchors, control, penalty, runtime))
  path <- file.path(store$path, "banks", paste0(job$name, ".rds"))
  if (file.exists(path)) {
    previous <- readRDS(path)
    if (!identical(previous$signature, signature)) stop("Incompatible bank checkpoint: ", job$name)
    return(previous)
  }
  set.seed(bands_seed(job$seed, "inner_bootstrap"))
  inner_indices <- sample.int(nrow(S), nrow(S), replace = TRUE)
  W <- S[inner_indices, , drop = FALSE]; inner_key <- digest::digest(W)
  cache <- bands_cache(control$cache_mb * 1024^2)
  candidates <- bands_baselines(S, store, data_key); logs <- list()
  for (k in seq_len(nrow(anchors))) {
    point <- anchors[k, ]; source <- pilot$candidates[[point$key]]
    initial <- bands_project(source, W, point$C, point$lambda, store, cache, inner_key)$labels
    for (branch in control$branches) {
      warm <- bands_descent(W, initial, point$C, point$lambda, branch, control, penalty,
        store, inner_key, k, "warm_inner", 0L, job$seed)
      logs[[length(logs) + 1L]] <- warm$log
      if (is.null(warm$candidate)) next
      bound <- if (branch == "free") Inf else point$C
      projected <- bands_project(warm$candidate, S, bound, point$lambda, store, cache, data_key)$labels
      target <- bands_descent(S, projected, point$C, point$lambda, branch, control, penalty,
        store, data_key, k, "warm_target", 0L, job$seed)
      logs[[length(logs) + 1L]] <- target$log
      candidates <- bands_add_candidate(candidates, target$candidate)
    }
    requested <- 2L * point$K
    effective <- min(requested, control$max_order, nrow(unique(S)))
    for (shot in seq_len(control$extra_starts)) {
      initial_seed <- bands_seed(job$seed, "extra", k, shot)
      set.seed(initial_seed); initial <- bands_initial(S, effective, shot)
      for (branch in control$branches) {
        result <- bands_descent(S, initial, point$C, point$lambda, branch, control, penalty,
          store, data_key, k, "extra_target", shot, initial_seed, requested)
        logs[[length(logs) + 1L]] <- result$log
        candidates <- bands_add_candidate(candidates, result$candidate)
      }
    }
  }
  log <- if (length(logs)) do.call(rbind, logs) else data.frame()
  result <- list(name = job$name, kind = job$kind, rep = job$rep, seed = job$seed,
    indices = job$indices, inner_indices = inner_indices, inner_data_key = inner_key,
    data_key = data_key, signature = signature, candidates = candidates, log = log,
    seconds = unname(proc.time()[3L] - started),
    limits = if (nrow(log) && any(log$requested_order > log$effective_order, na.rm = TRUE))
      "initial_order_capped_by_max_order_or_sample_capacity" else character())
  bands_save(result, path); result
}

bands_run_job <- function(job, data, pilot, anchors, control, penalty, store, runtime) {
  S <- if (job$kind == "bootstrap") data[job$indices, , drop = FALSE] else data
  bands_training_bank(S, job, pilot, anchors, control, penalty, store, runtime)
}

bands_run_jobs <- function(jobs, X, pilot, anchors, control, penalty, store, runtime, workers) {
  if (workers == 1L) return(lapply(jobs, bands_run_job, data = X, pilot = pilot,
    anchors = anchors, control = control, penalty = penalty, store = store, runtime = runtime))
  cluster <- parallel::makePSOCKcluster(min(workers, length(jobs)))
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster, bands_worker_setup, paths = .libPaths())
  parallel::parLapplyLB(cluster, jobs, bands_run_job, data = X, pilot = pilot,
    anchors = anchors, control = control, penalty = penalty, store = store, runtime = runtime)
}

bands_worker_setup <- function(paths) {
  .libPaths(paths); loadNamespace("CEClust"); NULL
}

bands_jobs <- function(n, L, B, seed, bootstrap_indices) {
  if (is.null(bootstrap_indices)) {
    set.seed(bands_seed(seed, "outer_bootstrap"))
    bootstrap_indices <- replicate(B, sample.int(n, n, replace = TRUE), simplify = FALSE)
  }
  if (!is.list(bootstrap_indices) || length(bootstrap_indices) != B)
    stop("Supply exactly B bootstrap index vectors.")
  for (index in bootstrap_indices) bands_check_indices(index, n)
  jobs <- vector("list", L + B)
  for (i in seq_len(L + B)) {
    kind <- if (i <= L) "original" else "bootstrap"; rep <- if (i <= L) i else i - L
    jobs[[i]] <- list(name = sprintf("%s_%02d", kind, rep), kind = kind, rep = rep,
      seed = bands_seed(seed, kind, rep),
      indices = if (kind == "original") NULL else as.integer(bootstrap_indices[[rep]]))
  }
  jobs
}

bands_check_indices <- function(indices, n) {
  if (length(indices) != n || anyNA(indices) ||
      any(indices != floor(indices) | indices < 1 | indices > n)) stop("Invalid bootstrap indices.")
  invisible(TRUE)
}
