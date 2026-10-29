# Post-fit removal only: independent banks and stability diagnostics are unchanged.
bands_cleanup_split <- function(X, labels, groups, mode, seed) {
  set.seed(seed); initial <- labels; next_label <- max(labels)
  for (group in groups) {
    rows <- which(labels == group); A <- X[rows, , drop = FALSE]
    unique_rows <- unique(A)
    if (nrow(unique_rows) < 2L) next
    child <- NULL
    if (mode == "kmeans") {
      centers <- unique_rows[sample.int(nrow(unique_rows), 2L), , drop = FALSE]
      fit <- tryCatch(suppressWarnings(stats::kmeans(A, centers, iter.max = 100L,
        algorithm = "Lloyd")), error = identity)
      if (!inherits(fit, "error")) child <- fit$cluster
    }
    if (is.null(child) || length(unique(child)) != 2L) {
      centered <- sweep(A, 2L, colMeans(A), "-")
      direction <- eigen(crossprod(centered), symmetric = TRUE)$vectors[, 1L]
      score <- as.vector(centered %*% direction); values <- sort(unique(score))
      if (length(values) < 2L) next
      cut <- sort(score)[max(1L, floor(length(score) / 2L))]
      if (cut >= max(score)) cut <- values[max(1L, floor(length(values) / 2L))]
      child <- ifelse(score <= cut, 1L, 2L)
    }
    next_label <- next_label + 1L
    initial[rows[child == 2L]] <- next_label
  }
  as.integer(match(initial, unique(initial)))
}

bands_cleanup_job <- function(job) {
  N <- asNamespace("CEClust"); started <- unname(proc.time()[3L])
  path <- file.path(.cleanup_directory, "jobs", paste0(job$id, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  parent <- readRDS(job$parent_file)
  initial <- bands_cleanup_split(.cleanup_X, parent, job$groups, job$mode, job$seed)
  fit <- NULL; reused_from <- NULL
  if (!is.null(job$reuse_file)) {
    old <- readRDS(job$reuse_file)
    if (identical(old$initial, initial) && identical(old$job$C, job$C) &&
        identical(old$job$lambda, job$lambda)) {
      fit <- old$fit; reused_from <- job$reuse_file
    }
  }
  if (is.null(fit)) fit <- tryCatch(N$bands_fit(.cleanup_X, initial, job$C,
    job$lambda, .cleanup_penalty, .cleanup_control), error = identity)
  valid <- !inherits(fit, "error") && isTRUE(fit$valid)
  score <- if (valid) N$bands_profile(N$bands_stats(.cleanup_X, fit$labels),
    job$C, job$lambda, ncol(.cleanup_X), .cleanup_penalty) else NULL
  result <- list(job = job, initial = initial, fit = fit, score = score,
    valid = valid, reused_from = reused_from, seconds = unname(proc.time()[3L]) - started)
  N$bands_save(result, path); result
}

bands_cleanup_reject <- function(reference, candidate, tolerance, saturation_max,
    entropy_tolerance = 1e-8) {
  improved <- reference$J - candidate$J > tolerance * max(1, abs(reference$J))
  higher_entropy <- candidate$entropy - reference$entropy > entropy_tolerance
  saturated <- candidate$saturation > saturation_max
  c(improved = improved, higher_entropy = higher_entropy, saturated = saturated,
    remove = improved && (higher_entropy || saturated))
}

bands_cleanup <- function(fit, directory, workers = 1L, max_rounds = 4L,
    trials = 4L, seed = 1L, score_tolerance = 1e-8, entropy_tolerance = 1e-8, resume = TRUE,
    reuse_directory = NULL) {
  N <- asNamespace("CEClust"); rng <- N$bands_rng_save()
  on.exit(N$bands_rng_restore(rng), add = TRUE)
  if (!inherits(fit, "cec_bands")) stop("Supply a cec_bands result.")
  workers <- N$bands_integer(workers, "workers")
  max_rounds <- N$bands_integer(max_rounds, "max_rounds")
  trials <- N$bands_integer(trials, "trials", 2L); seed <- N$bands_integer(seed, "seed", 0L)
  if (length(score_tolerance) != 1L || !is.finite(score_tolerance) || score_tolerance <= 0)
    stop("score_tolerance must be positive and finite.")
  if (length(entropy_tolerance) != 1L || !is.finite(entropy_tolerance) || entropy_tolerance < 0)
    stop("entropy_tolerance must be finite and nonnegative.")
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) stop("resume must be TRUE or FALSE.")
  started <- unname(proc.time()[3L]); store <- N$bands_store(directory)
  if (identical(store$path, normalizePath(fit$store$path, winslash = "/", mustWork = TRUE)))
    stop("Use a separate cleanup directory.")
  runtime <- N$bands_runtime()
  signature <- digest::digest(list("cleanup-removal-2", fit$signature, fit$tiles,
    runtime, max_rounds, trials, seed, score_tolerance, entropy_tolerance))
  metadata <- list(signature = signature, source_archive = fit$store$path,
    source_runtime = fit$runtime, runtime = runtime, control = fit$control,
    penalty = fit$penalty, max_rounds = max_rounds, trials = trials, seed = seed,
    score_tolerance = score_tolerance, entropy_tolerance = entropy_tolerance,
    data_hash = digest::digest(fit$X, algo = "sha256"),
    rule = "Lower criterion AND (higher H OR saturation above the retention limit), after split and CEM.")
  manifest <- file.path(store$path, "cleanup_manifest.rds")
  if (file.exists(manifest)) {
    if (!resume || readRDS(manifest)$signature != signature) stop("Incompatible cleanup archive.")
  } else N$bands_save(metadata, manifest)
  final <- file.path(store$path, "cleanup_result.rds")
  if (file.exists(final)) return(readRDS(final))
  reuse <- list()
  for (previous_directory in reuse_directory) {
    previous <- readRDS(file.path(previous_directory, "cleanup_manifest.rds"))
    fields <- c("source_archive", "source_runtime", "control", "penalty", "seed")
    if (!all(vapply(fields, function(k) identical(metadata[[k]], previous[[k]]), TRUE)))
      stop("Reuse requires the identical source archive, numerical settings and seed.")
    for (path in list.files(file.path(previous_directory, "jobs"), "[.]rds$", full.names = TRUE)) {
      job <- readRDS(path)$job
      reuse[[paste(job$key, job$round, job$trial)]] <- normalizePath(path, winslash = "/")
    }
  }
  X <- fit$X; parts <- N$bands_partitions(fit); logs <- list(); targets <- list()
  parts$removed <- rep(FALSE, nrow(parts)); parts$reason <- rep("not_refuted_within_budget", nrow(parts))
  parts$competitor_K <- rep(NA_integer_, nrow(parts)); parts$gain <- rep(NA_real_, nrow(parts))
  parts$competitor_H <- rep(NA_real_, nrow(parts)); parts$competitor_saturation <- rep(NA_real_, nrow(parts))
  cluster <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  setup <- function(paths, X, control, penalty, directory) {
    .libPaths(paths); Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1")
    loadNamespace("CEClust")
    assign(".cleanup_X", X, .GlobalEnv); assign(".cleanup_control", control, .GlobalEnv)
    assign(".cleanup_penalty", penalty, .GlobalEnv); assign(".cleanup_directory", directory, .GlobalEnv)
    NULL
  }
  environment(setup) <- baseenv()
  parallel::clusterCall(cluster, setup, .libPaths(), X, fit$control, fit$penalty, store$path)
  parallel::clusterExport(cluster, c("bands_cleanup_split", "bands_cleanup_job"), envir = environment(bands_cleanup))
  for (p in seq_len(nrow(parts))) {
    part <- parts[p, ]; candidate <- fit$reference$candidates[[part$key]]
    reference <- N$bands_candidate_score(candidate, part$C, part$lambda, fit$penalty)
    parent_file <- file.path(store$path, "seeds", paste0(part$key, "_0.rds"))
    N$bands_save(N$bands_load_labels(candidate, fit$store), parent_file)
    used <- 0L
    for (round in seq_len(max_rounds)) {
      parent <- readRDS(parent_file)
      groups <- as.integer(names(Filter(function(ids) nrow(unique(X[ids, , drop = FALSE])) >= 2L,
        split(seq_len(nrow(X)), parent))))
      if (!length(groups)) { parts$reason[p] <- "no_splittable_groups"; break }
      message("Cleanup ", p, "/", nrow(parts), " ", part$name, " K=", part$K,
        "; round ", round, "/", max_rounds, "; seed K=", max(parent))
      jobs <- vector("list", trials)
      for (trial in seq_len(trials)) {
        random_seed <- N$bands_seed(seed, part$key, round, trial); set.seed(random_seed)
        chosen <- if (trial >= trials - 1L) groups else {
          count <- if (trial == 1L) 1L else max(1L, ceiling(length(groups) / 2L))
          groups[sample.int(length(groups), count)]
        }
        jobs[[trial]] <- list(id = digest::digest(list(signature, part$key, round, trial), algo = "xxhash64"),
          part = part$name, key = part$key, round = round, trial = trial, seed = random_seed,
          groups = chosen, mode = if (trial %% 2L) "median" else "kmeans", parent_file = parent_file,
          C = part$C, lambda = part$lambda, reuse_file = reuse[[paste(part$key, round, trial)]])
      }
      values <- parallel::parLapplyLB(cluster, jobs, bands_cleanup_job); used <- round
      rejected <- integer()
      for (trial in seq_along(values)) {
        value <- values[[trial]]; score <- value$score
        reject <- if (value$valid) bands_cleanup_reject(reference, score, score_tolerance,
          fit$control$saturation_max, entropy_tolerance) else
            c(improved = FALSE, higher_entropy = FALSE, saturated = FALSE, remove = FALSE)
        if (reject["remove"]) rejected <- c(rejected, trial)
        logs[[length(logs) + 1L]] <- data.frame(partition = part$name, key = part$key,
          source_K = part$K, C = part$C, lambda = part$lambda, round = round, trial = trial,
          job = value$job$id, mode = value$job$mode, initial_K = max(value$initial), valid = value$valid,
          final_K = if (value$valid) score$K else NA_integer_, J = if (value$valid) score$J else NA_real_,
          H = if (value$valid) score$entropy else NA_real_,
          gain = if (value$valid) reference$J - score$J else NA_real_,
          delta_H = if (value$valid) score$entropy - reference$entropy else NA_real_,
          saturation = if (value$valid) score$saturation else NA_real_,
          improved = reject["improved"], higher_entropy = reject["higher_entropy"],
          saturated = reject["saturated"], removes = reject["remove"], reused = !is.null(value$reused_from),
          iterations = if (value$valid) value$fit$iterations else NA_integer_, seconds = value$seconds)
      }
      if (length(rejected)) {
        j <- rejected[which.min(vapply(values[rejected], function(v) v$score$J, 0))]
        best <- values[[j]]$score; why <- bands_cleanup_reject(reference, best, score_tolerance,
          fit$control$saturation_max, entropy_tolerance)
        parts$removed[p] <- TRUE; parts$reason[p] <- paste(names(why)[why & names(why) != "remove"], collapse = "+")
        parts$competitor_K[p] <- best$K; parts$gain[p] <- reference$J - best$J
        parts$competitor_H[p] <- best$entropy; parts$competitor_saturation[p] <- best$saturation
        break
      }
      all_trials <- seq.int(max(1L, trials - 1L), trials)
      j <- all_trials[which.min(vapply(values[all_trials], function(v) if(v$valid) v$score$J else Inf, 0))]
      parent_file <- file.path(store$path, "seeds", paste0(part$key, "_", round, ".rds"))
      N$bands_save(values[[j]]$initial, parent_file)
    }
    targets[[p]] <- data.frame(name = part$name, K = part$K, rounds = used, removed = parts$removed[p], reason = parts$reason[p])
    N$bands_save(list(partitions = parts, descents = logs, targets = targets), file.path(store$path, "cleanup_checkpoint.rds"))
    utils::write.csv(parts, file.path(store$path, "cleanup_partitions.csv"), row.names = FALSE)
  }
  result <- structure(list(partitions = parts, remaining = parts[!parts$removed, , drop = FALSE],
    descents = if(length(logs)) do.call(rbind, logs) else data.frame(), targets = targets,
    source_archive = fit$store$path, store = store, signature = signature, rule = metadata$rule,
    limits = "Finite subdivision search at original representatives only; not a global optimum test.",
    elapsed_seconds = unname(proc.time()[3L]) - started, diagnostics_recomputed = FALSE), class = "bands_cleanup")
  N$bands_save(result, final)
  utils::write.csv(result$descents, file.path(store$path, "cleanup_descents.csv"), row.names = FALSE)
  result
}
