# Grid experiments for the CEC lambda path as a function of the bound C.

####
# Fonctions
####
{

  ###
  # Fonctions utilitaires
  ###
  {

    cec_grid_script_start_dirs <- function() {
      frames <- sys.frames()
      dirs <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
      env_root <- Sys.getenv("CEC_ARTICLE_PROJECT_ROOT", "")
      known_roots <- character(0)

      for (i in rev(seq_along(frames))) {
        ofile <- frames[[i]]$ofile
        if (!is.null(ofile)) {
          dirs <- unique(c(
            dirs,
            dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
          ))
        }
      }

      if (nzchar(env_root) && dir.exists(env_root)) {
        dirs <- unique(c(dirs, normalizePath(env_root, winslash = "/", mustWork = TRUE)))
      }

      known_roots <- known_roots[dir.exists(known_roots)]
      if (length(known_roots) > 0L) {
        dirs <- unique(c(dirs, normalizePath(known_roots, winslash = "/", mustWork = TRUE)))
      }

      unique(dirs)
    }


    find_cec_grid_project_root <- function() {
      for (start_dir in cec_grid_script_start_dirs()) {
        candidates <- unique(c(
          start_dir,
          dirname(start_dir),
          dirname(dirname(start_dir))
        ))

        for (candidate in candidates) {
          marker <- file.path(candidate, "R", "simulation_article_similarity_1d.R")
          if (file.exists(marker)) {
            return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
          }
        }
      }

      normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    }


    cec_grid_project_path <- function(...) {
      file.path(find_cec_grid_project_root(), ...)
    }


    cec_grid_simulation_dir <- function(...) {
      cec_grid_project_path("R", "article_simulations", ...)
    }


    cec_grid_compact_token <- function(x, n = 16L) {
      raw <- serialize(x, connection = NULL, ascii = FALSE, version = 2)
      hex <- paste(sprintf("%02x", as.integer(raw)), collapse = "")
      substr(hex, 1L, min(nchar(hex), n))
    }


    cec_grid_clean_token <- function(x) {
      x <- format(x, scientific = FALSE, trim = TRUE)
      gsub("[^A-Za-z0-9]+", "p", x)
    }


    cec_grid_force_checkpoint_dir <- function(base_dir) {
      parent_dir <- dirname(base_dir)
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      file.path(parent_dir, paste0("fr_", timestamp))
    }


    cec_grid_checkpoint_dir <- function(
      n,
      data_seed,
      algo_seed,
      lambda_grid,
      C_grid,
      r0,
      k0,
      B,
      Nshots_fresh,
      Nshots_warm,
      Nloop
    ) {
      token <- cec_grid_compact_token(list(
        kind = "C_lambda_grid_1d",
        n = n,
        data_seed = data_seed,
        algo_seed = algo_seed,
        lambda_grid = lambda_grid,
        C_grid = C_grid,
        r0 = r0,
        k0 = k0,
        B = B,
        Nshots_fresh = Nshots_fresh,
        Nshots_warm = Nshots_warm,
        Nloop = Nloop
      ))

      cec_grid_simulation_dir(paste0(
        "grid_cp_n", n,
        "_a", algo_seed,
        "_", token
      ))
    }


    load_ceclust_for_grid <- function() {
      invisible(tryCatch(utils::packageVersion("CEClust"), error = function(e) "dev"))
    }


    ceclust_fun <- function(name) {
      get(name, envir = asNamespace("CEClust"), inherits = FALSE)
    }


    set_ceclust_fast_backend <- function(enabled = TRUE, rebuild = FALSE) {
      ns <- asNamespace("CEClust")
      if (exists("CECset_fast_backend", envir = ns, inherits = FALSE)) {
        fun <- get("CECset_fast_backend", envir = ns, inherits = FALSE)
        fun(enabled = enabled, rebuild = rebuild)
      }

      invisible(TRUE)
    }


    article_lambda_grid <- function(from = 0.05, to = 2, by = 0.05) {
      round(seq(from, to, by = by), 10)
    }


    simulate_gaussian_mixture_1d <- function(
      n = 5000,
      seed = 13,
      mu = 2 * sqrt(3) * c(
        -3,
        -1.75,
        -0.75,
        0,
        0.75,
        1.75,
        3
      ),
      sd_vec = c(0.2, 1, 0.5, 1, 0.5, 1, 0.2),
      pi_vec = c(0.5, 2, 0.5, 1, 0.5, 2, 0.5)
    ) {
      set.seed(seed)
      pi_vec <- pi_vec / sum(pi_vec)
      comp <- sample.int(length(mu), size = n, replace = TRUE, prob = pi_vec)
      Z <- stats::rnorm(n = n, mean = mu[comp], sd = sd_vec[comp])

      list(
        Z = Z,
        z_true = comp,
        mu = mu,
        sd_vec = sd_vec,
        pi_vec = pi_vec
      )
    }


    gaussian_mixture_1d_density <- function(x, mu, sd_vec, pi_vec) {
      pi_vec <- pi_vec / sum(pi_vec)
      dens <- numeric(length(x))

      for (j in seq_along(mu)) {
        dens <- dens + pi_vec[j] * stats::dnorm(x, mean = mu[j], sd = sd_vec[j])
      }

      dens
    }


    extract_cec_entropy <- function(obj) {
      candidates <- c("Hphi", "H", "H_projected")
      for (nm in candidates) {
        if (!is.null(obj[[nm]]) && length(obj[[nm]]) == 1L) {
          return(as.numeric(obj[[nm]]))
        }
      }

      NA_real_
    }


    extract_cec_partition <- function(obj) {
      if (!is.null(obj$clusters)) {
        return(obj$clusters)
      }
      if (!is.null(obj$partition)) {
        return(obj$partition)
      }
      if (!is.null(obj$phi)) {
        return(as.integer(as.factor(obj$phi)))
      }

      stop("Could not find a partition in the CECclassif output.")
    }


    cec_grid_fit_params <- function(obj) {
      candidates <- list(
        if (!is.null(obj$fit) && !is.null(obj$fit$params)) obj$fit$params else NULL,
        if (!is.null(obj$params)) obj$params else NULL,
        if (!is.null(obj$projection) && !is.null(obj$projection$params)) obj$projection$params else NULL,
        if (!is.null(obj$proj) && !is.null(obj$proj$params)) obj$proj$params else NULL
      )

      for (params in candidates) {
        if (!is.null(params)) {
          return(params)
        }
      }

      NULL
    }


    cec_grid_safe_col <- function(x, nm, default = NA_real_) {
      if (is.null(x) || !(nm %in% names(x))) {
        return(rep(default, if (is.null(x)) 0L else nrow(x)))
      }
      x[[nm]]
    }


    cec_grid_cell_bounds <- function(x) {
      x <- sort(unique(as.numeric(x)))
      if (length(x) == 1L) {
        width <- max(0.05, abs(x) * 0.05)
        return(data.frame(value = x, lower = x - width, upper = x + width))
      }

      mid <- (x[-1L] + x[-length(x)]) / 2
      lower <- c(x[1L] - (mid[1L] - x[1L]), mid)
      upper <- c(mid, x[length(x)] + (x[length(x)] - mid[length(mid)]))

      data.frame(value = x, lower = lower, upper = upper)
    }


    cec_grid_row_nearest <- function(grid_summary, C, lambda) {
      d <- abs(grid_summary$C - C) / max(diff(range(grid_summary$C)), .Machine$double.eps) +
        abs(grid_summary$lambda - lambda) / max(diff(range(grid_summary$lambda)), .Machine$double.eps)

      which.min(d)
    }


    cec_grid_format_duration_seconds <- function(secs) {
      secs <- max(0, round(secs))
      sprintf(
        "%02d:%02d:%02d",
        secs %/% 3600,
        (secs %% 3600) %/% 60,
        secs %% 60
      )
    }


    cec_grid_format_elapsed_time <- function(start_time) {
      secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      cec_grid_format_duration_seconds(secs)
    }


    cec_grid_batch_size_is_all <- function(batch_size) {
      if (is.null(batch_size)) {
        return(FALSE)
      }
      if (is.character(batch_size)) {
        return(tolower(batch_size[1L]) %in% c("all", "full", "queue"))
      }
      isTRUE(is.infinite(batch_size[1L]))
    }


    cec_grid_effective_batch_size <- function(
      batch_size,
      n_tasks,
      n_cores = 1L,
      dynamic_task_queue = TRUE
    ) {
      n_cores <- max(1L, as.integer(n_cores[1L]))
      if (isTRUE(dynamic_task_queue) && n_cores > 1L && n_tasks > 0L) {
        return(Inf)
      }
      batch_size
    }


    cec_grid_batches_per_C <- function(k0, B, batch_size = NULL, n_cores = 1L) {
      n_tasks <- as.integer(k0) + as.integer(B)
      if (n_tasks <= 0L) {
        return(0L)
      }
      n_cores <- max(1L, as.integer(n_cores[1L]))
      if (is.null(batch_size)) {
        batch_size <- min(n_cores, max(1L, n_tasks))
      } else if (cec_grid_batch_size_is_all(batch_size)) {
        batch_size <- n_tasks
      }
      batch_size <- suppressWarnings(as.integer(batch_size[1L]))
      if (length(batch_size) == 0L || is.na(batch_size) || batch_size < 1L) {
        batch_size <- min(n_cores, max(1L, n_tasks))
      }

      if (n_cores > 1L) {
        ceiling(n_tasks / batch_size)
      } else {
        n_tasks
      }
    }


    cec_grid_new_global_progress <- function(
      C_grid,
      k0,
      B,
      batch_size = NULL,
      n_cores = 1L,
      enabled = TRUE
    ) {
      state <- new.env(parent = emptyenv())
      state$enabled <- isTRUE(enabled)
      state$start_time <- Sys.time()
      state$batch_done <- 0L
      state$task_done <- 0L
      state$current_C_task_done <- 0L
      state$tasks_per_C <- as.integer(k0) + as.integer(B)
      state$batches_per_C <- cec_grid_batches_per_C(
        k0 = k0,
        B = B,
        batch_size = batch_size,
        n_cores = n_cores
      )
      state$total_C <- length(C_grid)
      state$total_batches <- state$batches_per_C * state$total_C
      state$total_tasks <- state$tasks_per_C * state$total_C
      state$current_C <- NA_real_
      state$current_C_index <- NA_integer_
      state
    }


    cec_grid_set_global_progress_C <- function(state, C_value, C_index) {
      if (is.null(state)) {
        return(invisible(NULL))
      }
      state$current_C <- C_value
      state$current_C_index <- C_index
      state$current_C_task_done <- 0L
      invisible(state)
    }


    cec_grid_print_global_batch_progress <- function(state, n_tasks = NA_integer_) {
      if (is.null(state) || !isTRUE(state$enabled)) {
        return(invisible(NULL))
      }

      state$batch_done <- state$batch_done + 1L
      task_txt <- if (is.finite(n_tasks)) {
        paste0(" (", n_tasks, " task(s))")
      } else {
        ""
      }
      C_txt <- if (is.finite(state$current_C) && is.finite(state$current_C_index)) {
        paste0(
          " | C ", format(state$current_C),
          " (", state$current_C_index, "/", state$total_C, ")"
        )
      } else {
        ""
      }

      cat(
        "CECdiagnose batch ",
        state$batch_done,
        "/",
        state$total_batches,
        task_txt,
        C_txt,
        " | elapsed ",
        cec_grid_format_elapsed_time(state$start_time),
        "\n",
        sep = ""
      )

      invisible(state)
    }


    cec_grid_task_label <- function(task) {
      bits <- c(task$type, task$idx, task$direction)
      bits <- bits[!vapply(bits, is.null, logical(1))]
      bits <- as.character(bits)
      bits <- bits[nzchar(bits)]
      paste(bits, collapse = " ")
    }


    cec_grid_print_global_task_progress <- function(state, task = NULL) {
      if (is.null(state) || !isTRUE(state$enabled)) {
        return(invisible(NULL))
      }

      state$task_done <- state$task_done + 1L
      state$current_C_task_done <- state$current_C_task_done + 1L

      elapsed_secs <- as.numeric(difftime(Sys.time(), state$start_time, units = "secs"))
      eta_txt <- ""
      if (state$task_done > 0L && is.finite(state$total_tasks) && state$total_tasks > state$task_done &&
          is.finite(elapsed_secs) && elapsed_secs > 0) {
        rate <- state$task_done / elapsed_secs
        remaining_secs <- (state$total_tasks - state$task_done) / max(rate, .Machine$double.eps)
        eta_txt <- paste0(" | eta ", cec_grid_format_duration_seconds(remaining_secs))
      }

      C_txt <- if (is.finite(state$current_C) && is.finite(state$current_C_index)) {
        paste0(
          " | C ", format(state$current_C),
          " (", state$current_C_index, "/", state$total_C, ")"
        )
      } else {
        ""
      }

      within_C_txt <- if (is.finite(state$tasks_per_C) && state$tasks_per_C > 0L) {
        paste0(" | C-task ", state$current_C_task_done, "/", state$tasks_per_C)
      } else {
        ""
      }

      last_txt <- if (!is.null(task)) {
        paste0(" | last=", cec_grid_task_label(task))
      } else {
        ""
      }

      cat(
        "CECdiagnose progress: task ",
        state$task_done,
        "/",
        state$total_tasks,
        C_txt,
        within_C_txt,
        last_txt,
        " | elapsed ",
        cec_grid_format_elapsed_time(state$start_time),
        eta_txt,
        "\n",
        sep = ""
      )

      invisible(state)
    }


    with_cec_grid_global_batch_progress <- function(expr, state = NULL) {
      if (is.null(state) || !isTRUE(state$enabled)) {
        return(eval.parent(substitute(expr)))
      }

      ns <- asNamespace("CEClust")
      parallel_name <- "CECrun_diagnostic_tasks_parallel_batch"
      sequential_name <- "CECrun_diagnostic_tasks_sequential"
      mark_name <- "CECdiagnostics_mark_task_completed"

      if (!exists(parallel_name, envir = ns, inherits = FALSE) ||
          !exists(mark_name, envir = ns, inherits = FALSE)) {
        return(eval.parent(substitute(expr)))
      }

      original_parallel <- get(parallel_name, envir = ns, inherits = FALSE)
      original_mark <- get(mark_name, envir = ns, inherits = FALSE)
      original_sequential <- if (exists(sequential_name, envir = ns, inherits = FALSE)) {
        get(sequential_name, envir = ns, inherits = FALSE)
      } else {
        NULL
      }

      wrapped_parallel <- function(tasks, context, n_cores, progress_state, ...) {
        cec_grid_print_global_batch_progress(state, n_tasks = length(tasks))
        original_parallel(
          tasks = tasks,
          context = context,
          n_cores = n_cores,
          progress_state = progress_state,
          ...
        )
      }

      wrapped_mark <- function(progress_state, task) {
        out <- original_mark(progress_state, task)
        cec_grid_print_global_task_progress(state, task = task)
        out
      }

      wrapped_sequential <- function(tasks, context, progress_state, ...) {
        original_sequential(
          tasks = tasks,
          context = context,
          progress_state = progress_state,
          ...
        )
      }

      unlockBinding(parallel_name, ns)
      assign(parallel_name, wrapped_parallel, envir = ns)
      lockBinding(parallel_name, ns)

      unlockBinding(mark_name, ns)
      assign(mark_name, wrapped_mark, envir = ns)
      lockBinding(mark_name, ns)

      on.exit({
        unlockBinding(mark_name, ns)
        assign(mark_name, original_mark, envir = ns)
        lockBinding(mark_name, ns)

        unlockBinding(parallel_name, ns)
        assign(parallel_name, original_parallel, envir = ns)
        lockBinding(parallel_name, ns)
      }, add = TRUE)

      eval.parent(substitute(expr))
    }

  }


  ###
  # Fonctions principales
  ###
  {

    cec_grid_bound_saturation_table <- function(best_parts) {
      out <- data.frame(
        lambda = best_parts$summary$lambda,
        sat = NA_real_,
        n_sat_clusters = NA_integer_,
        any_bound_saturation = NA,
        stringsAsFactors = FALSE
      )

      for (i in seq_along(best_parts$best)) {
        obj <- best_parts$best[[i]]
        if (is.null(obj)) {
          next
        }

        params <- cec_grid_fit_params(obj)
        if (is.null(params)) {
          next
        }

        nu <- suppressWarnings(as.numeric(params$nu))
        if (length(nu) == 0L || all(!is.finite(nu))) {
          next
        }

        nu[!is.finite(nu)] <- 0
        nu_total <- sum(nu, na.rm = TRUE)
        if (is.finite(nu_total) && nu_total > 0) {
          nu <- nu / nu_total
        }

        sat_idx <- unique(suppressWarnings(as.integer(params$functionBoundReached)))
        sat_idx <- sat_idx[is.finite(sat_idx)]
        sat_idx <- sat_idx[sat_idx >= 1L & sat_idx <= length(nu)]

        sat_share <- if (length(sat_idx) == 0L) 0 else sum(nu[sat_idx], na.rm = TRUE)
        sat_share <- min(1, max(0, sat_share))

        out$sat[i] <- sat_share
        out$n_sat_clusters[i] <- length(sat_idx)
        out$any_bound_saturation[i] <- length(sat_idx) > 0L
      }

      out
    }


    cec_grid_add_bound_saturation_to_best_parts <- function(best_parts) {
      sat_df <- cec_grid_bound_saturation_table(best_parts)
      idx <- match(best_parts$summary$lambda, sat_df$lambda)
      best_parts$summary$sat <- sat_df$sat[idx]
      best_parts$summary$n_sat_clusters <- sat_df$n_sat_clusters[idx]
      best_parts$summary$any_bound_saturation <- sat_df$any_bound_saturation[idx]
      best_parts
    }


    cec_grid_add_bound_saturation_to_lambda_diag <- function(lambda_diag, best_parts) {
      sat_df <- cec_grid_bound_saturation_table(best_parts)
      idx <- match(lambda_diag$summary$lambda, sat_df$lambda)
      lambda_diag$summary$sat <- sat_df$sat[idx]
      lambda_diag$summary$n_sat_clusters <- sat_df$n_sat_clusters[idx]
      lambda_diag$summary$any_bound_saturation <- sat_df$any_bound_saturation[idx]
      lambda_diag
    }


    cec_grid_lambda_runs_with_bounds <- function(mask, lambda) {
      run_mask <- ceclust_fun("CECstableRunMask")(
        mask,
        min_consecutive = 1L
      )
      runs <- run_mask$runs
      if (nrow(runs) > 0L) {
        runs$lambda_start <- lambda[runs$start]
        runs$lambda_end <- lambda[runs$end]
      } else {
        runs$lambda_start <- numeric(0)
        runs$lambda_end <- numeric(0)
      }
      runs
    }


    cec_grid_single_lambda_changes <- function(lambda) {
      list(
        change_indices = 1L,
        change_lambdas = lambda,
        summary = data.frame(
          lambda = lambda,
          is_change = TRUE,
          stringsAsFactors = FALSE
        )
      )
    }


    cec_grid_get_fit_param <- function(obj, name, default = NULL) {
      params <- cec_grid_fit_params(obj)
      if (!is.null(params) && !is.null(params[[name]])) {
        return(params[[name]])
      }
      default
    }


    cec_grid_best_index_for_lambda <- function(best_parts, lambda, tol = sqrt(.Machine$double.eps)) {
      if (is.null(best_parts) || is.null(best_parts$summary) || !("lambda" %in% names(best_parts$summary))) {
        return(NA_integer_)
      }
      lambdas <- best_parts$summary$lambda
      if (length(lambdas) == 0L || all(!is.finite(lambdas))) {
        return(NA_integer_)
      }
      idx <- which.min(abs(lambdas - lambda))
      if (length(idx) == 0L || !is.finite(lambdas[idx]) || abs(lambdas[idx] - lambda) > tol) {
        return(NA_integer_)
      }
      idx
    }


    cec_grid_rebuild_best_parts_summary <- function(best_parts) {
      summary_row_fun <- ceclust_fun("CECbestPartitionSummaryRow")
      rows <- lapply(best_parts$best, summary_row_fun)
      summary_out <- do.call(rbind, rows)
      summary_out$C_repaired <- vapply(
        best_parts$best,
        function(obj) if (is.null(obj)) FALSE else isTRUE(obj$C_repaired),
        logical(1)
      )
      summary_out$propagated_from_C <- vapply(
        best_parts$best,
        function(obj) {
          if (is.null(obj) || is.null(obj$propagated_from_C)) NA_real_ else as.numeric(obj$propagated_from_C[1L])
        },
        numeric(1)
      )
      summary_out$repair_origin_C <- vapply(
        best_parts$best,
        function(obj) {
          if (is.null(obj) || is.null(obj$repair_origin_C)) NA_real_ else as.numeric(obj$repair_origin_C[1L])
        },
        numeric(1)
      )
      best_parts$summary <- summary_out
      best_parts
    }


    cec_grid_refresh_result_after_best_parts <- function(
      result_obj,
      Z,
      lambda_grid,
      partition_metric = "match",
      change_sim_threshold = 0.95,
      change_detection_method = "regime_reference",
      detect_changes_on_kept_lambdas = TRUE
    ) {
      meta <- result_obj$meta
      best_parts <- cec_grid_rebuild_best_parts_summary(result_obj$best_parts)
      best_parts <- cec_grid_add_bound_saturation_to_best_parts(best_parts)
      lambda_diag <- cec_grid_add_bound_saturation_to_lambda_diag(result_obj$lambda_diag, best_parts)

      lambda_keep <- identify_cec_grid_lambda_keep(
        lambda_diag = lambda_diag,
        best_parts = best_parts,
        stab_algo_threshold = meta$stab_algo_threshold,
        rsi_threshold = meta$rsi_threshold,
        sat_threshold = meta$sat_threshold,
        use_smoothed = meta$use_smoothed_keep,
        min_consecutive = meta$min_consecutive
      )

      change_lambda_subset <- if (isTRUE(detect_changes_on_kept_lambdas)) {
        lambda_keep$kept_lambdas
      } else {
        NULL
      }

      valid_best <- vapply(best_parts$best, function(obj) !is.null(obj), logical(1))
      valid_idx <- which(valid_best)
      change_valid_idx <- valid_idx
      if (!is.null(change_lambda_subset)) {
        change_valid_idx <- valid_idx[best_parts$summary$lambda[valid_idx] %in% change_lambda_subset]
      }

      if (length(change_valid_idx) <= 1L) {
        only_idx <- if (length(change_valid_idx) == 1L) change_valid_idx else valid_idx[1L]
        only_lambda <- best_parts$summary$lambda[only_idx]
        if (length(only_lambda) == 0L || is.na(only_lambda)) {
          only_lambda <- lambda_grid[1L]
        }
        changes <- cec_grid_single_lambda_changes(only_lambda)
      } else {
        changes <- CEClust::CECdetectPartitionChanges(
          best_parts,
          partition_metric = partition_metric,
          sim_threshold = change_sim_threshold,
          criterion = "REO_or_threshold",
          include_first = TRUE,
          method = change_detection_method,
          lambda_subset = change_lambda_subset
        )
      }

      if (nrow(changes$summary) > 0L) {
        first_lambda <- changes$summary$lambda[1L]
        if (!(first_lambda %in% changes$change_lambdas)) {
          changes$summary$is_change[1L] <- TRUE
          changes$change_indices <- sort(unique(c(1L, changes$change_indices)))
          changes$change_lambdas <- changes$summary$lambda[changes$change_indices]
        }
      }

      result_obj$best_parts <- best_parts
      result_obj$lambda_diag <- lambda_diag
      result_obj$lambda_keep <- lambda_keep
      result_obj$changes <- changes
      result_obj
    }


    cec_grid_repair_C_trajectories <- function(
      runs_by_C,
      Z,
      C_grid,
      lambda_grid,
      familyType,
      tol = 1e-10,
      max_iter = 100L,
      verbose = TRUE
    ) {
      if (length(C_grid) <= 1L || length(lambda_grid) == 0L) {
        return(list(runs_by_C = runs_by_C, repair_log = data.frame(), n_repairs = 0L, converged = TRUE))
      }

      eval_fun <- ceclust_fun("CECevaluatePartitionOnLambda")
      build_fun <- ceclust_fun("CECbuildBestPartitionObjectAtLambda")

      repair_log <- list()
      converged_all <- TRUE

      get_slot_family <- function(C_idx, obj = NULL) {
        fam <- if (!is.null(obj)) cec_grid_get_fit_param(obj, "familyType", default = NULL) else NULL
        if (is.null(fam) || length(fam) == 0L || is.na(fam[1L])) {
          fam <- runs_by_C[[C_idx]]$meta$familyType
        }
        if (is.null(fam) || length(fam) == 0L || is.na(fam[1L])) {
          fam <- familyType
        }
        as.character(fam[1L])
      }

      evaluate_partition_at_slot <- function(partition, C_idx, lambda_value, family = NULL) {
        if (is.null(family)) {
          family <- get_slot_family(C_idx)
        }
        eval_fun(
          partition = partition,
          Z = Z,
          lambda = lambda_value,
          C = C_grid[C_idx],
          familyType = family
        )
      }

      current_slot_H <- function(obj, C_idx, lambda_value) {
        if (is.null(obj)) {
          return(Inf)
        }
        partition <- extract_cec_partition(obj)
        evaluate_partition_at_slot(
          partition = partition,
          C_idx = C_idx,
          lambda_value = lambda_value,
          family = get_slot_family(C_idx, obj)
        )$H_total
      }

      for (lambda_value in lambda_grid) {
        slot <- lapply(seq_along(C_grid), function(C_idx) {
          idx <- cec_grid_best_index_for_lambda(runs_by_C[[C_idx]]$best_parts, lambda_value)
          obj <- if (is.na(idx)) NULL else runs_by_C[[C_idx]]$best_parts$best[[idx]]
          list(C_idx = C_idx, best_idx = idx, obj = obj)
        })

        valid <- which(vapply(slot, function(x) !is.na(x$best_idx) && !is.null(x$obj), logical(1)))
        if (length(valid) <= 1L) {
          next
        }

        slot <- slot[valid]
        n_slot <- length(slot)
        converged_lambda <- TRUE

        for (iter in seq_len(max_iter)) {
          changed_iter <- FALSE

          for (k in seq_len(n_slot - 1L)) {
            from <- slot[[k]]
            to <- slot[[k + 1L]]
            partition_from <- extract_cec_partition(from$obj)
            eval_to <- evaluate_partition_at_slot(
              partition = partition_from,
              C_idx = to$C_idx,
              lambda_value = lambda_value,
              family = get_slot_family(to$C_idx, to$obj)
            )
            H_new <- eval_to$H_total
            H_old <- current_slot_H(to$obj, to$C_idx, lambda_value)

            if (isTRUE(H_new < H_old - tol)) {
              old_obj <- to$obj
              new_obj <- build_fun(
                partition = partition_from,
                lambda = lambda_value,
                Z = Z,
                C = C_grid[to$C_idx],
                familyType = get_slot_family(to$C_idx, to$obj),
                template_obj = from$obj,
                eval_obj = eval_to,
                criterion = if (!is.null(runs_by_C[[to$C_idx]]$best_parts$criterion)) runs_by_C[[to$C_idx]]$best_parts$criterion else NULL
              )
              new_obj$C_repaired <- TRUE
              new_obj$propagated_from_C <- C_grid[from$C_idx]
              new_obj$repair_origin_C <- if (!is.null(from$obj$repair_origin_C)) from$obj$repair_origin_C else C_grid[from$C_idx]
              new_obj$repair_direction_C <- "increasing_C"

              slot[[k + 1L]]$obj <- new_obj
              repair_log[[length(repair_log) + 1L]] <- data.frame(
                lambda = lambda_value,
                iteration = iter,
                direction = "increasing_C",
                from_C = C_grid[from$C_idx],
                to_C = C_grid[to$C_idx],
                old_H = H_old,
                new_H = H_new,
                improvement = H_old - H_new,
                old_REO = if (is.null(old_obj$REO)) NA_integer_ else as.integer(old_obj$REO),
                new_REO = as.integer(new_obj$REO),
                stringsAsFactors = FALSE
              )
              changed_iter <- TRUE
            }
          }

          for (k in n_slot:2L) {
            from <- slot[[k]]
            to <- slot[[k - 1L]]
            partition_from <- extract_cec_partition(from$obj)
            eval_to <- evaluate_partition_at_slot(
              partition = partition_from,
              C_idx = to$C_idx,
              lambda_value = lambda_value,
              family = get_slot_family(to$C_idx, to$obj)
            )
            H_new <- eval_to$H_total
            H_old <- current_slot_H(to$obj, to$C_idx, lambda_value)

            if (isTRUE(H_new < H_old - tol)) {
              old_obj <- to$obj
              new_obj <- build_fun(
                partition = partition_from,
                lambda = lambda_value,
                Z = Z,
                C = C_grid[to$C_idx],
                familyType = get_slot_family(to$C_idx, to$obj),
                template_obj = from$obj,
                eval_obj = eval_to,
                criterion = if (!is.null(runs_by_C[[to$C_idx]]$best_parts$criterion)) runs_by_C[[to$C_idx]]$best_parts$criterion else NULL
              )
              new_obj$C_repaired <- TRUE
              new_obj$propagated_from_C <- C_grid[from$C_idx]
              new_obj$repair_origin_C <- if (!is.null(from$obj$repair_origin_C)) from$obj$repair_origin_C else C_grid[from$C_idx]
              new_obj$repair_direction_C <- "decreasing_C"

              slot[[k - 1L]]$obj <- new_obj
              repair_log[[length(repair_log) + 1L]] <- data.frame(
                lambda = lambda_value,
                iteration = iter,
                direction = "decreasing_C",
                from_C = C_grid[from$C_idx],
                to_C = C_grid[to$C_idx],
                old_H = H_old,
                new_H = H_new,
                improvement = H_old - H_new,
                old_REO = if (is.null(old_obj$REO)) NA_integer_ else as.integer(old_obj$REO),
                new_REO = as.integer(new_obj$REO),
                stringsAsFactors = FALSE
              )
              changed_iter <- TRUE
            }
          }

          if (!changed_iter) {
            break
          }
          if (iter == max_iter) {
            converged_lambda <- FALSE
          }
        }

        converged_all <- converged_all && converged_lambda
        for (entry in slot) {
          runs_by_C[[entry$C_idx]]$best_parts$best[[entry$best_idx]] <- entry$obj
        }
      }

      log_df <- if (length(repair_log) > 0L) do.call(rbind, repair_log) else data.frame()
      if (isTRUE(verbose)) {
        message("C-trajectory repair: ", nrow(log_df), " replacement(s).")
      }

      list(
        runs_by_C = runs_by_C,
        repair_log = log_df,
        n_repairs = nrow(log_df),
        converged = converged_all,
        tol = tol,
        max_iter = max_iter
      )
    }


    cec_grid_repair_lambda_trajectories <- function(
      runs_by_C,
      Z,
      C_grid,
      tol = 1e-10,
      max_iter = 100L,
      verbose = TRUE
    ) {
      repair_log <- list()
      converged_all <- TRUE
      total_repairs <- 0L

      for (i in seq_along(runs_by_C)) {
        repaired <- CEClust::CECrepairBestPartitionTrajectory(
          runs_by_C[[i]]$best_parts,
          Z = Z,
          tol = tol,
          max_iter = max_iter
        )
        n_repairs_i <- if (is.null(repaired$n_repairs)) 0L else as.integer(repaired$n_repairs)
        total_repairs <- total_repairs + n_repairs_i
        converged_all <- converged_all && isTRUE(repaired$converged)
        runs_by_C[[i]]$best_parts <- repaired

        log_i <- repaired$repair_log
        if (!is.null(log_i) && nrow(log_i) > 0L) {
          log_i$C <- C_grid[i]
          repair_log[[length(repair_log) + 1L]] <- log_i
        }
      }

      log_df <- if (length(repair_log) > 0L) do.call(rbind, repair_log) else data.frame()
      if (isTRUE(verbose)) {
        message("Lambda-trajectory repair: ", total_repairs, " replacement(s).")
      }

      list(
        runs_by_C = runs_by_C,
        repair_log = log_df,
        n_repairs = total_repairs,
        converged = converged_all,
        tol = tol,
        max_iter = max_iter
      )
    }


    cec_grid_refresh_all_results <- function(
      runs_by_C,
      Z,
      C_grid,
      lambda_grid,
      partition_metric = "match",
      change_sim_threshold = 0.95,
      change_detection_method = "regime_reference",
      detect_changes_on_kept_lambdas = TRUE
    ) {
      summary_list <- vector("list", length(runs_by_C))
      for (i in seq_along(runs_by_C)) {
        runs_by_C[[i]] <- cec_grid_refresh_result_after_best_parts(
          result_obj = runs_by_C[[i]],
          Z = Z,
          lambda_grid = lambda_grid,
          partition_metric = partition_metric,
          change_sim_threshold = change_sim_threshold,
          change_detection_method = change_detection_method,
          detect_changes_on_kept_lambdas = detect_changes_on_kept_lambdas
        )
        summary_list[[i]] <- cec_grid_summary_for_single_C(runs_by_C[[i]], C_requested = C_grid[i])
      }

      list(runs_by_C = runs_by_C, summary_list = summary_list)
    }


    cec_grid_bind_rows_flexible <- function(rows) {
      rows <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0L, rows)
      if (length(rows) == 0L) {
        return(data.frame())
      }
      all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
      rows <- lapply(rows, function(x) {
        missing <- setdiff(all_names, names(x))
        for (nm in missing) {
          x[[nm]] <- NA
        }
        x[, all_names, drop = FALSE]
      })
      do.call(rbind, rows)
    }


    cec_grid_alternate_repairs <- function(
      runs_by_C,
      Z,
      C_grid,
      lambda_grid,
      familyType,
      C_repair_tol = 1e-10,
      C_repair_max_iter = 100L,
      lambda_repair_tol = 1e-10,
      lambda_repair_max_iter = 100L,
      alternate_max_iter = 20L,
      verbose = TRUE
    ) {
      all_logs <- list()
      summary <- data.frame(
        iteration = integer(0),
        axis = character(0),
        n_repairs = integer(0),
        converged = logical(0),
        stringsAsFactors = FALSE
      )
      converged <- TRUE
      stop_reason <- "max_iter"

      for (iter in seq_len(alternate_max_iter)) {
        C_rep <- cec_grid_repair_C_trajectories(
          runs_by_C = runs_by_C,
          Z = Z,
          C_grid = C_grid,
          lambda_grid = lambda_grid,
          familyType = familyType,
          tol = C_repair_tol,
          max_iter = C_repair_max_iter,
          verbose = verbose
        )
        runs_by_C <- C_rep$runs_by_C
        converged <- converged && isTRUE(C_rep$converged)
        summary <- rbind(summary, data.frame(
          iteration = iter,
          axis = "C",
          n_repairs = C_rep$n_repairs,
          converged = isTRUE(C_rep$converged),
          stringsAsFactors = FALSE
        ))
        if (nrow(C_rep$repair_log) > 0L) {
          log_i <- C_rep$repair_log
          log_i$alternate_iteration <- iter
          log_i$axis <- "C"
          all_logs[[length(all_logs) + 1L]] <- log_i
        }
        if (C_rep$n_repairs == 0L) {
          stop_reason <- "C_stable"
          break
        }

        lambda_rep <- cec_grid_repair_lambda_trajectories(
          runs_by_C = runs_by_C,
          Z = Z,
          C_grid = C_grid,
          tol = lambda_repair_tol,
          max_iter = lambda_repair_max_iter,
          verbose = verbose
        )
        runs_by_C <- lambda_rep$runs_by_C
        converged <- converged && isTRUE(lambda_rep$converged)
        summary <- rbind(summary, data.frame(
          iteration = iter,
          axis = "lambda",
          n_repairs = lambda_rep$n_repairs,
          converged = isTRUE(lambda_rep$converged),
          stringsAsFactors = FALSE
        ))
        if (nrow(lambda_rep$repair_log) > 0L) {
          log_i <- lambda_rep$repair_log
          log_i$alternate_iteration <- iter
          log_i$axis <- "lambda"
          all_logs[[length(all_logs) + 1L]] <- log_i
        }
        if (lambda_rep$n_repairs == 0L) {
          stop_reason <- "lambda_stable"
          break
        }
      }

      log_df <- cec_grid_bind_rows_flexible(all_logs)
      list(
        runs_by_C = runs_by_C,
        repair_log = log_df,
        summary = summary,
        n_repairs = sum(summary$n_repairs),
        converged = converged,
        stop_reason = stop_reason,
        alternate_max_iter = alternate_max_iter
      )
    }


    identify_cec_grid_lambda_keep <- function(
      lambda_diag,
      best_parts = NULL,
      stab_algo_threshold = 0.8,
      rsi_threshold = 0.8,
      sat_threshold = 0.1,
      use_smoothed = FALSE,
      min_consecutive = 1,
      stab_algo_col = "stab0_mean",
      rsi_col = "stab_ratio_mean",
      sat_col = "sat"
    ) {
      lambda_diag_local <- lambda_diag
      if (!(sat_col %in% names(lambda_diag_local$summary))) {
        if (!is.null(best_parts)) {
          best_parts <- cec_grid_add_bound_saturation_to_best_parts(best_parts)
          lambda_diag_local <- cec_grid_add_bound_saturation_to_lambda_diag(
            lambda_diag_local,
            best_parts
          )
        } else {
          lambda_diag_local$summary[[sat_col]] <- NA_real_
          lambda_diag_local$summary$n_sat_clusters <- NA_integer_
          lambda_diag_local$summary$any_bound_saturation <- NA
        }
      }

      out <- CEClust::CECidentifyStableLambdas(
        lambda_diag = lambda_diag_local,
        rule = "both",
        ratio_threshold = rsi_threshold,
        stabB_threshold = stab_algo_threshold,
        sat_threshold = sat_threshold,
        use_smoothed = use_smoothed,
        min_consecutive = min_consecutive,
        best_partitions_obj = best_parts,
        ratio_col = rsi_col,
        stabB_col = stab_algo_col,
        sat_col = sat_col
      )

      summary <- out$summary
      names(summary)[names(summary) == "ratio"] <- "rsi"
      names(summary)[names(summary) == "stabB"] <- "stab_algo"
      names(summary)[names(summary) == "stable"] <- "keep"

      out$kept_lambdas <- out$stable_lambda
      out$rejected_lambdas <- summary$lambda[!summary$keep]
      out$stability_rejected_lambdas <- summary$lambda[!summary$keep_stability]
      out$saturation_rejected_lambdas <- summary$lambda[summary$exclude_saturation]
      out$stable <- out$stable
      out$unstable <- !out$stable
      out$stable_stability <- out$stable_stability
      out$unstable_stability <- !out$stable_stability
      out$stability_runs <- if (!is.null(out$stability_runs)) {
        out$stability_runs
      } else {
        cec_grid_lambda_runs_with_bounds(summary$keep_stability, summary$lambda)
      }
      out$saturation_runs <- cec_grid_lambda_runs_with_bounds(
        !summary$exclude_saturation,
        summary$lambda
      )
      out$summary <- summary
      out$stab_algo_threshold <- stab_algo_threshold
      out$rsi_threshold <- rsi_threshold
      out$sat_threshold <- sat_threshold
      out
    }


    run_cec_grid_single_C <- function(
      Z,
      sim,
      n,
      data_seed,
      algo_seed,
      lambda_grid,
      C,
      r0 = 10,
      k0 = 5,
      B = 5,
      Nshots_fresh = 50,
      Nshots_warm = 50,
      Nloop = 100,
      familyType = "gaussUniv",
      use_cpp = TRUE,
      partition_metric = "match",
      stability_approach = "reference_bestH",
      extract_source = "fixed",
      stab_algo_threshold = 0.8,
      rsi_threshold = 0.8,
      sat_threshold = 0.1,
      use_smoothed_keep = FALSE,
      smooth_k = 3,
      min_consecutive = 1,
      change_detection_method = "regime_reference",
      change_sim_threshold = 0.95,
      detect_changes_on_kept_lambdas = TRUE,
      n_cores = 1L,
      batch_size = NULL,
      checkpoint_dir = NULL,
      auto_checkpoint = FALSE,
      resume = FALSE,
      show_progress = interactive(),
      verbose = TRUE
    ) {
      if (isTRUE(verbose)) {
        message(
          "Running lambda-path for C = ", C,
          " (", length(lambda_grid), " lambdas, k0 = ", k0, ", B = ", B, ")."
        )
      }

      lambda_diag <- CEClust::CECfitLambdaGrid(
        Z = Z,
        lambda_grid = lambda_grid,
        k0 = k0,
        B = B,
        C = C,
        r0 = r0,
        Nshots_fresh = Nshots_fresh,
        Nshots_warm = Nshots_warm,
        Nloop = Nloop,
        familyType = familyType,
        seed = algo_seed,
        silent = TRUE,
        verbose = verbose,
        n_cores = n_cores,
        checkpoint_dir = checkpoint_dir,
        auto_checkpoint = auto_checkpoint,
        resume = resume,
        batch_size = batch_size,
        show_progress = show_progress,
        rerun_failed_serial = TRUE,
        partition_metric = partition_metric,
        stability_approach = stability_approach
      )

      lambda_diag <- CEClust::CECaddSmoothedDiagnostics(
        lambda_diag,
        k = smooth_k
      )

      best_parts <- CEClust::CECextractBestPartitions(
        lambda_diag,
        source = extract_source,
        criterion = "projected_H",
        Z = Z
      )

      best_parts <- CEClust::CECrepairBestPartitionTrajectory(
        best_parts,
        Z
      )
      best_parts <- cec_grid_add_bound_saturation_to_best_parts(best_parts)
      lambda_diag <- cec_grid_add_bound_saturation_to_lambda_diag(lambda_diag, best_parts)

      lambda_keep <- identify_cec_grid_lambda_keep(
        lambda_diag = lambda_diag,
        best_parts = best_parts,
        stab_algo_threshold = stab_algo_threshold,
        rsi_threshold = rsi_threshold,
        sat_threshold = sat_threshold,
        use_smoothed = use_smoothed_keep,
        min_consecutive = min_consecutive
      )

      change_lambda_subset <- if (isTRUE(detect_changes_on_kept_lambdas)) {
        lambda_keep$kept_lambdas
      } else {
        NULL
      }

      valid_best <- vapply(best_parts$best, function(obj) !is.null(obj), logical(1))
      valid_idx <- which(valid_best)
      change_valid_idx <- valid_idx
      if (!is.null(change_lambda_subset)) {
        change_valid_idx <- valid_idx[best_parts$summary$lambda[valid_idx] %in% change_lambda_subset]
      }

      if (length(change_valid_idx) <= 1L) {
        only_idx <- if (length(change_valid_idx) == 1L) change_valid_idx else valid_idx[1L]
        only_lambda <- best_parts$summary$lambda[only_idx]
        if (length(only_lambda) == 0L || is.na(only_lambda)) {
          only_lambda <- lambda_grid[1L]
        }
        changes <- cec_grid_single_lambda_changes(only_lambda)
      } else {
        changes <- CEClust::CECdetectPartitionChanges(
          best_parts,
          partition_metric = partition_metric,
          sim_threshold = change_sim_threshold,
          criterion = "REO_or_threshold",
          include_first = TRUE,
          method = change_detection_method,
          lambda_subset = change_lambda_subset
        )
      }

      if (nrow(changes$summary) > 0L) {
        first_lambda <- changes$summary$lambda[1L]
        if (!(first_lambda %in% changes$change_lambdas)) {
          changes$summary$is_change[1L] <- TRUE
          changes$change_indices <- sort(unique(c(1L, changes$change_indices)))
          changes$change_lambdas <- changes$summary$lambda[changes$change_indices]
        }
      }

      density <- NULL
      if (!is.null(sim$density)) {
        density <- as.data.frame(sim$density)
      } else if (!is.null(sim$mu) && !is.null(sim$sd_vec) && !is.null(sim$pi_vec) &&
          is.numeric(Z) && is.null(dim(Z))) {
        density_grid <- seq(
          min(Z, sim$mu - 4 * sim$sd_vec),
          max(Z, sim$mu + 4 * sim$sd_vec),
          length.out = 2000
        )
        true_density <- gaussian_mixture_1d_density(
          density_grid,
          mu = sim$mu,
          sd_vec = sim$sd_vec,
          pi_vec = sim$pi_vec
        )
        density <- data.frame(x = density_grid, y = true_density)
      }

      list(
        data = sim,
        Z = Z,
        lambda_diag = lambda_diag,
        best_parts = best_parts,
        lambda_keep = lambda_keep,
        changes = changes,
        density = density,
        meta = list(
          n = n,
          data_seed = data_seed,
          algo_seed = algo_seed,
          lambda_grid = lambda_grid,
          C = C,
          C_source = C,
          C_reused = FALSE,
          r0 = r0,
          k0 = k0,
          B = B,
          Nshots_fresh = Nshots_fresh,
          Nshots_warm = Nshots_warm,
          Nloop = Nloop,
          n_cores = n_cores,
          batch_size = batch_size,
          familyType = familyType,
          extract_source = extract_source,
          stab_algo_threshold = stab_algo_threshold,
          rsi_threshold = rsi_threshold,
          sat_threshold = sat_threshold,
          use_smoothed_keep = use_smoothed_keep,
          smooth_k = smooth_k,
          min_consecutive = min_consecutive,
          change_detection_method = change_detection_method,
          change_sim_threshold = change_sim_threshold,
          detect_changes_on_kept_lambdas = detect_changes_on_kept_lambdas
        )
      )
    }


    cec_grid_alias_result_for_C <- function(result_obj, C_requested, C_source) {
      out <- result_obj
      out$meta$C <- C_requested
      out$meta$C_requested <- C_requested
      out$meta$C_source <- C_source
      out$meta$C_reused <- TRUE
      out
    }


    cec_grid_summary_for_single_C <- function(result_obj, C_requested = NULL) {
      if (is.null(C_requested)) {
        C_requested <- result_obj$meta$C
      }

      keep <- result_obj$lambda_keep$summary
      best_summary <- result_obj$best_parts$summary
      diag_summary <- result_obj$lambda_diag$summary
      lambda <- keep$lambda
      idx_best <- match(lambda, best_summary$lambda)
      idx_diag <- match(lambda, diag_summary$lambda)

      REO <- cec_grid_safe_col(best_summary, "REO")[idx_best]
      if (all(is.na(REO))) {
        REO <- vapply(result_obj$best_parts$best, function(obj) {
          if (is.null(obj) || is.null(obj$REO)) NA_real_ else as.numeric(obj$REO)
        }, numeric(1))[idx_best]
      }

      k_hat <- cec_grid_safe_col(best_summary, "k_hat", NA_integer_)[idx_best]
      if (all(is.na(k_hat))) {
        k_hat <- vapply(result_obj$best_parts$best, function(obj) {
          if (is.null(obj) || is.null(obj$partition)) NA_integer_ else length(unique(obj$partition))
        }, integer(1))[idx_best]
      }

      sat <- cec_grid_safe_col(keep, "sat")[match(lambda, keep$lambda)]
      if (all(is.na(sat))) {
        sat <- cec_grid_safe_col(best_summary, "sat")[idx_best]
      }

      n_sat_clusters <- cec_grid_safe_col(keep, "n_sat_clusters", NA_integer_)[match(lambda, keep$lambda)]
      if (all(is.na(n_sat_clusters))) {
        n_sat_clusters <- cec_grid_safe_col(best_summary, "n_sat_clusters", NA_integer_)[idx_best]
      }

      any_bound_saturation <- cec_grid_safe_col(
        best_summary,
        "any_bound_saturation",
        NA
      )[idx_best]

      keep_stability <- if ("keep_stability" %in% names(keep)) {
        keep$keep_stability
      } else {
        cec_grid_safe_col(keep, "keep", FALSE)
      }
      exclude_saturation <- if ("exclude_saturation" %in% names(keep)) {
        keep$exclude_saturation
      } else {
        sat > result_obj$meta$sat_threshold
      }
      keep_all <- if ("keep" %in% names(keep)) {
        keep$keep
      } else {
        keep_stability & !exclude_saturation
      }

      stab_algo <- cec_grid_safe_col(keep, "stab_algo")[match(lambda, keep$lambda)]
      if (all(is.na(stab_algo))) {
        stab_algo <- cec_grid_safe_col(diag_summary, "stab0_mean")[idx_diag]
      }

      rsi <- cec_grid_safe_col(keep, "rsi")[match(lambda, keep$lambda)]
      if (all(is.na(rsi))) {
        rsi <- cec_grid_safe_col(diag_summary, "stab_ratio_mean")[idx_diag]
      }

      fail_reason <- rep("ok", length(lambda))
      fail_reason[!keep_stability] <- "stability_or_rsi"
      fail_reason[exclude_saturation] <- "saturation"
      fail_reason[!keep_stability & exclude_saturation] <- "stability_rsi_and_saturation"

      data.frame(
        C = C_requested,
        C_source = result_obj$meta$C_source,
        C_reused = isTRUE(result_obj$meta$C_reused),
        lambda = lambda,
        REO = REO,
        Hphi = cec_grid_safe_col(best_summary, "Hphi")[idx_best],
        k_hat = k_hat,
        stab_algo = stab_algo,
        rsi = rsi,
        sat = sat,
        n_sat_clusters = n_sat_clusters,
        any_bound_saturation = any_bound_saturation,
        keep_stability = keep_stability,
        exclude_saturation = exclude_saturation,
        keep = keep_all,
        fail_reason = fail_reason,
        stringsAsFactors = FALSE
      )
    }


    run_cec_grid <- function(
      n = 100,
      Z = NULL,
      display_data = NULL,
      density = NULL,
      labels = NULL,
      label_name = "Label",
      dataset_name = NULL,
      data_seed = 13,
      algo_seed = 2026030,
      lambda_grid = article_lambda_grid(0.5, 2, 0.05),
      C_grid = seq(0.2, 10, by = 0.2),
      r0 = 10,
      k0 = 5,
      B = 5,
      Nshots_fresh = 50,
      Nshots_warm = 50,
      Nloop = 100,
      stab_algo_threshold = 0.8,
      rsi_threshold = 0.8,
      sat_threshold = 0.1,
      familyType = "gaussUniv",
      use_cpp = TRUE,
      partition_metric = "match",
      stability_approach = "reference_bestH",
      extract_source = "fixed",
      use_smoothed_keep = FALSE,
      smooth_k = 3,
      min_consecutive = 1,
      change_detection_method = "regime_reference",
      change_sim_threshold = 0.95,
      detect_changes_on_kept_lambdas = TRUE,
      n_cores = 8,
      batch_size = NULL,
      dynamic_task_queue = TRUE,
      checkpoint_dir = FALSE,
      auto_checkpoint = FALSE,
      resume = FALSE,
      force_recompute = FALSE,
      show_progress = interactive(),
      global_progress = TRUE,
      repair_C_trajectories = TRUE,
      C_repair_tol = 1e-10,
      C_repair_max_iter = 100L,
      repair_alternate_trajectories = TRUE,
      lambda_repair_tol = 1e-10,
      lambda_repair_max_iter = 100L,
      repair_alternate_max_iter = 20L,
      output_dir = cec_grid_simulation_dir(),
      save_results = FALSE,
      verbose = TRUE
    ) {
      ceclust_version <- load_ceclust_for_grid()
      set_ceclust_fast_backend(enabled = isTRUE(use_cpp), rebuild = FALSE)

      C_grid <- sort(unique(round(as.numeric(C_grid), 10)))
      lambda_grid <- sort(unique(round(as.numeric(lambda_grid), 10)))
      if (length(C_grid) == 0L || any(!is.finite(C_grid))) {
        stop("C_grid must contain at least one finite value.")
      }
      if (length(lambda_grid) == 0L || any(!is.finite(lambda_grid))) {
        stop("lambda_grid must contain at least one finite value.")
      }

      if (is.null(Z)) {
        sim <- simulate_gaussian_mixture_1d(n = n, seed = data_seed)
        Z <- sim$Z
        if (is.null(dataset_name)) {
          dataset_name <- "gaussian_mixture_1d"
        }
      } else {
        Z_input <- Z
        if (is.numeric(Z_input) && is.null(dim(Z_input))) {
          Z <- as.numeric(Z_input)
          n <- length(Z)
        } else {
          Z_df <- as.data.frame(Z_input)
          if (identical(familyType, "gaussUniv") && ncol(Z_df) == 1L && is.numeric(Z_df[[1L]])) {
            Z <- as.numeric(Z_df[[1L]])
            n <- length(Z)
          } else {
            Z <- Z_df
            n <- nrow(Z)
          }
        }
        sim <- list(Z = Z)
        if (is.null(dataset_name)) {
          dataset_name <- "custom"
        }
      }
      if (!is.null(density)) {
        density <- as.data.frame(density)
        if (!all(c("x", "y") %in% names(density))) {
          stop("density must be a data frame with columns 'x' and 'y'.", call. = FALSE)
        }
        sim$density <- density[, c("x", "y"), drop = FALSE]
      }
      if (is.null(display_data)) {
        display_data <- if (is.data.frame(Z) || is.matrix(Z)) as.data.frame(Z) else data.frame(Z = as.numeric(Z))
      } else {
        display_data <- as.data.frame(display_data)
      }
      if (!is.null(labels) && !(label_name %in% names(display_data))) {
        display_data[[label_name]] <- factor(labels)
      }

      if (is.null(checkpoint_dir) && isTRUE(auto_checkpoint)) {
        checkpoint_dir <- cec_grid_checkpoint_dir(
          n = n,
          data_seed = data_seed,
          algo_seed = algo_seed,
          lambda_grid = lambda_grid,
          C_grid = C_grid,
          r0 = r0,
          k0 = k0,
          B = B,
          Nshots_fresh = Nshots_fresh,
          Nshots_warm = Nshots_warm,
          Nloop = Nloop
        )
      }
      if (isTRUE(force_recompute)) {
        resume <- FALSE
        if (isTRUE(auto_checkpoint) && !isFALSE(checkpoint_dir) && !is.null(checkpoint_dir)) {
          checkpoint_dir <- cec_grid_force_checkpoint_dir(checkpoint_dir)
        }
      }

      if (isTRUE(verbose)) {
        message(
          "Running C x lambda grid ",
          "(n = ", n,
          ", C = ", length(C_grid),
          ", lambdas = ", length(lambda_grid),
          ", k0 = ", k0,
          ", B = ", B, ")."
        )
      }

      effective_batch_size <- cec_grid_effective_batch_size(
        batch_size = batch_size,
        n_tasks = as.integer(k0) + as.integer(B),
        n_cores = n_cores,
        dynamic_task_queue = dynamic_task_queue
      )

      show_global_progress <- isTRUE(global_progress) && isTRUE(show_progress)
      internal_show_progress <- if (isTRUE(show_global_progress)) FALSE else show_progress
      global_progress_state <- cec_grid_new_global_progress(
        C_grid = C_grid,
        k0 = k0,
        B = B,
        batch_size = effective_batch_size,
        n_cores = n_cores,
        enabled = show_global_progress
      )
      if (isTRUE(show_global_progress) && !is.null(checkpoint_dir) && !isFALSE(checkpoint_dir)) {
        cat("CECdiagnose checkpoints root: ", checkpoint_dir, "\n", sep = "")
      }

      runs_by_C <- vector("list", length(C_grid))
      names(runs_by_C) <- cec_grid_clean_token(C_grid)
      summary_list <- vector("list", length(C_grid))
      stop_C <- NA_real_
      stop_index <- NA_integer_
      stop_result <- NULL

      with_cec_grid_global_batch_progress({
        for (i in seq_along(C_grid)) {
          C_value <- C_grid[i]

          if (!is.null(stop_result)) {
            aliased <- cec_grid_alias_result_for_C(
              stop_result,
              C_requested = C_value,
              C_source = stop_C
            )
            runs_by_C[[i]] <- aliased
            summary_list[[i]] <- cec_grid_summary_for_single_C(aliased, C_requested = C_value)
            next
          }

          cec_grid_set_global_progress_C(global_progress_state, C_value, i)

          c_checkpoint_dir <- if (is.null(checkpoint_dir) || isFALSE(checkpoint_dir)) {
            checkpoint_dir
          } else {
            file.path(checkpoint_dir, paste0("C_", cec_grid_clean_token(C_value)))
          }

          single <- run_cec_grid_single_C(
            Z = Z,
            sim = sim,
            n = n,
            data_seed = data_seed,
            algo_seed = algo_seed + i - 1L,
            lambda_grid = lambda_grid,
            C = C_value,
            r0 = r0,
            k0 = k0,
            B = B,
            Nshots_fresh = Nshots_fresh,
            Nshots_warm = Nshots_warm,
            Nloop = Nloop,
            familyType = familyType,
            use_cpp = use_cpp,
            partition_metric = partition_metric,
            stability_approach = stability_approach,
            extract_source = extract_source,
            stab_algo_threshold = stab_algo_threshold,
            rsi_threshold = rsi_threshold,
            sat_threshold = sat_threshold,
            use_smoothed_keep = use_smoothed_keep,
            smooth_k = smooth_k,
            min_consecutive = min_consecutive,
            change_detection_method = change_detection_method,
            change_sim_threshold = change_sim_threshold,
            detect_changes_on_kept_lambdas = detect_changes_on_kept_lambdas,
            n_cores = n_cores,
            batch_size = effective_batch_size,
            checkpoint_dir = c_checkpoint_dir,
            auto_checkpoint = auto_checkpoint,
            resume = resume,
            show_progress = internal_show_progress,
            verbose = verbose
          )
          single$meta$ceclust_version <- as.character(ceclust_version)
          runs_by_C[[i]] <- single
          summary_list[[i]] <- cec_grid_summary_for_single_C(single, C_requested = C_value)

          any_saturated <- any(
            (summary_list[[i]]$any_bound_saturation %in% TRUE) |
              (is.finite(summary_list[[i]]$n_sat_clusters) & summary_list[[i]]$n_sat_clusters > 0L),
            na.rm = TRUE
          )

          if (!isTRUE(any_saturated)) {
            stop_C <- C_value
            stop_index <- i
            stop_result <- single
            if (isTRUE(verbose) && i < length(C_grid)) {
              message(
                "No saturated cluster for C = ", C_value,
                ". Reusing this lambda-path for larger C values."
              )
            }
          }
        }
      }, state = global_progress_state)

      C_repair <- list(repair_log = data.frame(), n_repairs = 0L, converged = TRUE)
      if (isTRUE(repair_C_trajectories)) {
        C_repair <- if (isTRUE(repair_alternate_trajectories)) {
          cec_grid_alternate_repairs(
            runs_by_C = runs_by_C,
            Z = Z,
            C_grid = C_grid,
            lambda_grid = lambda_grid,
            familyType = familyType,
            C_repair_tol = C_repair_tol,
            C_repair_max_iter = C_repair_max_iter,
            lambda_repair_tol = lambda_repair_tol,
            lambda_repair_max_iter = lambda_repair_max_iter,
            alternate_max_iter = repair_alternate_max_iter,
            verbose = verbose
          )
        } else {
          cec_grid_repair_C_trajectories(
            runs_by_C = runs_by_C,
            Z = Z,
            C_grid = C_grid,
            lambda_grid = lambda_grid,
            familyType = familyType,
            tol = C_repair_tol,
            max_iter = C_repair_max_iter,
            verbose = verbose
          )
        }
        runs_by_C <- C_repair$runs_by_C
        refreshed <- cec_grid_refresh_all_results(
          runs_by_C = runs_by_C,
          Z = Z,
          C_grid = C_grid,
          lambda_grid = lambda_grid,
          partition_metric = partition_metric,
          change_sim_threshold = change_sim_threshold,
          change_detection_method = change_detection_method,
          detect_changes_on_kept_lambdas = detect_changes_on_kept_lambdas
        )
        runs_by_C <- refreshed$runs_by_C
        summary_list <- refreshed$summary_list
      }

      grid_summary <- do.call(rbind, summary_list)
      row.names(grid_summary) <- NULL

      out <- list(
        data = sim,
        Z = Z,
        display_data = display_data,
        labels = if (is.null(labels)) NULL else factor(labels),
        label_name = label_name,
        dataset_name = dataset_name,
        C_grid = C_grid,
        lambda_grid = lambda_grid,
        runs_by_C = runs_by_C,
        summary = grid_summary,
        C_repair = C_repair,
        stop_C = stop_C,
        stop_index = stop_index,
        meta = list(
          n = n,
          dataset = dataset_name,
          data_seed = data_seed,
          algo_seed = algo_seed,
          C_grid = C_grid,
          lambda_grid = lambda_grid,
          r0 = r0,
          k0 = k0,
          B = B,
          Nshots_fresh = Nshots_fresh,
          Nshots_warm = Nshots_warm,
          Nloop = Nloop,
          n_cores = n_cores,
          batch_size = batch_size,
          effective_batch_size = effective_batch_size,
          dynamic_task_queue = isTRUE(dynamic_task_queue),
          familyType = familyType,
          ceclust_version = as.character(ceclust_version),
          stab_algo_threshold = stab_algo_threshold,
          rsi_threshold = rsi_threshold,
          sat_threshold = sat_threshold,
          use_smoothed_keep = use_smoothed_keep,
          smooth_k = smooth_k,
          min_consecutive = min_consecutive,
          change_detection_method = change_detection_method,
          change_sim_threshold = change_sim_threshold,
          detect_changes_on_kept_lambdas = detect_changes_on_kept_lambdas,
          global_progress = global_progress,
          repair_C_trajectories = repair_C_trajectories,
          C_repair_tol = C_repair_tol,
          C_repair_max_iter = C_repair_max_iter,
          repair_alternate_trajectories = repair_alternate_trajectories,
          lambda_repair_tol = lambda_repair_tol,
          lambda_repair_max_iter = lambda_repair_max_iter,
          repair_alternate_max_iter = repair_alternate_max_iter
        )
      )

      if (isTRUE(save_results)) {
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }

        stem <- paste0(
          "cec_grid_", cec_grid_clean_token(dataset_name), "_n", n,
          "_k0", k0,
          "_B", B,
          "_", format(Sys.time(), "%Y%m%d_%H%M%S")
        )
        saveRDS(out, file.path(output_dir, paste0(stem, ".rds")))
        utils::write.csv(
          grid_summary,
          file.path(output_dir, paste0(stem, "_summary.csv")),
          row.names = FALSE
        )
      }

      out
    }


    cec_grid_prepare_quantitative_input_data <- function(
      Z,
      dataset_name = "data"
    ) {
      Z_df <- as.data.frame(Z)
      out <- data.frame(row.names = seq_len(nrow(Z_df)))

      for (nm in names(Z_df)) {
        x <- Z_df[[nm]]
        if (is.numeric(x)) {
          out[[nm]] <- x
        } else if (is.ordered(x)) {
          out[[nm]] <- as.numeric(x)
        } else if (is.factor(x) || is.character(x) || is.logical(x)) {
          x_num <- suppressWarnings(as.numeric(as.character(x)))
          bad <- !is.na(x) & is.na(x_num)
          if (any(bad)) {
            stop(
              dataset_name, " contains a non-numeric quantitative column: '", nm,
              "'.",
              call. = FALSE
            )
          }
          out[[nm]] <- x_num
        }
      }

      if (ncol(out) == 0L) {
        stop(dataset_name, " data must contain at least one numeric or coercible quantitative column.", call. = FALSE)
      }

      out
    }


    cec_grid_prepare_quantitative_representation <- function(
      Z,
      quantitative_representation = "raw",
      dataset_name = "data"
    ) {
      rep_arg <- as.character(quantitative_representation)[1L]
      if (length(rep_arg) == 0L || is.na(rep_arg) || !nzchar(rep_arg)) {
        rep_arg <- "raw"
      }
      rep_lower <- tolower(rep_arg)
      Z <- cec_grid_prepare_quantitative_input_data(Z, dataset_name = dataset_name)
      num_cols <- vapply(Z, is.numeric, logical(1))
      if (!any(num_cols)) {
        stop(dataset_name, " data must contain at least one numeric column.", call. = FALSE)
      }

      if (identical(rep_lower, "raw")) {
        out <- Z[, num_cols, drop = FALSE]
      } else if (identical(rep_lower, "normalized")) {
        out <- as.data.frame(scale(Z[, num_cols, drop = FALSE]))
        names(out) <- names(Z)[num_cols]
      } else if (grepl("^pc[0-9]+$", rep_lower)) {
        n_pc <- as.integer(sub("^pc", "", rep_lower))
        pca_fit <- stats::prcomp(Z[, num_cols, drop = FALSE], center = TRUE, scale. = TRUE)
        if (n_pc < 1L || n_pc > ncol(pca_fit$x)) {
          stop("quantitative_representation = '", rep_arg, "' is not available for ", dataset_name, ".", call. = FALSE)
        }
        out <- as.data.frame(pca_fit$x[, seq_len(n_pc), drop = FALSE])
        names(out) <- paste0("PC", seq_len(n_pc))
      } else if (grepl("^rawpc[0-9]+$", rep_lower) || grepl("^pc[0-9]+_raw$", rep_lower)) {
        n_pc <- if (grepl("^rawpc[0-9]+$", rep_lower)) {
          as.integer(sub("^rawpc", "", rep_lower))
        } else {
          as.integer(sub("^pc([0-9]+)_raw$", "\\1", rep_lower))
        }
        pca_fit <- stats::prcomp(Z[, num_cols, drop = FALSE], center = TRUE, scale. = FALSE)
        if (n_pc < 1L || n_pc > ncol(pca_fit$x)) {
          stop("quantitative_representation = '", rep_arg, "' is not available for ", dataset_name, ".", call. = FALSE)
        }
        out <- as.data.frame(pca_fit$x[, seq_len(n_pc), drop = FALSE])
        names(out) <- paste0("rawPC", seq_len(n_pc))
      } else {
        stop(
          "Unknown quantitative_representation: ", rep_arg,
          ". Use 'raw', 'normalized', 'PCk', or 'rawPCk'.",
          call. = FALSE
        )
      }

      list(data = out, representation = rep_lower)
    }


    load_cec_grid_breast_cancer_data <- function() {
      if (!requireNamespace("mlbench", quietly = TRUE)) {
        stop("Package 'mlbench' is required for BreastCancer. Install it with install.packages('mlbench').", call. = FALSE)
      }
      data("BreastCancer", package = "mlbench")
      bc <- get("BreastCancer", envir = environment())
      bc <- bc[stats::complete.cases(bc), , drop = FALSE]
      raw_Z <- bc[, setdiff(names(bc), c("Id", "Class")), drop = FALSE]
      Z <- cec_grid_prepare_quantitative_input_data(raw_Z, dataset_name = "BreastCancer")

      list(
        Z = Z,
        labels = factor(bc$Class),
        label_name = "Class"
      )
    }


    run_cec_grid_iris <- function(
      include_species_as_feature = FALSE,
      quantitative_representation = c("raw", "normalized", "PC2", "PC3", "PC4", "rawPC2", "rawPC3", "rawPC4"),
      algo_seed = 2026030,
      lambda_grid = article_lambda_grid(0.4, 1.8, 0.1),
      C_grid = seq(1.5, 3, by = 0.25),
      r0 = 10,
      k0 = 20,
      B = 20,
      Nshots_fresh = 50,
      Nshots_warm = 50,
      Nloop = 100,
      stab_algo_threshold = 0.8,
      rsi_threshold = 0.8,
      sat_threshold = 0.05,
      familyType = NULL,
      n_cores = 8,
      batch_size = NULL,
      checkpoint_dir = FALSE,
      auto_checkpoint = FALSE,
      resume = FALSE,
      force_recompute = FALSE,
      save_results = FALSE,
      verbose = TRUE,
      ...
    ) {
      iris_data <- datasets::iris
      quantitative_representation <- as.character(quantitative_representation)[1L]
      prep <- cec_grid_prepare_quantitative_representation(
        iris_data[, 1:4, drop = FALSE],
        quantitative_representation,
        dataset_name = "Iris"
      )
      Z <- if (isTRUE(include_species_as_feature)) {
        if (is.null(familyType)) {
          familyType <- "gaussAndDiscreteVector"
        }
        cbind(prep$data, Species = iris_data$Species)
      } else {
        if (is.null(familyType)) {
          familyType <- "gaussVector"
        }
        prep$data
      }
      if (is.null(batch_size)) {
        batch_size <- n_cores
      }

      out <- run_cec_grid(
        n = nrow(iris_data),
        Z = Z,
        display_data = iris_data,
        labels = iris_data$Species,
        label_name = "Species",
        dataset_name = paste0(
          if (isTRUE(include_species_as_feature)) "iris_with_species_" else "iris_",
          cec_grid_clean_token(prep$representation)
        ),
        data_seed = NA_integer_,
        algo_seed = algo_seed,
        lambda_grid = lambda_grid,
        C_grid = C_grid,
        r0 = r0,
        k0 = k0,
        B = B,
        Nshots_fresh = Nshots_fresh,
        Nshots_warm = Nshots_warm,
        Nloop = Nloop,
        stab_algo_threshold = stab_algo_threshold,
        rsi_threshold = rsi_threshold,
        sat_threshold = sat_threshold,
        familyType = familyType,
        n_cores = n_cores,
        batch_size = batch_size,
        checkpoint_dir = checkpoint_dir,
        auto_checkpoint = auto_checkpoint,
        resume = resume,
        force_recompute = force_recompute,
        save_results = save_results,
        verbose = verbose,
        ...
      )
      out$meta$quantitative_representation <- prep$representation
      out$quantitative_representation <- prep$representation
      out$clustering_data <- Z
      out
    }


    run_cec_grid_breast_cancer <- function(
      include_class_as_feature = FALSE,
      quantitative_representation = c("normalized", "raw", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "rawPC2", "rawPC3", "rawPC4", "rawPC5", "rawPC6", "rawPC7", "rawPC8", "rawPC9"),
      algo_seed = 2026030,
      lambda_grid = article_lambda_grid(0.4, 1.8, 0.1),
      C_grid = seq(1.5, 3, by = 0.25),
      r0 = 10,
      k0 = 20,
      B = 20,
      Nshots_fresh = 50,
      Nshots_warm = 50,
      Nloop = 100,
      stab_algo_threshold = 0.8,
      rsi_threshold = 0.8,
      sat_threshold = 0.05,
      familyType = NULL,
      n_cores = 8,
      batch_size = NULL,
      checkpoint_dir = FALSE,
      auto_checkpoint = FALSE,
      resume = FALSE,
      force_recompute = FALSE,
      save_results = FALSE,
      verbose = TRUE,
      ...
    ) {
      bc <- load_cec_grid_breast_cancer_data()
      quantitative_representation <- as.character(quantitative_representation)[1L]
      prep <- cec_grid_prepare_quantitative_representation(
        bc$Z,
        quantitative_representation,
        dataset_name = "BreastCancer"
      )

      display_data <- cbind(bc$Z, Class = bc$labels)
      Z <- if (isTRUE(include_class_as_feature)) {
        if (is.null(familyType)) {
          familyType <- "gaussAndDiscreteVector"
        }
        cbind(prep$data, Class = bc$labels)
      } else {
        if (is.null(familyType)) {
          familyType <- "gaussVector"
        }
        prep$data
      }
      if (is.null(batch_size)) {
        batch_size <- n_cores
      }

      out <- run_cec_grid(
        n = nrow(bc$Z),
        Z = Z,
        display_data = display_data,
        labels = bc$labels,
        label_name = bc$label_name,
        dataset_name = paste0(
          if (isTRUE(include_class_as_feature)) "breast_cancer_with_class_" else "breast_cancer_",
          cec_grid_clean_token(prep$representation)
        ),
        data_seed = NA_integer_,
        algo_seed = algo_seed,
        lambda_grid = lambda_grid,
        C_grid = C_grid,
        r0 = r0,
        k0 = k0,
        B = B,
        Nshots_fresh = Nshots_fresh,
        Nshots_warm = Nshots_warm,
        Nloop = Nloop,
        stab_algo_threshold = stab_algo_threshold,
        rsi_threshold = rsi_threshold,
        sat_threshold = sat_threshold,
        familyType = familyType,
        n_cores = n_cores,
        batch_size = batch_size,
        checkpoint_dir = checkpoint_dir,
        auto_checkpoint = auto_checkpoint,
        resume = resume,
        force_recompute = force_recompute,
        save_results = save_results,
        verbose = verbose,
        ...
      )
      out$meta$quantitative_representation <- prep$representation
      out$quantitative_representation <- prep$representation
      out$clustering_data <- Z
      out
    }


    run_cec_grid_icl <- function(
      dataset_name,
      Z,
      display_data = NULL,
      density = NULL,
      labels = NULL,
      label_name = "Label",
      G = 1:10,
      modelNames = NULL,
      scale_clustering_data = TRUE,
      output_dir = cec_grid_simulation_dir(),
      save_results = TRUE,
      save_pca_plot = FALSE,
      verbose = TRUE
    ) {
      if (!requireNamespace("mclust", quietly = TRUE)) {
        stop("Package 'mclust' is required for ICL comparisons.", call. = FALSE)
      }
      suppressPackageStartupMessages(library(mclust))

      Z <- cec_grid_prepare_quantitative_input_data(Z, dataset_name = dataset_name)
      if (is.null(display_data)) {
        display_data <- Z
      } else {
        display_data <- as.data.frame(display_data)
      }
      if (!is.null(labels) && !(label_name %in% names(display_data))) {
        display_data[[label_name]] <- factor(labels)
      }

      X <- as.matrix(Z)
      ok <- stats::complete.cases(X)
      if (!all(ok)) {
        stop("ICL input contains missing or non-finite rows.", call. = FALSE)
      }
      if (isTRUE(scale_clustering_data)) {
        X <- scale(X)
      }

      if (isTRUE(verbose)) {
        message(
          "Running ", dataset_name, " ICL clustering ",
          "(n = ", nrow(Z), ", p = ", ncol(Z), ", G = ",
          paste(range(G), collapse = ":"), ")."
        )
      }

      icl_grid <- mclust::mclustICL(
        data = X,
        G = G,
        modelNames = modelNames
      )
      icl_summary <- summary(icl_grid)
      best_name <- names(icl_summary)[1L]
      best_tokens <- strsplit(best_name, ",", fixed = TRUE)[[1L]]
      best_modelName <- best_tokens[1L]
      best_G <- as.integer(best_tokens[2L])

      fit <- mclust::Mclust(
        data = X,
        G = best_G,
        modelNames = best_modelName
      )
      partition <- as.integer(fit$classification)
      cluster_weights <- fit$parameters$pro
      if (is.null(cluster_weights)) {
        cluster_weights <- as.numeric(table(partition)) / length(partition)
      }

      out <- list(
        dataset_name = dataset_name,
        method = "ICL",
        data = Z,
        Z = Z,
        display_data = display_data,
        density = if (is.null(density)) NULL else as.data.frame(density),
        labels = if (is.null(labels)) NULL else factor(labels),
        label_name = label_name,
        icl_grid = icl_grid,
        icl_summary = icl_summary,
        icl_fit = fit,
        partition = partition,
        cluster_weights = cluster_weights,
        meta = list(
          dataset = dataset_name,
          method = "ICL",
          n = nrow(Z),
          p = ncol(Z),
          G = G,
          best_G = best_G,
          modelName = best_modelName,
          icl = as.numeric(icl_summary[1L]),
          bic = fit$bic,
          loglik = fit$loglik,
          label_name = label_name,
          scale_clustering_data = isTRUE(scale_clustering_data)
        )
      )

      if (isTRUE(save_results)) {
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }
        stem <- paste0(
          "cec_grid_", cec_grid_clean_token(dataset_name),
          "_icl_", format(Sys.time(), "%Y%m%d_%H%M%S")
        )
        saveRDS(out, file.path(output_dir, paste0(stem, ".rds")))
        utils::write.csv(
          data.frame(
            dataset = dataset_name,
            modelName = best_modelName,
            G = best_G,
            icl = as.numeric(icl_summary[1L]),
            bic = fit$bic,
            loglik = fit$loglik,
            stringsAsFactors = FALSE
          ),
          file.path(output_dir, paste0(stem, "_summary.csv")),
          row.names = FALSE
        )
      }

      out
    }

  }


  ###
  # Fonctions de plots
  ###
  {

    cec_grid_color_palette <- function(n = 100) {
      grDevices::colorRampPalette(c(
        "#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#B2182B"
      ))(n)
    }


    cec_grid_cluster_palette <- function(k) {
      base <- c(
        "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
        "#56B4E9", "#F0E442", "#000000", "#999999", "#882255",
        "#44AA99", "#AA4499", "#117733", "#332288", "#DDCC77"
      )
      if (k <= length(base)) {
        return(base[seq_len(max(1L, k))])
      }
      grDevices::colorRampPalette(base)(k)
    }


    true_uniform_gaussian_lambda_path <- function(
      C = 1,
      lambda_grid = article_lambda_grid(0.05, 2, 0.05),
      m_max = NULL,
      max_m_hard = 10000L
    ) {
      C <- as.numeric(C)[1L]
      lambda_grid <- sort(unique(as.numeric(lambda_grid)))
      if (!is.finite(C) || C <= 0) {
        stop("C must be a positive finite value.", call. = FALSE)
      }
      if (length(lambda_grid) == 0L || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
        stop("lambda_grid must contain positive finite values.", call. = FALSE)
      }

      sigma_C <- 1 / (sqrt(2 * pi) * C)
      if (is.null(m_max)) {
        m_max <- max(
          50L,
          ceiling(10 * C * sqrt(pi / (6 * min(lambda_grid))))
        )
      }
      m_max <- max(1L, as.integer(m_max[1L]))
      max_m_hard <- max(m_max, as.integer(max_m_hard[1L]))

      compute_for_mmax <- function(current_m_max) {
        m_vals <- seq_len(current_m_max)
        ell <- 1 / m_vals
        sigma_hat <- pmax(ell / sqrt(12), sigma_C)
        cross_entropy <- 0.5 * log(2 * pi * sigma_hat^2) +
          ell^2 / (24 * sigma_hat^2)

        rows <- lapply(lambda_grid, function(lambda) {
          objective <- lambda * log(m_vals) + cross_entropy
          idx <- which.min(objective)
          data.frame(
            C = C,
            lambda = lambda,
            m = m_vals[idx],
            H_lambda = objective[idx],
            sigma_C = sigma_C,
            interval_length = 1 / m_vals[idx],
            saturated = (1 / m_vals[idx]) < sqrt(12) * sigma_C,
            stringsAsFactors = FALSE
          )
        })
        do.call(rbind, rows)
      }

      repeat {
        out <- compute_for_mmax(m_max)
        if (!any(out$m >= m_max) || m_max >= max_m_hard) {
          break
        }
        m_max <- min(max_m_hard, 2L * m_max)
      }

      row.names(out) <- NULL
      out$m_max <- m_max
      out
    }


    plot_true_uniform_gaussian_lambda_path <- function(
      C = 1,
      lambda_grid = article_lambda_grid(0.05, 2, 0.05),
      Z = NULL,
      n = 1000,
      seed = 13,
      m_max = NULL,
      show_density_panel = FALSE,
      shade_lambda_below = NULL,
      shade_col = grDevices::adjustcolor("grey40", alpha.f = 0.48),
      shade_hatch = TRUE,
      shade_hatch_density = 12,
      shade_hatch_angle = 45,
      shade_hatch_col = grDevices::adjustcolor("grey10", alpha.f = 0.75),
      histogram_col = grDevices::adjustcolor("grey75", alpha.f = 0.45),
      density_col = "#222222",
      cluster_alpha = 0.62,
      main_top = "Uniform sample and density",
      main_path = expression("True " * lambda * "-path")
    ) {
      path <- true_uniform_gaussian_lambda_path(
        C = C,
        lambda_grid = lambda_grid,
        m_max = m_max
      )

      if (is.null(Z)) {
        set.seed(seed)
        Z <- stats::runif(n, min = 0, max = 1)
      }
      Z <- as.numeric(Z)
      Z <- Z[is.finite(Z)]

      l_bounds <- cec_grid_cell_bounds(path$lambda)
      pal <- cec_grid_cluster_palette(max(path$m, na.rm = TRUE))

      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(oldpar), add = TRUE)
      if (isTRUE(show_density_panel)) {
        graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(1.05, 3.2))

        graphics::par(mar = c(0.6, 4.3, 2.8, 1.1))
        hist_obj <- graphics::hist(
          Z,
          breaks = max(8L, min(80L, ceiling(sqrt(length(Z))))),
          plot = FALSE
        )
        ymax <- max(c(hist_obj$density, 1), na.rm = TRUE)
        graphics::plot(
          NA,
          NA,
          xlim = c(0, 1),
          ylim = c(0, 1.08 * ymax),
          axes = FALSE,
          xlab = "",
          ylab = "density",
          main = main_top
        )
        graphics::plot(hist_obj, freq = FALSE, add = TRUE, col = histogram_col, border = "white")
        graphics::segments(0, 1, 1, 1, col = density_col, lwd = 2)
        graphics::rug(Z, col = grDevices::adjustcolor("grey20", alpha.f = 0.55))
        graphics::axis(2, las = 2)
        graphics::box()
      }

      graphics::par(mar = c(4.3, 4.3, 2.8, 1.1))
      graphics::plot(
        NA,
        NA,
        xlim = c(0, 1),
        ylim = range(l_bounds$lower, l_bounds$upper),
        xlab = "z",
        ylab = expression(lambda),
        main = main_path,
        axes = FALSE
      )

      for (i in seq_len(nrow(path))) {
        lb <- l_bounds[l_bounds$value == path$lambda[i], , drop = FALSE]
        xs <- seq(0, 1, length.out = path$m[i] + 1L)
        for (j in seq_len(path$m[i])) {
          graphics::rect(
            xleft = xs[j],
            ybottom = lb$lower,
            xright = xs[j + 1L],
            ytop = lb$upper,
            col = grDevices::adjustcolor(pal[j], alpha.f = cluster_alpha),
            border = NA
          )
        }
      }

      if (!is.null(shade_lambda_below)) {
        shade_lambda_below <- as.numeric(shade_lambda_below)[1L]
        if (is.finite(shade_lambda_below)) {
          next_band_lower <- l_bounds$lower[l_bounds$value > shade_lambda_below]
          shade_top <- if (length(next_band_lower) > 0L) {
            min(next_band_lower, na.rm = TRUE)
          } else {
            max(l_bounds$upper, na.rm = TRUE)
          }
          graphics::rect(
            xleft = 0,
            ybottom = min(l_bounds$lower, na.rm = TRUE),
            xright = 1,
            ytop = min(shade_top, max(l_bounds$upper, na.rm = TRUE)),
            col = shade_col,
            border = NA
          )
          if (isTRUE(shade_hatch)) {
            graphics::rect(
              xleft = 0,
              ybottom = min(l_bounds$lower, na.rm = TRUE),
              xright = 1,
              ytop = min(shade_top, max(l_bounds$upper, na.rm = TRUE)),
              density = shade_hatch_density,
              angle = shade_hatch_angle,
              col = shade_hatch_col,
              border = NA
            )
          }
        }
      }

      graphics::axis(1)
      graphics::axis(2, at = pretty(path$lambda), las = 2)
      graphics::box()
      invisible(path)
    }


    save_true_uniform_gaussian_lambda_path_plot <- function(
      C = 1,
      plot_file = file.path(
        cec_grid_simulation_dir(),
        paste0("lambda_path_uniform_true_C", cec_grid_clean_token(C), ".png")
      ),
      width = 1600,
      height = 1000,
      res = 140,
      ...
    ) {
      plot_dir <- dirname(plot_file)
      if (!dir.exists(plot_dir)) {
        dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
      }
      grDevices::png(plot_file, width = width, height = height, res = res)
      on.exit(grDevices::dev.off(), add = TRUE)
      path <- plot_true_uniform_gaussian_lambda_path(C = C, ...)
      invisible(list(plot_file = plot_file, path = path))
    }


    cec_grid_display_data <- function(grid_obj) {
      candidates <- list(
        grid_obj$display_data,
        grid_obj$Z_with_species,
        grid_obj$Z,
        grid_obj$data
      )

      for (x in candidates) {
        if (is.data.frame(x)) {
          return(as.data.frame(x))
        }
        if (is.matrix(x)) {
          return(as.data.frame(x))
        }
        if (is.numeric(x) && is.null(dim(x))) {
          return(data.frame(Z = as.numeric(x)))
        }
      }

      if (!is.null(grid_obj$data$Z)) {
        return(data.frame(Z = as.numeric(grid_obj$data$Z)))
      }

      data.frame()
    }


    cec_grid_quant_vars <- function(display_data) {
      names(display_data)[vapply(display_data, is.numeric, logical(1))]
    }


    cec_grid_qual_vars <- function(display_data) {
      names(display_data)[vapply(display_data, function(x) {
        is.factor(x) || is.character(x) || is.logical(x) || is.ordered(x)
      }, logical(1))]
    }


    cec_grid_pca_scores <- function(display_data) {
      qvars <- cec_grid_quant_vars(display_data)
      if (length(qvars) < 2L) {
        return(NULL)
      }

      X <- as.matrix(display_data[, qvars, drop = FALSE])
      ok <- stats::complete.cases(X)
      scores <- matrix(NA_real_, nrow = nrow(display_data), ncol = 2L)
      colnames(scores) <- c("PCA1", "PCA2")
      if (sum(ok) < 2L) {
        return(as.data.frame(scores))
      }

      pca <- stats::prcomp(X[ok, , drop = FALSE], center = TRUE, scale. = TRUE)
      scores[ok, seq_len(min(2L, ncol(pca$x)))] <- pca$x[, seq_len(min(2L, ncol(pca$x))), drop = FALSE]
      as.data.frame(scores)
    }


    cec_grid_display_choices <- function(display_data) {
      qvars <- cec_grid_quant_vars(display_data)
      pca_choices <- if (length(qvars) >= 2L) c("PCA1", "PCA2") else character(0)
      c(pca_choices, qvars)
    }


    cec_grid_display_vector <- function(display_data, var_name, pca_scores = NULL) {
      if (var_name %in% c("PCA1", "PCA2")) {
        if (is.null(pca_scores) || !(var_name %in% names(pca_scores))) {
          return(rep(NA_real_, nrow(display_data)))
        }
        return(pca_scores[[var_name]])
      }
      if (var_name %in% names(display_data)) {
        return(display_data[[var_name]])
      }
      rep(NA_real_, nrow(display_data))
    }


    cec_grid_relabel_partition_by_weight <- function(partition, best_obj = NULL) {
      partition <- as.integer(as.factor(partition))
      levs <- sort(unique(partition))
      params <- if (is.null(best_obj)) NULL else cec_grid_fit_params(best_obj)
      nu <- if (!is.null(params)) suppressWarnings(as.numeric(params$nu)) else numeric(0)

      weights <- if (length(nu) >= max(levs)) {
        nu[levs]
      } else {
        as.numeric(table(factor(partition, levels = levs))) / length(partition)
      }
      weights[!is.finite(weights)] <- 0

      ord <- order(weights, decreasing = TRUE, seq_along(weights))
      map <- setNames(seq_along(levs), levs[ord])
      relabeled <- as.integer(map[as.character(partition)])

      list(
        partition = relabeled,
        weights = weights[ord],
        original_cluster = levs[ord],
        cluster = seq_along(levs)
      )
    }


    cec_grid_selected_partition_info <- function(grid_obj, C, lambda) {
      result_obj <- cec_grid_result_for_C(grid_obj, C)
      selected <- cec_grid_best_obj_for_lambda(result_obj, lambda)
      partition <- extract_cec_partition(selected$obj)
      relabeled <- cec_grid_relabel_partition_by_weight(partition, selected$obj)
      selected$partition <- relabeled$partition
      selected$cluster_weights <- relabeled$weights
      selected$original_cluster <- relabeled$original_cluster
      selected$cluster <- relabeled$cluster
      selected
    }


    cec_grid_best_row <- function(grid_obj) {
      s <- cec_grid_summary_data(grid_obj)
      candidate <- which(isTRUE(s$keep) | s$keep)
      if (length(candidate) == 0L) {
        candidate <- which(is.finite(s$REO))
      }
      if (length(candidate) == 0L) {
        return(s[1L, , drop = FALSE])
      }
      s[candidate[which.min(s$REO[candidate])], , drop = FALSE]
    }


    cec_grid_summary_data <- function(grid_obj) {
      s <- grid_obj$summary
      if (is.data.frame(s)) {
        return(s)
      }

      if (!is.null(grid_obj$runs_by_C) && length(grid_obj$runs_by_C) > 0L) {
        C_grid <- grid_obj$C_grid
        if (is.null(C_grid) || length(C_grid) != length(grid_obj$runs_by_C)) {
          C_grid <- vapply(grid_obj$runs_by_C, function(x) x$meta$C, numeric(1))
        }
        summary_list <- vector("list", length(grid_obj$runs_by_C))
        for (i in seq_along(grid_obj$runs_by_C)) {
          summary_list[[i]] <- cec_grid_summary_for_single_C(
            grid_obj$runs_by_C[[i]],
            C_requested = C_grid[i]
          )
        }
        out <- do.call(rbind, summary_list)
        row.names(out) <- NULL
        return(out)
      }

      stop(
        "launch_cec_grid_shiny() expects an object returned by run_cec_grid(). ",
        "This object has no tabular summary and cannot be displayed as a C x lambda grid.",
        call. = FALSE
      )
    }


    cec_grid_result_for_C <- function(grid_obj, C) {
      idx <- which.min(abs(grid_obj$C_grid - C))
      grid_obj$runs_by_C[[idx]]
    }


    cec_grid_best_obj_for_lambda <- function(result_obj, lambda) {
      lambdas <- result_obj$best_parts$summary$lambda
      idx <- which.min(abs(lambdas - lambda))
      list(
        index = idx,
        lambda = lambdas[idx],
        obj = result_obj$best_parts$best[[idx]],
        summary = result_obj$best_parts$summary[idx, , drop = FALSE]
      )
    }


    cec_grid_partition_for_cell <- function(grid_obj, C, lambda) {
      result_obj <- cec_grid_result_for_C(grid_obj, C)
      selected <- cec_grid_best_obj_for_lambda(result_obj, lambda)
      extract_cec_partition(selected$obj)
    }


    cec_grid_reo_same <- function(x, y, tol = sqrt(.Machine$double.eps)) {
      if (!is.finite(x) && !is.finite(y)) {
        return(TRUE)
      }
      if (!is.finite(x) || !is.finite(y)) {
        return(FALSE)
      }
      abs(x - y) <= tol
    }


    cec_grid_partition_similarity <- function(x, y) {
      ceclust_fun("partition_similarity")(x, y, method = "match")
    }


    cec_grid_apply_thresholds <- function(
      grid_obj,
      stab_algo_threshold = NULL,
      rsi_threshold = NULL,
      sat_threshold = NULL
    ) {
      s <- cec_grid_summary_data(grid_obj)
      if (!is.data.frame(s) || nrow(s) == 0L) {
        stop("Empty C x lambda summary: nothing to display.", call. = FALSE)
      }
      if (is.null(stab_algo_threshold)) {
        stab_algo_threshold <- if (!is.null(grid_obj$meta$stab_algo_threshold)) grid_obj$meta$stab_algo_threshold else 0.8
      }
      if (is.null(rsi_threshold)) {
        rsi_threshold <- if (!is.null(grid_obj$meta$rsi_threshold)) grid_obj$meta$rsi_threshold else 0.8
      }
      if (is.null(sat_threshold)) {
        sat_threshold <- if (!is.null(grid_obj$meta$sat_threshold)) grid_obj$meta$sat_threshold else 0.1
      }
      keep_stability <- is.finite(s$stab_algo) &
        is.finite(s$rsi) &
        s$stab_algo >= stab_algo_threshold &
        s$rsi >= rsi_threshold
      exclude_saturation <- is.finite(s$sat) & s$sat > sat_threshold
      keep_all <- keep_stability & !exclude_saturation

      fail_reason <- rep("ok", nrow(s))
      fail_reason[!keep_stability] <- "stability_or_rsi"
      fail_reason[exclude_saturation] <- "saturation"
      fail_reason[!keep_stability & exclude_saturation] <- "stability_rsi_and_saturation"

      s$keep_stability <- keep_stability
      s$exclude_saturation <- exclude_saturation
      s$keep <- keep_all
      s$fail_reason <- fail_reason
      s$stab_algo_threshold <- stab_algo_threshold
      s$rsi_threshold <- rsi_threshold
      s$sat_threshold <- sat_threshold
      s
    }


    cec_grid_reo_range <- function(grid_summary) {
      kept_vals <- grid_summary$REO[grid_summary$keep & is.finite(grid_summary$REO)]
      vals <- if (length(kept_vals) > 0L) {
        kept_vals
      } else {
        grid_summary$REO[is.finite(grid_summary$REO)]
      }

      if (length(vals) == 0L) {
        return(c(0, 1))
      }

      out <- range(vals)
      if (!all(is.finite(out)) || diff(out) == 0) {
        center <- if (is.finite(out[1L])) out[1L] else 0
        pad <- max(1, abs(center) * 0.05)
        out <- c(center - pad, center + pad)
      }

      out
    }


    cec_grid_zone_ids <- function(
      grid_obj,
      grid_summary,
      sensibility = 0.95,
      reo_tol = sqrt(.Machine$double.eps)
    ) {
      s <- grid_summary
      C_vals <- sort(unique(s$C))
      lambda_vals <- sort(unique(s$lambda))
      zone_id <- rep(NA_integer_, nrow(s))
      zone_rep <- integer(0)
      next_zone <- 1L

      row_at <- function(C, lambda) {
        which(abs(s$C - C) <= sqrt(.Machine$double.eps) &
                abs(s$lambda - lambda) <= sqrt(.Machine$double.eps))[1L]
      }

      compatible <- function(row_idx, rep_idx) {
        if (is.na(row_idx) || is.na(rep_idx)) {
          return(FALSE)
        }
        if (!cec_grid_reo_same(s$REO[row_idx], s$REO[rep_idx], tol = reo_tol)) {
          return(FALSE)
        }
        p <- cec_grid_partition_for_cell(grid_obj, s$C[row_idx], s$lambda[row_idx])
        p_ref <- cec_grid_partition_for_cell(grid_obj, s$C[rep_idx], s$lambda[rep_idx])
        isTRUE(cec_grid_partition_similarity(p, p_ref) >= sensibility)
      }

      merge_zones <- function(target, source) {
        if (target == source) {
          return(invisible(NULL))
        }
        zone_id[zone_id == source] <<- target
        if (zone_rep[source] < zone_rep[target]) {
          zone_rep[target] <<- zone_rep[source]
        }
        zone_rep[source] <<- NA_integer_
        invisible(NULL)
      }

      for (C_value in C_vals) {
        for (lambda_value in lambda_vals) {
          idx <- row_at(C_value, lambda_value)
          if (is.na(idx) || !isTRUE(s$keep[idx])) {
            next
          }

          candidates <- integer(0)
          C_pos <- match(C_value, C_vals)
          lambda_pos <- match(lambda_value, lambda_vals)

          if (C_pos > 1L) {
            left_idx <- row_at(C_vals[C_pos - 1L], lambda_value)
            if (!is.na(left_idx) && isTRUE(s$keep[left_idx]) && !is.na(zone_id[left_idx])) {
              candidates <- c(candidates, zone_id[left_idx])
            }
          }
          if (lambda_pos > 1L) {
            lower_idx <- row_at(C_value, lambda_vals[lambda_pos - 1L])
            if (!is.na(lower_idx) && isTRUE(s$keep[lower_idx]) && !is.na(zone_id[lower_idx])) {
              candidates <- c(candidates, zone_id[lower_idx])
            }
          }

          candidates <- unique(candidates)
          candidates <- candidates[!is.na(candidates)]
          candidates <- candidates[vapply(
            candidates,
            function(z) compatible(idx, zone_rep[z]),
            logical(1)
          )]

          if (length(candidates) == 0L) {
            zone_id[idx] <- next_zone
            zone_rep[next_zone] <- idx
            next_zone <- next_zone + 1L
          } else {
            reps <- zone_rep[candidates]
            chosen <- candidates[which.min(reps)]
            zone_id[idx] <- chosen
            for (z in candidates) {
              merge_zones(chosen, z)
            }
            if (idx < zone_rep[chosen]) {
              zone_rep[chosen] <- idx
            }
          }
        }
      }

      zone_id
    }


    cec_grid_draw_zone_borders <- function(
      grid_summary,
      zone_id,
      c_bounds,
      l_bounds,
      col = "black",
      lwd = 2.2
    ) {
      s <- grid_summary
      row_at <- function(C, lambda) {
        which(abs(s$C - C) <= sqrt(.Machine$double.eps) &
                abs(s$lambda - lambda) <= sqrt(.Machine$double.eps))[1L]
      }
      zone_at <- function(C, lambda) {
        idx <- row_at(C, lambda)
        if (is.na(idx)) NA_integer_ else zone_id[idx]
      }

      C_vals <- sort(unique(s$C))
      lambda_vals <- sort(unique(s$lambda))

      for (idx in which(!is.na(zone_id))) {
        C_value <- s$C[idx]
        lambda_value <- s$lambda[idx]
        cb <- c_bounds[c_bounds$value == C_value, , drop = FALSE]
        lb <- l_bounds[l_bounds$value == lambda_value, , drop = FALSE]
        z <- zone_id[idx]
        C_pos <- match(C_value, C_vals)
        lambda_pos <- match(lambda_value, lambda_vals)

        left_zone <- if (C_pos > 1L) zone_at(C_vals[C_pos - 1L], lambda_value) else NA_integer_
        right_zone <- if (C_pos < length(C_vals)) zone_at(C_vals[C_pos + 1L], lambda_value) else NA_integer_
        lower_zone <- if (lambda_pos > 1L) zone_at(C_value, lambda_vals[lambda_pos - 1L]) else NA_integer_
        upper_zone <- if (lambda_pos < length(lambda_vals)) zone_at(C_value, lambda_vals[lambda_pos + 1L]) else NA_integer_

        if (is.na(left_zone) || left_zone != z) {
          graphics::segments(cb$lower, lb$lower, cb$lower, lb$upper, col = col, lwd = lwd)
        }
        if (is.na(right_zone) || right_zone != z) {
          graphics::segments(cb$upper, lb$lower, cb$upper, lb$upper, col = col, lwd = lwd)
        }
        if (is.na(lower_zone) || lower_zone != z) {
          graphics::segments(cb$lower, lb$lower, cb$upper, lb$lower, col = col, lwd = lwd)
        }
        if (is.na(upper_zone) || upper_zone != z) {
          graphics::segments(cb$lower, lb$upper, cb$upper, lb$upper, col = col, lwd = lwd)
        }
      }

      invisible(zone_id)
    }


    cec_grid_partition_bands_1d <- function(Z, partition) {
      ord <- order(Z)
      z_ord <- Z[ord]
      p_ord <- as.integer(as.factor(partition[ord]))
      r <- rle(p_ord)
      ends <- cumsum(r$lengths)
      starts <- c(1L, head(ends, -1L) + 1L)

      lower_point <- function(i) {
        if (i <= 1L) {
          return(min(z_ord))
        }
        (z_ord[i - 1L] + z_ord[i]) / 2
      }
      upper_point <- function(i) {
        if (i >= length(z_ord)) {
          return(max(z_ord))
        }
        (z_ord[i] + z_ord[i + 1L]) / 2
      }

      data.frame(
        xmin = vapply(starts, lower_point, numeric(1)),
        xmax = vapply(ends, upper_point, numeric(1)),
        cluster = r$values,
        stringsAsFactors = FALSE
      )
    }


    plot_cec_grid_reo_map <- function(
      grid_obj,
      selected_C = NULL,
      selected_lambda = NULL,
      stab_algo_threshold = grid_obj$meta$stab_algo_threshold,
      rsi_threshold = grid_obj$meta$rsi_threshold,
      sat_threshold = grid_obj$meta$sat_threshold,
      sensibility = 0.95,
      grid_summary = NULL,
      main = "REO on the C x lambda grid",
      reo_col = cec_grid_color_palette(100),
      bad_fill = grDevices::adjustcolor("grey70", alpha.f = 0.72),
      bad_hatch_col = grDevices::adjustcolor("grey25", alpha.f = 0.75),
      draw_zones = TRUE,
      show_legend = FALSE
    ) {
      s <- if (is.null(grid_summary)) {
        cec_grid_apply_thresholds(
          grid_obj,
          stab_algo_threshold = stab_algo_threshold,
          rsi_threshold = rsi_threshold,
          sat_threshold = sat_threshold
        )
      } else {
        grid_summary
      }
      c_bounds <- cec_grid_cell_bounds(s$C)
      l_bounds <- cec_grid_cell_bounds(s$lambda)
      reo_range <- cec_grid_reo_range(s)

      oldpar <- graphics::par(c("mar", "xpd"))
      on.exit(graphics::par(mar = oldpar$mar, xpd = oldpar$xpd), add = TRUE)
      graphics::par(mar = c(4.2, 4.3, 2.7, 1.1))

      graphics::plot(
        NA,
        NA,
        xlim = range(c_bounds$lower, c_bounds$upper),
        ylim = range(l_bounds$lower, l_bounds$upper),
        xlab = "C",
        ylab = expression(lambda),
        main = main,
        axes = FALSE
      )

      for (i in seq_len(nrow(s))) {
        cb <- c_bounds[c_bounds$value == s$C[i], , drop = FALSE]
        lb <- l_bounds[l_bounds$value == s$lambda[i], , drop = FALSE]
        if (isTRUE(s$keep[i]) && is.finite(s$REO[i])) {
          col_idx <- 1L + floor(
            (length(reo_col) - 1L) * (s$REO[i] - reo_range[1L]) / diff(reo_range)
          )
          col_idx <- max(1L, min(length(reo_col), col_idx))
          cell_col <- reo_col[col_idx]
        } else {
          cell_col <- bad_fill
        }

        graphics::rect(
          xleft = cb$lower,
          ybottom = lb$lower,
          xright = cb$upper,
          ytop = lb$upper,
          col = cell_col,
          border = "white",
          lwd = 0.7
        )

        if (!isTRUE(s$keep[i])) {
          if (!isTRUE(s$keep_stability[i])) {
            graphics::rect(
              xleft = cb$lower,
              ybottom = lb$lower,
              xright = cb$upper,
              ytop = lb$upper,
              density = 10,
              angle = 45,
              col = bad_hatch_col,
              border = NA
            )
          }
          if (isTRUE(s$exclude_saturation[i])) {
            graphics::rect(
              xleft = cb$lower,
              ybottom = lb$lower,
              xright = cb$upper,
              ytop = lb$upper,
              density = 10,
              angle = -45,
              col = bad_hatch_col,
              border = NA
            )
          }
        }
      }

      if (isTRUE(draw_zones)) {
        zone_id <- cec_grid_zone_ids(
          grid_obj = grid_obj,
          grid_summary = s,
          sensibility = sensibility
        )
        cec_grid_draw_zone_borders(
          grid_summary = s,
          zone_id = zone_id,
          c_bounds = c_bounds,
          l_bounds = l_bounds
        )
      }

      if (!is.null(selected_C) && !is.null(selected_lambda)) {
        cb <- c_bounds[which.min(abs(c_bounds$value - selected_C)), , drop = FALSE]
        lb <- l_bounds[which.min(abs(l_bounds$value - selected_lambda)), , drop = FALSE]
        graphics::rect(
          xleft = cb$lower,
          ybottom = lb$lower,
          xright = cb$upper,
          ytop = lb$upper,
          border = "black",
          lwd = 2.2
        )
      }

      graphics::axis(1, at = c_bounds$value, labels = c_bounds$value)
      graphics::axis(2, at = l_bounds$value, labels = l_bounds$value)
      graphics::box()

      if (isTRUE(show_legend)) {
        plot_cec_grid_reo_legend(
          grid_obj = grid_obj,
          grid_summary = s,
          reo_col = reo_col,
          bad_fill = bad_fill,
          bad_hatch_col = bad_hatch_col
        )
      }

      invisible(s)
    }


    plot_cec_grid_reo_legend <- function(
      grid_obj,
      grid_summary = NULL,
      reo_col = cec_grid_color_palette(100),
      bad_fill = grDevices::adjustcolor("grey70", alpha.f = 0.72),
      bad_hatch_col = grDevices::adjustcolor("grey25", alpha.f = 0.75)
    ) {
      s <- if (is.null(grid_summary)) cec_grid_summary_data(grid_obj) else grid_summary
      reo_range <- cec_grid_reo_range(s)

      oldpar <- graphics::par(c("mar", "xpd"))
      on.exit(graphics::par(mar = oldpar$mar, xpd = oldpar$xpd), add = TRUE)
      graphics::par(mar = c(0.5, 0.5, 0.5, 0.5), xpd = NA)
      graphics::plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")

      n_col <- length(reo_col)
      x0 <- 0.06
      x1 <- 0.58
      y0 <- 0.62
      y1 <- 0.82
      xs <- seq(x0, x1, length.out = n_col + 1L)
      for (i in seq_len(n_col)) {
        graphics::rect(xs[i], y0, xs[i + 1L], y1, col = reo_col[i], border = NA)
      }
      graphics::rect(x0, y0, x1, y1, border = "grey25", lwd = 0.8)
      graphics::text(x0, y0 - 0.10, labels = signif(reo_range[1L], 4), adj = c(0, 0.5), cex = 0.82)
      graphics::text(x1, y0 - 0.10, labels = signif(reo_range[2L], 4), adj = c(1, 0.5), cex = 0.82)
      graphics::text(x0, y1 + 0.11, labels = "REO (kept cells)", adj = c(0, 0.5), cex = 0.9)

      graphics::rect(0.66, 0.67, 0.74, 0.82, col = bad_fill, border = "grey25")
      graphics::rect(0.66, 0.67, 0.74, 0.82, density = 10, angle = 45, col = bad_hatch_col, border = NA)
      graphics::text(0.78, 0.745, labels = "rejected: stability/RSI", adj = c(0, 0.5), cex = 0.82)

      graphics::rect(0.66, 0.40, 0.74, 0.55, col = bad_fill, border = "grey25")
      graphics::rect(0.66, 0.40, 0.74, 0.55, density = 10, angle = -45, col = bad_hatch_col, border = NA)
      graphics::text(0.78, 0.475, labels = "rejected: saturation", adj = c(0, 0.5), cex = 0.82)

      graphics::rect(0.66, 0.13, 0.74, 0.28, col = bad_fill, border = "grey25")
      graphics::rect(0.66, 0.13, 0.74, 0.28, density = 10, angle = 45, col = bad_hatch_col, border = NA)
      graphics::rect(0.66, 0.13, 0.74, 0.28, density = 10, angle = -45, col = bad_hatch_col, border = NA)
      graphics::text(0.78, 0.205, labels = "rejected: both", adj = c(0, 0.5), cex = 0.82)

      invisible(s)
    }


    plot_cec_grid_selected_partition_1d <- function(
      grid_obj,
      C,
      lambda,
      hist_breaks = NULL,
      histogram_col = grDevices::adjustcolor("grey75", alpha.f = 0.45),
      density_col = "#222222",
      main = NULL
    ) {
      result_obj <- cec_grid_result_for_C(grid_obj, C)
      selected <- cec_grid_best_obj_for_lambda(result_obj, lambda)
      obj <- selected$obj
      Z <- result_obj$Z
      if (is.data.frame(Z) || is.matrix(Z)) {
        Z <- as.numeric(as.data.frame(Z)[[1L]])
      }
      partition <- extract_cec_partition(obj)
      k <- length(unique(partition))
      palette <- cec_grid_cluster_palette(k)
      bands <- cec_grid_partition_bands_1d(Z, partition)
      density_df <- result_obj$density
      has_density <- is.data.frame(density_df) &&
        all(c("x", "y") %in% names(density_df)) &&
        any(is.finite(density_df$x)) &&
        any(is.finite(density_df$y))

      if (is.null(hist_breaks)) {
        hist_breaks <- max(6L, min(80L, ceiling(sqrt(length(Z)))))
      }
      hist_obj <- graphics::hist(Z, breaks = hist_breaks, plot = FALSE)
      ymax_values <- hist_obj$density
      if (isTRUE(has_density)) {
        ymax_values <- c(ymax_values, density_df$y)
      }
      ymax <- max(ymax_values, na.rm = TRUE)
      if (!is.finite(ymax) || ymax <= 0) {
        ymax <- 1
      }
      xlim_values <- Z
      if (isTRUE(has_density)) {
        xlim_values <- c(xlim_values, density_df$x)
      }
      xlim <- range(xlim_values, finite = TRUE)

      if (is.null(main)) {
        main <- paste0(
          "Best partition for C = ", format(C),
          ", lambda = ", format(selected$lambda)
        )
      }

      graphics::plot(
        NA,
        NA,
        xlim = xlim,
        ylim = c(0, 1.08 * ymax),
        xlab = "Z",
        ylab = "density",
        main = main,
        type = "n"
      )
      for (i in seq_len(nrow(bands))) {
        graphics::rect(
          bands$xmin[i],
          0,
          bands$xmax[i],
          1.08 * ymax,
          col = grDevices::adjustcolor(palette[bands$cluster[i]], alpha.f = 0.18),
          border = NA
        )
      }
      graphics::plot(hist_obj, freq = FALSE, add = TRUE, col = histogram_col, border = "white")
      if (isTRUE(has_density)) {
        graphics::lines(density_df$x, density_df$y, col = density_col, lwd = 2)
      }
      graphics::rug(Z, col = grDevices::adjustcolor("grey20", alpha.f = 0.55))
      invisible(selected)
    }


    plot_cec_grid_selected_lambda_path_1d <- function(
      grid_obj,
      C,
      selected_lambda = NULL,
      main = NULL
    ) {
      result_obj <- cec_grid_result_for_C(grid_obj, C)
      s <- cec_grid_summary_for_single_C(result_obj, C_requested = C)
      if (is.null(main)) {
        main <- paste0("Diagnostics for C = ", format(C))
      }

      y <- s$REO
      yr <- range(y, finite = TRUE)
      if (!all(is.finite(yr)) || diff(yr) == 0) {
        yr <- c(0, 1)
      }

      graphics::plot(
        s$lambda,
        y,
        type = "b",
        pch = ifelse(s$keep, 19, 1),
        col = ifelse(s$keep, "#2166AC", "grey45"),
        xlab = expression(lambda),
        ylab = "REO",
        main = main,
        ylim = yr
      )
      graphics::grid(col = "grey90")
      if (!is.null(selected_lambda)) {
        graphics::abline(v = selected_lambda, col = "#B2182B", lwd = 2)
      }
      invisible(s)
    }


    plot_cec_grid_selected_scatter <- function(
      grid_obj,
      C,
      lambda,
      x_var,
      y_var,
      shape_var = "None",
      cluster_selection = "All",
      point_cex = 1.05
    ) {
      display_data <- cec_grid_display_data(grid_obj)
      if (nrow(display_data) == 0L) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No display data available.")
        return(invisible(NULL))
      }

      pca_scores <- cec_grid_pca_scores(display_data)
      x <- cec_grid_display_vector(display_data, x_var, pca_scores)
      y <- cec_grid_display_vector(display_data, y_var, pca_scores)
      selected <- cec_grid_selected_partition_info(grid_obj, C, lambda)
      partition <- selected$partition
      if (length(partition) != nrow(display_data)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Partition and display data lengths differ.")
        return(invisible(NULL))
      }

      k <- max(partition, na.rm = TRUE)
      pal <- cec_grid_cluster_palette(k)
      keep <- rep(TRUE, length(partition))
      if (!identical(cluster_selection, "All")) {
        keep <- partition %in% as.integer(cluster_selection)
      }

      shape <- rep(21, length(partition))
      shape_legend <- NULL
      if (!is.null(shape_var) && !identical(shape_var, "None") && shape_var %in% names(display_data)) {
        f <- factor(display_data[[shape_var]])
        shape_values <- c(21:25, 0:14)
        shape_map <- setNames(shape_values[((seq_along(levels(f)) - 1L) %% length(shape_values)) + 1L], levels(f))
        shape <- unname(shape_map[as.character(f)])
        shape_legend <- list(levels = levels(f), pch = unname(shape_map))
      }

      ok <- is.finite(x) & is.finite(y) & keep
      plot_xlim <- range(x[is.finite(x)], finite = TRUE)
      plot_ylim <- range(y[is.finite(y)], finite = TRUE)
      if (!all(is.finite(plot_xlim))) {
        plot_xlim <- range(x[ok], finite = TRUE)
      }
      if (!all(is.finite(plot_ylim))) {
        plot_ylim <- range(y[ok], finite = TRUE)
      }
      graphics::plot(
        x[ok],
        y[ok],
        xlab = x_var,
        ylab = y_var,
        main = paste0("Partition for C = ", format(C), ", lambda = ", format(selected$lambda)),
        type = "n",
        xlim = plot_xlim,
        ylim = plot_ylim
      )
      graphics::grid(col = "grey90")
      graphics::points(
        x[ok],
        y[ok],
        pch = shape[ok],
        bg = pal[partition[ok]],
        col = "grey20",
        cex = point_cex
      )

      cluster_levels <- sort(unique(partition[ok]))
      if (length(cluster_levels) > 0L) {
        graphics::legend(
          "topright",
          legend = paste0("cluster ", cluster_levels),
          pt.bg = pal[cluster_levels],
          pch = 21,
          bty = "n",
          cex = 0.8
        )
      }
      if (!is.null(shape_legend) && length(shape_legend$levels) > 0L) {
        graphics::legend(
          "bottomright",
          legend = shape_legend$levels,
          pch = shape_legend$pch,
          col = "grey20",
          bty = "n",
          cex = 0.75,
          title = shape_var
        )
      }

      invisible(selected)
    }


    plot_cec_grid_qual_barplot <- function(
      grid_obj,
      C,
      lambda,
      shape_var = "None",
      cluster_selection = "All"
    ) {
      display_data <- cec_grid_display_data(grid_obj)
      if (is.null(shape_var) || identical(shape_var, "None") || !(shape_var %in% names(display_data))) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Choose a qualitative variable for the barplot.")
        return(invisible(NULL))
      }

      selected <- cec_grid_selected_partition_info(grid_obj, C, lambda)
      partition <- selected$partition
      if (length(partition) != nrow(display_data)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Partition and display data lengths differ.")
        return(invisible(NULL))
      }

      clusters <- sort(unique(partition))
      if (!identical(cluster_selection, "All")) {
        clusters <- intersect(clusters, as.integer(cluster_selection))
      }
      keep <- partition %in% clusters
      if (!any(keep)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No point in selected cluster(s).")
        return(invisible(NULL))
      }

      qual <- factor(display_data[[shape_var]])
      tab <- table(
        Level = qual[keep],
        Cluster = factor(partition[keep], levels = clusters)
      )
      pal <- cec_grid_cluster_palette(max(partition, na.rm = TRUE))[clusters]

      if (nrow(tab) == 1L) {
        graphics::barplot(
          tab[1L, ],
          col = pal,
          border = "white",
          las = 2,
          ylab = "count",
          main = paste0(shape_var, " by selected clusters")
        )
      } else {
        graphics::barplot(
          t(tab),
          col = pal,
          border = "white",
          las = 2,
          ylab = "count",
          main = paste0(shape_var, " by selected clusters")
        )
      }

      invisible(tab)
    }


    cec_grid_icl_partition_info <- function(icl_obj) {
      partition <- if (!is.null(icl_obj$partition)) {
        icl_obj$partition
      } else {
        icl_obj$icl_fit$classification
      }
      best_obj <- list(params = list(nu = icl_obj$cluster_weights))
      relabeled <- cec_grid_relabel_partition_by_weight(partition, best_obj)
      list(
        partition = relabeled$partition,
        cluster_weights = relabeled$weights,
        original_cluster = relabeled$original_cluster,
        cluster = relabeled$cluster
      )
    }


    plot_cec_icl_partition_1d <- function(
      icl_obj,
      cluster_selection = "All",
      hist_breaks = NULL,
      histogram_col = grDevices::adjustcolor("grey75", alpha.f = 0.45),
      main = NULL
    ) {
      display_data <- cec_grid_display_data(icl_obj)
      qvars <- cec_grid_quant_vars(display_data)
      if (length(qvars) == 0L) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No quantitative variable available.")
        return(invisible(NULL))
      }

      Z <- as.numeric(display_data[[qvars[1L]]])
      info <- cec_grid_icl_partition_info(icl_obj)
      partition <- info$partition
      if (length(partition) != length(Z)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Partition and display data lengths differ.")
        return(invisible(NULL))
      }

      keep <- is.finite(Z)
      if (!identical(cluster_selection, "All")) {
        keep <- keep & partition %in% as.integer(cluster_selection)
      }
      if (!any(keep)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No point in selected cluster.")
        return(invisible(NULL))
      }

      k <- max(partition, na.rm = TRUE)
      pal <- cec_grid_cluster_palette(k)
      bands <- cec_grid_partition_bands_1d(Z, partition)
      density_df <- icl_obj$density
      has_density <- is.data.frame(density_df) &&
        all(c("x", "y") %in% names(density_df)) &&
        any(is.finite(density_df$x)) &&
        any(is.finite(density_df$y))
      if (is.null(hist_breaks)) {
        hist_breaks <- max(6L, min(80L, ceiling(sqrt(length(Z)))))
      }
      hist_obj <- graphics::hist(Z[keep], breaks = hist_breaks, plot = FALSE)
      ymax_values <- hist_obj$density
      if (isTRUE(has_density)) {
        ymax_values <- c(ymax_values, density_df$y)
      }
      ymax <- max(ymax_values, na.rm = TRUE)
      if (!is.finite(ymax) || ymax <= 0) {
        ymax <- 1
      }
      if (is.null(main)) {
        main <- paste0(
          "ICL partition: ",
          icl_obj$meta$modelName,
          ", G = ",
          icl_obj$meta$best_G
        )
      }

      graphics::plot(
        NA,
        NA,
        xlim = range(Z, finite = TRUE),
        ylim = c(0, 1.08 * ymax),
        xlab = qvars[1L],
        ylab = "density",
        main = main,
        type = "n"
      )
      for (i in seq_len(nrow(bands))) {
        graphics::rect(
          bands$xmin[i],
          0,
          bands$xmax[i],
          1.08 * ymax,
          col = grDevices::adjustcolor(pal[bands$cluster[i]], alpha.f = 0.18),
          border = NA
        )
      }
      graphics::plot(hist_obj, freq = FALSE, add = TRUE, col = histogram_col, border = "white")
      if (isTRUE(has_density)) {
        graphics::lines(density_df$x, density_df$y, col = "#222222", lwd = 2)
      }
      graphics::rug(Z[keep], col = grDevices::adjustcolor("grey20", alpha.f = 0.55))
      invisible(info)
    }


    plot_cec_icl_selected_scatter <- function(
      icl_obj,
      x_var,
      y_var,
      shape_var = "None",
      cluster_selection = "All",
      point_cex = 1.05
    ) {
      display_data <- cec_grid_display_data(icl_obj)
      pca_scores <- cec_grid_pca_scores(display_data)
      x <- cec_grid_display_vector(display_data, x_var, pca_scores)
      y <- cec_grid_display_vector(display_data, y_var, pca_scores)
      info <- cec_grid_icl_partition_info(icl_obj)
      partition <- info$partition
      if (length(partition) != nrow(display_data)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Partition and display data lengths differ.")
        return(invisible(NULL))
      }

      k <- max(partition, na.rm = TRUE)
      pal <- cec_grid_cluster_palette(k)
      keep <- rep(TRUE, length(partition))
      if (!identical(cluster_selection, "All")) {
        keep <- partition %in% as.integer(cluster_selection)
      }

      shape <- rep(21, length(partition))
      shape_legend <- NULL
      if (!is.null(shape_var) && !identical(shape_var, "None") && shape_var %in% names(display_data)) {
        f <- factor(display_data[[shape_var]])
        shape_values <- c(21:25, 0:14)
        shape_map <- setNames(shape_values[((seq_along(levels(f)) - 1L) %% length(shape_values)) + 1L], levels(f))
        shape <- unname(shape_map[as.character(f)])
        shape_legend <- list(levels = levels(f), pch = unname(shape_map))
      }

      ok <- is.finite(x) & is.finite(y) & keep
      plot_xlim <- range(x[is.finite(x)], finite = TRUE)
      plot_ylim <- range(y[is.finite(y)], finite = TRUE)
      graphics::plot(
        x[ok],
        y[ok],
        xlab = x_var,
        ylab = y_var,
        main = paste0("ICL partition: ", icl_obj$meta$modelName, ", G = ", icl_obj$meta$best_G),
        type = "n",
        xlim = plot_xlim,
        ylim = plot_ylim
      )
      graphics::grid(col = "grey90")
      graphics::points(
        x[ok],
        y[ok],
        pch = shape[ok],
        bg = pal[partition[ok]],
        col = "grey20",
        cex = point_cex
      )

      cluster_levels <- sort(unique(partition[ok]))
      if (length(cluster_levels) > 0L) {
        graphics::legend(
          "topright",
          legend = paste0("cluster ", cluster_levels),
          pt.bg = pal[cluster_levels],
          pch = 21,
          bty = "n",
          cex = 0.8
        )
      }
      if (!is.null(shape_legend) && length(shape_legend$levels) > 0L) {
        graphics::legend(
          "bottomright",
          legend = shape_legend$levels,
          pch = shape_legend$pch,
          col = "grey20",
          bty = "n",
          cex = 0.75,
          title = shape_var
        )
      }

      invisible(info)
    }


    plot_cec_icl_qual_barplot <- function(
      icl_obj,
      shape_var = "None",
      cluster_selection = "All"
    ) {
      display_data <- cec_grid_display_data(icl_obj)
      if (is.null(shape_var) || identical(shape_var, "None") || !(shape_var %in% names(display_data))) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Choose a qualitative variable for the barplot.")
        return(invisible(NULL))
      }

      info <- cec_grid_icl_partition_info(icl_obj)
      partition <- info$partition
      clusters <- sort(unique(partition))
      if (!identical(cluster_selection, "All")) {
        clusters <- intersect(clusters, as.integer(cluster_selection))
      }
      keep <- partition %in% clusters
      if (!any(keep)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No point in selected cluster(s).")
        return(invisible(NULL))
      }

      qual <- factor(display_data[[shape_var]])
      tab <- table(
        Level = qual[keep],
        Cluster = factor(partition[keep], levels = clusters)
      )
      pal <- cec_grid_cluster_palette(max(partition, na.rm = TRUE))[clusters]

      if (nrow(tab) == 1L) {
        graphics::barplot(
          tab[1L, ],
          col = pal,
          border = "white",
          las = 2,
          ylab = "count",
          main = paste0(shape_var, " by selected clusters")
        )
      } else {
        graphics::barplot(
          t(tab),
          col = pal,
          border = "white",
          las = 2,
          ylab = "count",
          main = paste0(shape_var, " by selected clusters")
        )
      }

      invisible(tab)
    }


    launch_cec_icl_shiny <- function(
      icl_obj,
      launch.browser = TRUE
    ) {
      if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Package 'shiny' is required for launch_cec_icl_shiny().", call. = FALSE)
      }

      display_data <- cec_grid_display_data(icl_obj)
      display_choices <- cec_grid_display_choices(display_data)
      qual_choices <- c("None", cec_grid_qual_vars(display_data))
      default_x <- if ("PCA1" %in% display_choices) "PCA1" else display_choices[1L]
      default_y <- if ("PCA2" %in% display_choices) {
        "PCA2"
      } else if (length(display_choices) >= 2L) {
        display_choices[2L]
      } else {
        display_choices[1L]
      }
      info <- cec_grid_icl_partition_info(icl_obj)
      cluster_choices <- c("All", as.character(info$cluster))

      ui <- shiny::fluidPage(
        shiny::titlePanel(paste0("ICL - ", icl_obj$dataset_name)),
        shiny::fluidRow(
          shiny::column(
            width = 3,
            shiny::wellPanel(
              shiny::tags$b("Selected ICL model"),
              shiny::tableOutput("icl_table")
            ),
            shiny::selectInput("cluster_select", "cluster", choices = cluster_choices, selected = "All")
          ),
          shiny::column(
            width = 9,
            shiny::uiOutput("display_controls"),
            shiny::uiOutput("detail_plots")
          )
        )
      )

      server <- function(input, output, session) {
        output$icl_table <- shiny::renderTable({
          data.frame(
            modelName = icl_obj$meta$modelName,
            G = icl_obj$meta$best_G,
            icl = icl_obj$meta$icl,
            bic = icl_obj$meta$bic,
            loglik = icl_obj$meta$loglik,
            stringsAsFactors = FALSE
          )
        }, digits = 4)

        output$display_controls <- shiny::renderUI({
          if (length(display_choices) < 2L) {
            return(NULL)
          }
          shiny::wellPanel(
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::selectInput("x_var", "x", choices = display_choices, selected = default_x)
              ),
              shiny::column(
                width = 4,
                shiny::selectInput("y_var", "y", choices = display_choices, selected = default_y)
              ),
              shiny::column(
                width = 4,
                shiny::selectInput("shape_var", "forme", choices = qual_choices, selected = qual_choices[1L])
              )
            )
          )
        })

        output$detail_plots <- shiny::renderUI({
          if (length(display_choices) >= 2L) {
            shiny::tagList(
              shiny::plotOutput("scatter_plot", height = "470px"),
              shiny::plotOutput("qual_barplot", height = "260px")
            )
          } else {
            shiny::plotOutput("partition_plot", height = "520px")
          }
        })

        output$partition_plot <- shiny::renderPlot({
          cluster_select <- if (is.null(input$cluster_select)) "All" else input$cluster_select
          plot_cec_icl_partition_1d(
            icl_obj,
            cluster_selection = cluster_select
          )
        })

        output$scatter_plot <- shiny::renderPlot({
          shiny::req(length(display_choices) >= 2L)
          x_var <- if (is.null(input$x_var)) default_x else input$x_var
          y_var <- if (is.null(input$y_var)) default_y else input$y_var
          shape_var <- if (is.null(input$shape_var)) "None" else input$shape_var
          cluster_select <- if (is.null(input$cluster_select)) "All" else input$cluster_select
          plot_cec_icl_selected_scatter(
            icl_obj,
            x_var = x_var,
            y_var = y_var,
            shape_var = shape_var,
            cluster_selection = cluster_select
          )
        })

        output$qual_barplot <- shiny::renderPlot({
          shiny::req(length(display_choices) >= 2L)
          shape_var <- if (is.null(input$shape_var)) "None" else input$shape_var
          cluster_select <- if (is.null(input$cluster_select)) "All" else input$cluster_select
          plot_cec_icl_qual_barplot(
            icl_obj,
            shape_var = shape_var,
            cluster_selection = cluster_select
          )
        })
      }

      shiny::runApp(
        shiny::shinyApp(ui = ui, server = server),
        launch.browser = launch.browser
      )
    }


    launch_cec_grid_shiny <- function(
      grid_obj,
      launch.browser = TRUE
    ) {
      if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Package 'shiny' is required for launch_cec_grid_shiny().")
      }

      default_summary <- cec_grid_apply_thresholds(grid_obj)
      default_row <- cec_grid_best_row(list(summary = default_summary))
      default_stab_algo_threshold <- if (!is.null(grid_obj$meta$stab_algo_threshold)) grid_obj$meta$stab_algo_threshold else 0.8
      default_rsi_threshold <- if (!is.null(grid_obj$meta$rsi_threshold)) grid_obj$meta$rsi_threshold else 0.8
      default_sat_threshold <- if (!is.null(grid_obj$meta$sat_threshold)) grid_obj$meta$sat_threshold else 0.1
      display_data <- cec_grid_display_data(grid_obj)
      display_choices <- cec_grid_display_choices(display_data)
      qual_choices <- c("None", cec_grid_qual_vars(display_data))
      default_x <- if ("PCA1" %in% display_choices) "PCA1" else display_choices[1L]
      default_y <- if ("PCA2" %in% display_choices) {
        "PCA2"
      } else if (length(display_choices) >= 2L) {
        display_choices[2L]
      } else {
        display_choices[1L]
      }
      selected <- shiny::reactiveVal(list(
        C = default_row$C[1L],
        lambda = default_row$lambda[1L]
      ))
      params <- shiny::reactiveVal(list(
        stab_algo_threshold = default_stab_algo_threshold,
        rsi_threshold = default_rsi_threshold,
        sat_threshold = default_sat_threshold,
        sensibility = 0.95
      ))

      ui <- shiny::fluidPage(
        shiny::titlePanel("CEC grid tests"),
        shiny::fluidRow(
          shiny::column(
            width = 4,
            shiny::wellPanel(
              shiny::sliderInput(
                "stab_algo_threshold",
                "stability",
                min = 0,
                max = 1,
                value = default_stab_algo_threshold,
                step = 0.01
              ),
              shiny::sliderInput(
                "rsi_threshold",
                "rsi",
                min = 0,
                max = 1,
                value = default_rsi_threshold,
                step = 0.01
              ),
              shiny::sliderInput(
                "sat_threshold",
                "saturation",
                min = 0,
                max = 1,
                value = default_sat_threshold,
                step = 0.01
              ),
              shiny::sliderInput(
                "sensibility",
                "sensibility",
                min = 0,
                max = 1,
                value = 0.95,
                step = 0.01
              ),
              shiny::actionButton("refresh_grid", "Actualiser")
            ),
            shiny::plotOutput("grid_plot", height = "470px", click = "grid_click"),
            shiny::plotOutput("legend_plot", height = "125px"),
            shiny::tableOutput("selected_table")
          ),
          shiny::column(
            width = 8,
            shiny::uiOutput("display_controls"),
            shiny::uiOutput("detail_plots"),
            shiny::plotOutput("lambda_plot", height = "300px")
          )
        )
      )

      server <- function(input, output, session) {
        current_summary <- shiny::reactive({
          p <- params()
          cec_grid_apply_thresholds(
            grid_obj,
            stab_algo_threshold = p$stab_algo_threshold,
            rsi_threshold = p$rsi_threshold,
            sat_threshold = p$sat_threshold
          )
        })

        shiny::observeEvent(input$refresh_grid, {
          params(list(
            stab_algo_threshold = input$stab_algo_threshold,
            rsi_threshold = input$rsi_threshold,
            sat_threshold = input$sat_threshold,
            sensibility = input$sensibility
          ))
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$grid_click, {
          s <- current_summary()
          idx <- cec_grid_row_nearest(
            s,
            C = input$grid_click$x,
            lambda = input$grid_click$y
          )
          row <- s[idx, , drop = FALSE]
          selected(list(C = row$C[1L], lambda = row$lambda[1L]))
        })

        cluster_choices <- shiny::reactive({
          sel <- selected()
          out <- tryCatch(
            {
              info <- cec_grid_selected_partition_info(grid_obj, sel$C, sel$lambda)
              c("All", as.character(info$cluster))
            },
            error = function(e) "All"
          )
          out
        })

        output$display_controls <- shiny::renderUI({
          if (length(display_choices) < 2L) {
            return(NULL)
          }
          shiny::wellPanel(
            shiny::fluidRow(
              shiny::column(
                width = 3,
                shiny::selectInput("x_var", "x", choices = display_choices, selected = default_x)
              ),
              shiny::column(
                width = 3,
                shiny::selectInput("y_var", "y", choices = display_choices, selected = default_y)
              ),
              shiny::column(
                width = 3,
                shiny::selectInput("shape_var", "forme", choices = qual_choices, selected = qual_choices[1L])
              ),
              shiny::column(
                width = 3,
                shiny::selectInput("cluster_select", "cluster", choices = "All", selected = "All")
              )
            )
          )
        })

        shiny::observeEvent(cluster_choices(), {
          if (length(display_choices) < 2L) {
            return(NULL)
          }
          choices <- cluster_choices()
          shiny::updateSelectInput(
            session,
            "cluster_select",
            choices = choices,
            selected = "All"
          )
        }, ignoreInit = FALSE)

        output$grid_plot <- shiny::renderPlot({
          sel <- selected()
          p <- params()
          plot_cec_grid_reo_map(
            grid_obj,
            selected_C = sel$C,
            selected_lambda = sel$lambda,
            grid_summary = current_summary(),
            sensibility = p$sensibility
          )
        })

        output$legend_plot <- shiny::renderPlot({
          plot_cec_grid_reo_legend(
            grid_obj,
            grid_summary = current_summary()
          )
        })

        output$detail_plots <- shiny::renderUI({
          if (length(display_choices) >= 2L) {
            shiny::tagList(
              shiny::plotOutput("scatter_plot", height = "410px"),
              shiny::plotOutput("qual_barplot", height = "250px")
            )
          } else {
            shiny::plotOutput("partition_plot", height = "430px")
          }
        })

        output$partition_plot <- shiny::renderPlot({
          sel <- selected()
          plot_cec_grid_selected_partition_1d(
            grid_obj,
            C = sel$C,
            lambda = sel$lambda
          )
        })

        output$scatter_plot <- shiny::renderPlot({
          shiny::req(length(display_choices) >= 2L)
          sel <- selected()
          x_var <- if (is.null(input$x_var)) default_x else input$x_var
          y_var <- if (is.null(input$y_var)) default_y else input$y_var
          shape_var <- if (is.null(input$shape_var)) "None" else input$shape_var
          cluster_select <- if (is.null(input$cluster_select)) "All" else input$cluster_select
          plot_cec_grid_selected_scatter(
            grid_obj,
            C = sel$C,
            lambda = sel$lambda,
            x_var = x_var,
            y_var = y_var,
            shape_var = shape_var,
            cluster_selection = cluster_select
          )
        })

        output$qual_barplot <- shiny::renderPlot({
          shiny::req(length(display_choices) >= 2L)
          sel <- selected()
          shape_var <- if (is.null(input$shape_var)) "None" else input$shape_var
          cluster_select <- if (is.null(input$cluster_select)) "All" else input$cluster_select
          plot_cec_grid_qual_barplot(
            grid_obj,
            C = sel$C,
            lambda = sel$lambda,
            shape_var = shape_var,
            cluster_selection = cluster_select
          )
        })

        output$lambda_plot <- shiny::renderPlot({
          sel <- selected()
          plot_cec_grid_selected_lambda_path_1d(
            grid_obj,
            C = sel$C,
            selected_lambda = sel$lambda
          )
        })

        output$selected_table <- shiny::renderTable({
          sel <- selected()
          s <- current_summary()
          idx <- cec_grid_row_nearest(s, sel$C, sel$lambda)
          row <- s[idx, , drop = FALSE]
          row[, c(
            "C", "C_source", "C_reused", "lambda", "REO", "k_hat",
            "stab_algo", "rsi", "sat", "n_sat_clusters", "keep",
            "fail_reason"
          ), drop = FALSE]
        }, digits = 4)
      }

      shiny::runApp(
        shiny::shinyApp(ui = ui, server = server),
        launch.browser = launch.browser
      )
    }

  }

}


