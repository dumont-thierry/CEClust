# Full-rank determinant-constrained fits and bounded affine-support scores.
# Singular empirical covariances are kept, with explicit support metadata.

CECdet_covariance_error <- function(message, rank = NA_integer_, dimension = NA_integer_) {
  stop(structure(
    list(message = message, call = NULL, rank = rank, dimension = dimension),
    class = c("CEC_determinant_fit_error", "error", "condition")
  ))
}

CECfit_covariance_determinant <- function(Sigma_emp, sigmaZmin) {
  S <- as.matrix(Sigma_emp)
  d <- nrow(S)
  if (d < 1L || ncol(S) != d || !is.numeric(S) || any(!is.finite(S))) {
    stop("Sigma_emp must be a finite, nonempty square numeric matrix.", call. = FALSE)
  }
  if (length(sigmaZmin) == 1L) sigmaZmin <- rep(sigmaZmin, d)
  if (length(sigmaZmin) != d || any(!is.finite(sigmaZmin)) || any(sigmaZmin < 0)) {
    stop("sigmaZmin must contain d finite, nonnegative scales.", call. = FALSE)
  }
  log_det_min <- 2 * sum(log(sigmaZmin))
  magnitude <- max(abs(S))
  if (magnitude > 0 && max(abs(S - t(S))) > 64 * d * .Machine$double.eps * magnitude) {
    stop("Sigma_emp must be symmetric.", call. = FALSE)
  }
  S <- (S + t(S)) / 2

  # Numerical rank detection does not inflate eigenvalues. Small eigenvalues
  # are treated as null for scoring; the stored empirical matrix is unchanged.
  eig <- eigen(S, symmetric = TRUE)
  spectrum <- eig$values
  rank_tolerance <- 8 * d * .Machine$double.eps * max(abs(spectrum))
  if (min(spectrum) < -rank_tolerance) {
    stop("Sigma_emp is not positive semidefinite.", call. = FALSE)
  }
  numerical_rank <- sum(spectrum > rank_tolerance)
  if (numerical_rank < d) {
    log_C <- -d / 2 * log(2 * pi) - log_det_min / 2
    if (!is.finite(log_C)) {
      stop("Singular affine-support scores require a finite positive C.", call. = FALSE)
    }
    attr(S, "CEC_singular_support") <- list(
      rank = numerical_rank,
      null_basis = eig$vectors[, spectrum <= rank_tolerance, drop = FALSE],
      rank_tolerance = rank_tolerance,
      support_tolerance = sqrt(rank_tolerance),
      log_C = log_C)
    return(list(Sigma_corrected = S, boundReached = TRUE,
                covariance_status = "singular_affine_support",
                log_scale = NA_real_, empirical_rank = numerical_rank,
                log_det_min = log_det_min, log_det_fitted = -Inf))
  }

  R <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(R)) {
    CECdet_covariance_error("Covariance factorization failed; no diagonal jitter was applied.",
                           rank = numerical_rank, dimension = d)
  }
  log_det_emp <- 2 * sum(log(diag(R)))
  log_scale <- max(0, (log_det_min - log_det_emp) / d)
  multiplier <- exp(log_scale)
  if (!is.finite(multiplier)) {
    stop("The required covariance inflation is outside floating-point range.", call. = FALSE)
  }
  out <- S * multiplier
  if (any(!is.finite(out))) {
    stop("The fitted covariance is outside floating-point range.", call. = FALSE)
  }
  R_fit <- tryCatch(chol(out), error = function(e) NULL)
  if (is.null(R_fit)) {
    CECdet_covariance_error("Fitted covariance factorization failed; no diagonal jitter was applied.",
                           rank = numerical_rank, dimension = d)
  }
  log_det_fit <- 2 * sum(log(diag(R_fit)))
  boundary_tolerance <- 1e-10 * max(1, abs(log_det_min))
  reached <- is.finite(log_det_min) &&
    (log_det_emp <= log_det_min + boundary_tolerance)
  list(Sigma_corrected = out, boundReached = reached,
       covariance_status = if (log_scale > 0) "inflated_exact" else "unchanged",
       log_scale = log_scale, empirical_rank = numerical_rank,
       log_det_min = log_det_min, log_det_fitted = log_det_fit)
}

weighted_cov <- function(res, w) {
  X <- as.matrix(res)
  if (!is.numeric(X) || nrow(X) == 0L || ncol(X) == 0L || any(!is.finite(X))) {
    stop("res must be a finite, nonempty numeric matrix.", call. = FALSE)
  }
  if (length(w) != nrow(X) || any(!is.finite(w)) || any(w < 0) ||
      !is.finite(sum(w)) || sum(w) <= 0) {
    stop("w must contain nonnegative finite weights with positive finite sum.", call. = FALSE)
  }
  keep <- w > 0
  X <- X[keep, , drop = FALSE]
  w <- w[keep] / sum(w[keep])
  # Center first around an observed row to avoid cancellation for translated
  # data, and to preserve exact zero covariance for identical observations.
  shifted <- sweep(X, 2, X[1L, ], "-")
  centered <- sweep(shifted, 2, colSums(shifted * w), "-")
  crossprod(centered, centered * w)
}

adjust_covariance <- function(res, sigmaZmin, w = NULL, sizeMinSigmaComp = 3) {
  X <- as.matrix(res)
  if (is.null(w)) w <- rep(1, nrow(X))
  # sizeMinSigmaComp is retained for call compatibility, not as a fit rule.
  out <- CECfit_covariance_determinant(weighted_cov(X, w), sigmaZmin)
  keep <- w > 0
  Xk <- X[keep, , drop = FALSE]
  wk <- w[keep] / sum(w[keep])
  m <- Xk[1L, ] + colSums(sweep(Xk, 2, Xk[1L, ], "-") * wk)
  out$Sigma_corrected <- CECset_singular_support_data(out$Sigma_corrected, Xk, m)
  out
}

CECadjust_covariance_from_matrix <- function(
  Sigma_emp, sigmaZmin, count, sizeMinSigmaComp = 3
) {
  # count and sizeMinSigmaComp no longer select an arbitrary small-group fit.
  CECfit_covariance_determinant(Sigma_emp, sigmaZmin)
}


# Only actual factorization failures reach this compatibility wrapper;
# singular covariances are now valid fits, not rejected shots.
CECrun_one_shot_with_regroup <- function(...) {
  tryCatch(CECrun_one_shot_with_regroup_strict(...),
    CEC_determinant_fit_error = function(e) {
      list(Hphi = Inf, phi = NULL,
           determinant_fit_issue = list(message = conditionMessage(e),
                                        rank = e$rank, dimension = e$dimension))
    })
}

CECattach_det_diagnostics <- function(out, issues, attempted) {
  issues <- Filter(Negate(is.null), issues)
  if (is.null(out)) {
    if (length(issues)) {
      warning(sprintf("No finite fit retained; %d of %d shots had covariance factorization failures.",
                      length(issues), attempted), call. = FALSE)
    }
    return(NULL)
  }
  out$determinant_diagnostics <- list(attempted_shots = attempted,
                                    rejected_shots = length(issues),
                                    issues = issues)
  out$singular_clusters <- which(vapply(out$params$Sigma,
    function(S) !is.null(attr(S, "CEC_singular_support", exact = TRUE)), logical(1)))
  out
}

CECset_singular_support_data <- function(S, X, mean) {
  support <- attr(S, "CEC_singular_support", exact = TRUE)
  if (is.null(support)) return(S)
  residual <- sweep(as.matrix(X), 2, mean, "-") %*% support$null_basis
  # Include the fitted observations despite numerical rank/centering error.
  support$support_tolerance <- max(support$support_tolerance,
                                   sqrt(rowSums(residual^2)))
  attr(S, "CEC_singular_support") <- support
  S
}

CECgaussian_reference_logdens <- function(X, mean, Sigma) {
  X <- as.matrix(X)
  support <- attr(Sigma, "CEC_singular_support", exact = TRUE)
  if (is.null(support)) {
    return(mvtnorm::dmvnorm(X, mean = mean, sigma = Sigma, log = TRUE))
  }
  centered <- sweep(X, 2, mean, "-")
  residual <- sqrt(rowSums((centered %*% support$null_basis)^2))
  scale <- pmax(apply(abs(X), 1, max), max(abs(mean)))
  tolerance <- support$support_tolerance + 64 * ncol(X) * .Machine$double.eps * scale
  on_support <- residual <= tolerance
  # Avoid roundoff in C reconstructed from delta at C=1.
  log_C <- if (abs(support$log_C) <= 32 * .Machine$double.eps) 0 else support$log_C
  ifelse(on_support, log_C, -Inf)
}

# The floor changes reassignment proposals only. Never use this matrix to
# evaluate the entropy: refit each proposed partition and use the raw scores.
CECreassignment_log_floor <- function(C, divisor = 1000) {
  if (length(C) != 1L || is.na(C) || C <= 0 ||
      length(divisor) != 1L || !is.finite(divisor) || divisor <= 1) {
    stop("Reassignment needs C > 0 and a finite floor divisor > 1.", call. = FALSE)
  }
  min(-log(C), log(C) - log(divisor))
}

CECproposal_logdens <- function(logdens, params) {
  if (!params$familyType %in% c("gaussVector", "gaussAndDiscreteVector")) {
    return(logdens)
  }
  if (params$familyType == "gaussAndDiscreteVector" && is.null(params$Sigma)) {
    return(logdens)
  }
  # Apply the same floor to every component, regardless of covariance rank.
  # For mixed data this is the complete conditional score, including factors.
  pmax(logdens, CECreassignment_log_floor(params$C))
}

CECcheck_reference_support <- function(logdens, nu) {
  active <- which(nu > 0)
  unsupported <- which(rowSums(is.finite(logdens[, active, drop = FALSE])) == 0L)
  if (length(unsupported)) {
    stop(structure(list(
      message = paste0("No positive component score for observation(s) ",
                       paste(head(unsupported, 10L), collapse = ", "),
                       "; posterior normalization is undefined."),
      call = NULL, observations = unsupported),
      class = c("CEC_no_supported_component", "error", "condition")))
  }
  invisible(NULL)
}

# Reassignment uses floored scores, but acceptance uses the original profiled
# determinant-constrained entropy, without a lower floor. A singular cluster
# contributes -log(C); a regular cluster uses its fitted Gaussian log-density.
CECclassifOneShot_affine_guarded <- function(
  Z, lambda, C, Cquali, r0, Nloop, phi0, familyType,
  displayPlotEntropy, backend_data
) {
  n <- CECget_n_obs(Z)
  use_fast <- isTRUE(backend_data$optimized)
  if (is.null(phi0)) {
    r <- if (is.null(r0)) floor(log(n) + 1) else r0
    phi0 <- phiInit(n, r)
  }
  evaluate <- function(phi, clusters = NULL) {
    if (is.null(clusters)) {
      clusters <- CECclusters_from_phi_if_hard(phi, n = n, r = length(phi) / n)
    }
    if (!is.null(clusters)) {
      clusters <- CECcompress_clusters(clusters)
      phi <- CECclusters_to_phi(clusters, length(unique(clusters)))
    }
    params <- if (use_fast && !is.null(clusters)) {
      CECoptParam_fast_from_clusters(backend_data, clusters, lambda, C, Cquali)
    } else {
      optParam(Z, phi, lambda, C, Cquali, familyType)
    }
    if (is.null(clusters) && !is.null(params$phi)) phi <- params$phi
    logdens <- if (use_fast) {
      CECcompute_logdens_matrix_fast(backend_data, params, lambda)
    } else {
      vectToMat(CECdens(Z, params, lambda, applyLog = TRUE)$density,
                length(params$nu))
    }
    phiM <- vectToMat(phi, length(params$nu))
    positive <- phiM > 0
    H <- -sum(params$nu[params$nu > 0] * log(params$nu[params$nu > 0])) -
      sum(phiM[positive] * logdens[positive]) / (n * lambda)
    params$phi <- phi
    list(phi = phi, clusters = clusters, params = params, logdens = logdens, Hphi = H)
  }
  current <- evaluate(phi0)
  # Random soft weights initialize parameters. Establish a hard partition
  # before recording the descent, rather than discarding a soft-start shot.
  if (is.null(current$clusters)) {
    proposed_logdens <- CECproposal_logdens(current$logdens, current$params)
    CECcheck_reference_support(proposed_logdens, current$params$nu)
    scores <- sweep(proposed_logdens / lambda, 2, log(current$params$nu), "+")
    current <- evaluate(NULL, CECcompress_clusters(max.col(scores, ties.method = "first")))
  }
  history <- current$Hphi
  rejected <- 0L
  stop_reason <- "iteration_limit"
  fingerprints <- if (is.null(current$clusters)) numeric() else
    CECcluster_fingerprint(current$clusters)
  for (iteration in seq_len(Nloop)) {
    proposed_logdens <- CECproposal_logdens(current$logdens, current$params)
    CECcheck_reference_support(proposed_logdens, current$params$nu)
    scores <- sweep(proposed_logdens / lambda, 2, log(current$params$nu), "+")
    next_clusters <- CECcompress_clusters(max.col(scores, ties.method = "first"))
    if (identical(next_clusters, current$clusters)) {
      stop_reason <- "unchanged_partition"
      break
    }
    candidate <- evaluate(NULL, next_clusters)
    if (!is.finite(candidate$Hphi) || candidate$Hphi > current$Hphi) {
      rejected <- rejected + 1L
      stop_reason <- "refitted_objective_increased"
      break
    }
    current <- candidate
    history <- c(history, current$Hphi)
    if (displayPlotEntropy) plot(history, main = "entropy", xlab = "step", ylab = "H")
    fp <- CECcluster_fingerprint(current$clusters)
    if (fp %in% fingerprints) {
      stop_reason <- "repeated_partition"
      break
    }
    fingerprints <- c(fingerprints, fp)
  }
  assigned_raw <- current$logdens[cbind(seq_len(n), current$clusters)]
  below_floor <- sum(assigned_raw < CECreassignment_log_floor(C))
  list(phi = current$phi, params = current$params, Hphi = current$Hphi,
       descent_diagnostics = list(H_trace = history,
         rejected_steps = rejected, stop_reason = stop_reason,
         floor_divisor = 1000, n_assigned_below_floor = below_floor,
         fraction_assigned_below_floor = below_floor / n,
         acceptance_criterion = "original_entropy_without_lower_floor"))
}
