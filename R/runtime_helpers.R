# Runtime helpers ---------------------------------------------------------

CECdefault_n_cores <- function(max_cores = 8L, detected_cores = NULL) {
	max_cores <- suppressWarnings(as.integer(max_cores)[1])
	if (is.na(max_cores) || max_cores < 1L) {
		max_cores <- 1L
	}

	if (length(detected_cores) == 0L || is.null(detected_cores)) {
		detected_cores <- parallel::detectCores(logical = FALSE)
	}

	detected_cores <- suppressWarnings(as.integer(detected_cores)[1])
	if (is.na(detected_cores) || detected_cores < 1L) {
		return(1L)
	}

	max(1L, min(max_cores, detected_cores))
}

#' Configure the runtime mode used by CEClust.
#'
#' This helper provides an easy entry point for users who want either a
#' conservative pure-R setup or an automatic "use the fast path when available"
#' setup. The function immediately applies the backend choice through
#' [CECset_fast_backend()] and returns the number of cores that should be passed
#' to [CECfitLambdaGrid()].
#'
#' The `"base"` profile disables the optional compiled backend and forces
#' single-core execution. The `"auto"` profile activates the compiled backend
#' when available and chooses a conservative number of worker processes. The
#' `"fast"` profile requests the same backend explicitly and warns when the
#' compiled code is unavailable on the current machine.
#'
#' @param profile One of `"auto"`, `"base"`, or `"fast"`.
#' @param n_cores Optional manual number of cores. If `NULL`, a conservative
#'   automatic choice is made.
#' @param max_cores Upper bound used when `n_cores` is not supplied.
#' @param prefer_multicore Logical. If `FALSE`, the returned configuration uses
#'   one core even when more are available.
#' @param rebuild_cpp Logical. Kept for compatibility. If `TRUE`, the function
#'   refreshes the detection of the compiled backend that was built when the
#'   package was installed.
#'
#' @return A list with class `"CEC_runtime_config"` containing:
#' - `profile`: the selected execution profile;
#' - `use_cpp`: whether the optional compiled backend is active;
#' - `n_cores`: the recommended number of workers for
#'   [CECfitLambdaGrid()];
#' - `max_cores` and `prefer_multicore`: the decision inputs used by the helper;
#' - `performance`: the output of [CECperformanceInfo()].
#'
#' @seealso [CECperformanceInfo()], [CECset_fast_backend()],
#'   [CECreset_fast_backend()]
#' @export
CECconfigure_runtime <- function(profile = c("auto", "base", "fast"),
								 n_cores = NULL,
								 max_cores = 8L,
								 prefer_multicore = TRUE,
								 rebuild_cpp = FALSE) {
	profile <- match.arg(profile)

	if (identical(profile, "base")) {
		CECset_fast_backend(enabled = FALSE, rebuild = FALSE)
		info <- CECperformanceInfo()
		chosen_cores <- 1L
		use_cpp <- FALSE
	} else {
		CECset_fast_backend(enabled = TRUE, rebuild = rebuild_cpp)
		info <- CECperformanceInfo()
		use_cpp <- isTRUE(info$fast_backend_available)

		if (!is.null(n_cores)) {
			chosen_cores <- max(1L, suppressWarnings(as.integer(n_cores)[1]))
		} else if (!isTRUE(prefer_multicore)) {
			chosen_cores <- 1L
		} else {
			chosen_cores <- CECdefault_n_cores(
				max_cores = max_cores,
				detected_cores = info$detected_cores
			)
		}

		if (identical(profile, "fast") && !use_cpp) {
			warning(
				"The optional C++ backend is not available on this machine. ",
				"Falling back to pure R while keeping the requested core setting.",
				call. = FALSE
			)
		}
	}

	out <- list(
		profile = profile,
		use_cpp = use_cpp,
		n_cores = chosen_cores,
		max_cores = max_cores,
		prefer_multicore = isTRUE(prefer_multicore),
		performance = info
	)
	class(out) <- "CEC_runtime_config"
	out
}

#' @export
print.CEC_runtime_config <- function(x, ...) {
	cat("CEClust runtime configuration\n")
	cat("  profile :", x$profile, "\n")
	cat("  use_cpp :", x$use_cpp, "\n")
	cat("  n_cores :", x$n_cores, "\n")
	cat("  detected cores :", x$performance$detected_cores, "\n")
	invisible(x)
}

#' Simulate mixed multidimensional benchmark data.
#'
#' The returned data frame combines Gaussian numeric variables and categorical
#' variables sampled from cluster-specific multinomial distributions.
#'
#' @param n Number of rows to generate.
#' @param p_num Number of numeric variables.
#' @param p_fac Number of categorical variables.
#' @param K Number of latent groups used for simulation.
#' @param n_levels Number of levels for each categorical variable.
#' @param seed Random seed used for reproducibility.
#'
#' @return A data frame mixing numeric and factor columns. The returned data is
#'   designed for lightweight examples, tests, and documentation.
#'
#' @seealso [CECclassif()], [CECfitLambdaGrid()]
#' @export
simulate_multidim_benchmark_data <- function(n = 2500,
											 p_num = 8,
											 p_fac = 4,
											 K = 4,
											 n_levels = 4,
											 seed = 20260407) {
	set.seed(seed)

	z <- sample.int(K, size = n, replace = TRUE)
	mu <- matrix(rnorm(K * p_num, sd = 2), nrow = K, ncol = p_num)

	X_num <- matrix(0, nrow = n, ncol = p_num)
	for (k in seq_len(K)) {
		idx <- which(z == k)
		if (length(idx) == 0) {
			next
		}
		X_block <- matrix(rnorm(length(idx) * p_num), nrow = length(idx), ncol = p_num)
		X_block <- sweep(X_block, 2, mu[k, ], "+")
		X_num[idx, ] <- X_block
	}
	colnames(X_num) <- paste0("x", seq_len(p_num))

	level_names <- paste0("L", seq_len(n_levels))
	X_fac <- vector("list", p_fac)
	for (j in seq_len(p_fac)) {
		probs_j <- matrix(rexp(K * n_levels), nrow = K, ncol = n_levels)
		probs_j <- probs_j / rowSums(probs_j)

		values_j <- character(n)
		for (k in seq_len(K)) {
			idx <- which(z == k)
			if (length(idx) == 0) {
				next
			}
			values_j[idx] <- sample(
				level_names,
				size = length(idx),
				replace = TRUE,
				prob = probs_j[k, ]
			)
		}
		X_fac[[j]] <- factor(values_j, levels = level_names)
	}

	X_fac <- as.data.frame(X_fac, stringsAsFactors = TRUE)
	names(X_fac) <- paste0("f", seq_len(p_fac))

	data.frame(
		as.data.frame(X_num),
		X_fac,
		check.names = FALSE
	)
}

#' Sample a small diamonds data set ready for mixed-data clustering.
#'
#' This helper loads `ggplot2::diamonds`, randomly samples at most `n_max` rows,
#' and converts the ordered quality columns to plain factors.
#'
#' @param n_max Maximum number of rows to return.
#' @param seed Optional random seed.
#'
#' @return A data frame with four numeric variables and three factor variables.
#'
#' @details This helper requires the optional package \pkg{ggplot2}.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' CECconfigure_runtime("base")
#' Z_diamonds <- CECsample_diamonds(n_max = 80, seed = 1)
#' str(Z_diamonds)
#'
#' diag_diamonds <- CECfitLambdaGrid(
#'   Z = Z_diamonds,
#'   lambda_grid = seq(0.2, 0.8, by = 0.2),
#'   k0 = 2,
#'   B = 2,
#'   C = 10,
#'   r0 = 5,
#'   Nshots_fresh = 1,
#'   Nshots_warm = 1,
#'   Nloop = 10,
#'   familyType = "gaussAndDiscreteVector",
#'   seed = 1,
#'   silent = TRUE,
#'   verbose = FALSE,
#'   n_cores = 1L,
#'   checkpoint_dir = FALSE,
#'   auto_checkpoint = FALSE,
#'   show_progress = FALSE
#' )
#' head(diag_diamonds$summary)
#' @export
CECsample_diamonds <- function(n_max = 100, seed = NULL) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		stop("Package 'ggplot2' is required for CECsample_diamonds().")
	}

	diamonds <- ggplot2::diamonds

	if (!is.null(seed)) {
		set.seed(seed)
	}

	n_keep <- min(as.integer(n_max)[1], nrow(diamonds))
	idx <- sample.int(nrow(diamonds), size = n_keep)
	out <- diamonds[idx, c("carat", "depth", "table", "price", "cut", "color", "clarity")]
	out$cut <- factor(out$cut)
	out$color <- factor(out$color)
	out$clarity <- factor(out$clarity)
	out
}
