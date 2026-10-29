library(CEClust)
N <- asNamespace("CEClust")

near <- function(a, b, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(a, b, tolerance = tolerance,
                             check.attributes = FALSE)))
}

# Exhaustive agreements have known values for these small partitions.
allocations <- lapply(list(c(1, 1, 2, 2), c(2, 2, 1, 1), c(1, 1, 1, 2),
                         c(1, 1, 2, 2), c(1, 2, 2, 2)), N$bands_labels)
keys <- vapply(allocations, digest::digest, "", algo = "xxhash64")
agreements <- N$bands_agreements(allocations, keys, 3, 2, N$bands_cache(4096))
stopifnot(agreements$original_pairs == 3, agreements$bootstrap_pairs == 6)
near(agreements$stab_algo, 5 / 6)
near(agreements$stab_boot, 19 / 24)
near(agreements$RSI, .95)
near(agreements, N$bands_agreements(allocations, keys, 3, 2, N$bands_cache(0)))
forty <- rep(allocations[1], 40)
forty_keys <- rep(keys[1], 40)
twenty <- N$bands_agreements(forty, forty_keys, 20, 20, N$bands_cache(4096))
stopifnot(twenty$original_pairs == 190, twenty$bootstrap_pairs == 400,
          twenty$stab_algo == 1, twenty$RSI == 1)

# The new policy has no champion-agreement gate. The candidate chosen for the
# reference is deliberately different from all independent-run winners here.
store <- N$bands_store(tempfile("bands-gate-test-"))
X <- matrix(c(-2, -1, 1, 2), ncol = 1)
one <- N$bands_candidate(X, rep(1L, 4), store, "independent")
two <- N$bands_candidate(X, c(1L, 1L, 2L, 2L), store, "reference")
bank <- setNames(list(one), one$key)
runs <- lapply(seq_len(4), function(i) {
  list(name = paste0("run", i), kind = if (i <= 2) "original" else "bootstrap",
       signature = paste0("signature", i), candidates = bank)
})
control <- bands_control(orders = 1:2, starts = 1, extra_starts = 0,
                         max_order = 4, max_iter = 20, max_anchors = 2,
                         C_points = 3, max_C_points = 7, refine_rounds = 0,
                         projection_rounds = 0, max_projection_queries = 10)
point <- N$bands_evaluate_point(X, runs, setNames(list(two), two$key), 10, 1,
                                2, 2, control, "none", store, N$bands_cache(0))
stopifnot(point$valid, point$passed, !point$conditional,
          point$reference_key == two$key,
          all(point$winners$candidate_key == one$key),
          point$diagnostics$stab_algo == 1, point$diagnostics$RSI == 1)

# A projected bootstrap candidate is refitted on Z, not scored using the
# parameters or objective from its bootstrap training sample.
Z <- matrix(c(-4, -3, -2, 0, 1, 3, 4, 5), ncol = 1)
training <- Z[c(1, 1, 2, 3, 6, 7, 8, 8), , drop = FALSE]
original_one <- N$bands_candidate(Z, rep(1L, 8), store, "original-reference")
bootstrap_two <- N$bands_candidate(training, rep(1:2, each = 4), store,
                                    "bootstrap-reference")
original_bank <- setNames(list(original_one), original_one$key)
bootstrap_bank <- setNames(list(bootstrap_two), bootstrap_two$key)
reference_runs <- list(
  list(name = "original_1", kind = "original", signature = "o1", candidates = original_bank),
  list(name = "original_2", kind = "original", signature = "o2", candidates = original_bank),
  list(name = "bootstrap_1", kind = "bootstrap", signature = "b1", candidates = bootstrap_bank),
  list(name = "bootstrap_2", kind = "bootstrap", signature = "b2", candidates = bootstrap_bank))
reference_control <- control
reference_control$projection_rounds <- 1L
reference_control$max_projection_queries <- 2L
reference <- N$bands_reference(Z, reference_runs, data.frame(C = 10, lambda = 1),
                               reference_control, "none", store)
from_bootstrap <- reference$provenance[reference$provenance$origin == "bootstrap_projection", ]
stopifnot(nrow(from_bootstrap) == 2L,
          length(unique(from_bootstrap$key)) == 1L,
          all(from_bootstrap$source_key == bootstrap_two$key))
refitted <- reference$candidates[[from_bootstrap$key[1]]]
stopifnot(identical(refitted$data_key, digest::digest(Z)))
refitted_labels <- N$bands_load_labels(refitted, store)
near(N$bands_candidate_score(refitted, 10, 1, "none")$J,
     bands_score(Z, refitted_labels, 10, 1)$J)
stopifnot(length(reference_runs[[3]]$candidates) == 1L,
          reference_runs[[3]]$candidates[[1]]$key == bootstrap_two$key)

# A query budget must never stop partway through the B banks of one point.
for (budget in c(1L, 3L)) {
  batched_control <- reference_control
  batched_control$max_projection_queries <- budget
  batched <- N$bands_reference(Z, reference_runs, data.frame(C = 10, lambda = 1),
                               batched_control, "none", store)
  expected_queries <- if (budget == 1L) 0L else 2L
  stopifnot(batched$queries %% 2L == 0L, batched$queries <= budget,
            batched$queries == expected_queries,
            "projection_query_budget" %in% batched$limits)
}

# Small complete runs exercise RNG preservation, nested indices, checkpoints,
# compact output and the distinction between measured points and region display.
X <- matrix(c(seq(-3, -1, length.out = 12), seq(1, 3, length.out = 12)), ncol = 1)
indices <- list(rep(seq_len(12), each = 2), rep(13:24, each = 2))
set.seed(910)
before <- .Random.seed
directory <- tempfile("bands-pipeline-test-")
fit <- cec_bands(X, L = 2, B = 2, workers = 1, seed = 19,
                 pilot_indices = seq_len(nrow(X)), bootstrap_indices = indices,
                 control = control, directory = directory)
stopifnot(identical(before, .Random.seed), length(fit$runs) == 4,
          identical(lapply(fit$runs[3:4], `[[`, "indices"), indices),
          all(vapply(fit$runs, function(run) length(run$inner_indices), 0L) == nrow(X)))
stopifnot(all(is.na(fit$tiles$stab_algo[!fit$tiles$measured])),
          all(is.na(fit$tiles$RSI[!fit$tiles$measured])))
for (run in fit$runs) for (candidate in run$candidates)
  stopifnot(is.null(candidate$labels), !is.null(candidate$label_file))
has_environment <- function(value) {
  if (is.environment(value)) return(TRUE)
  if (is.list(value)) return(any(vapply(value, has_environment, TRUE)))
  FALSE
}
stopifnot(!has_environment(fit))
validation <- bands_validate(fit)
stopifnot(validation$passed)
partition_table <- bands_partitions(fit)
stopifnot(all(c("witness_saturation", "min_retained_saturation",
                "max_retained_saturation") %in% names(partition_table)),
          all(partition_table$max_retained_saturation <= control$saturation_max))
again <- cec_bands(X, L = 2, B = 2, workers = 1, seed = 19,
                   pilot_indices = seq_len(nrow(X)), bootstrap_indices = indices,
                   control = control, directory = directory)
stopifnot(identical(fit, again), identical(before, .Random.seed))
bad <- try(cec_bands(X + .1, L = 2, B = 2, workers = 1, seed = 19,
                     pilot_indices = seq_len(nrow(X)), bootstrap_indices = indices,
                     control = control, directory = directory), silent = TRUE)
stopifnot(inherits(bad, "try-error"), identical(before, .Random.seed))

# Every archived point agrees with a direct evaluation using cache disabled.
for (i in seq_len(nrow(fit$points))) {
  q <- fit$points[i, ]
  saved <- readRDS(file.path(fit$store$path, q$checkpoint))$value
  direct <- N$bands_evaluate_point(fit$X, fit$runs, fit$reference$candidates,
                                   q$C, q$lambda, fit$L, fit$B, fit$control,
                                   fit$penalty, fit$store, N$bands_cache(0))
  stopifnot(identical(saved$valid, direct$valid),
            identical(saved$passed, direct$passed))
  if (saved$valid) near(saved$diagnostics, direct$diagnostics)
}
if (nrow(fit$points)) {
  q <- fit$points[1, ]
  direct <- bands_at(fit, q$C, q$lambda)
  stopifnot(identical(direct$valid, q$valid))
}

all_control <- control
all_control$diag_mode <- "all"
all_fit <- cec_bands(X, L = 2, B = 2, workers = 1, seed = 19,
                     pilot_indices = seq_len(nrow(X)), bootstrap_indices = indices,
                     control = all_control, directory = tempfile("bands-all-test-"))
eligible <- all_fit$tiles$saturation <= all_control$measurement_saturation
stopifnot(all(all_fit$tiles$measured[eligible]),
          !any(all_fit$tiles$measured[!eligible]),
          all(all_fit$tiles$witness[all_fit$tiles$keep] == all_fit$tiles$id[all_fit$tiles$keep]))
stopifnot(identical(names(fit$reference$candidates), names(all_fit$reference$candidates)))

figure <- tempfile(fileext = ".pdf")
grDevices::pdf(figure)
bands_plot(fit, filtered = FALSE)
bands_plot(fit, filtered = TRUE)
grDevices::dev.off()
stopifnot(file.exists(figure), file.info(figure)$size > 0)
