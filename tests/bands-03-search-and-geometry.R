library(CEClust)
N <- asNamespace("CEClust")

# A thin winning lambda interval must not disappear between grid values.
region <- N$bands_hull(c(0, .99999, 2), c(2, 1, 0), c(.9, 1.1))
stopifnot(nrow(region) == 3L,
          abs(region$lower[2] - .99999) < 1e-12,
          abs(region$upper[2] - 1.00001) < 1e-12)

# Direct finite-bank minima and the hull must agree away from crossings.
A <- c(-.2, .1, .9, 2)
H <- c(2, 1.5, .5, 0)
hull <- N$bands_hull(A, H, c(.01, 3))
for (i in seq_len(nrow(hull))) {
  lambda <- (hull$lower[i] + hull$upper[i]) / 2
  stopifnot(hull$partition[i] == which.min(A + lambda * H))
}

X <- matrix(c(-2.6, -2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2,
               1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6), ncol = 1)
store <- N$bands_store(tempfile("bands-search-test-"))
control <- bands_control(orders = 1:2, starts = 1, extra_starts = 1,
                         max_order = 4, max_iter = 20, max_anchors = 2,
                         C_points = 7, max_C_points = 15, refine_rounds = 0,
                         projection_rounds = 0, max_projection_queries = 10)
data_key <- digest::digest(X)
candidates <- N$bands_baselines(X, store, data_key)
two <- N$bands_candidate(X, rep(1:2, each = 8), store, "fixture", data_key)
candidates <- N$bands_add_candidate(candidates, two)
geometry <- N$bands_geometry(candidates, control, "none")
left <- geometry$tiles[geometry$tiles$column == 1L, ]
stopifnot(all(left$saturation == 1))
right <- N$bands_lines(candidates, max(geometry$C), "none")
upper <- geometry$window$lambda[2]
winner <- which.min(right$A + upper * right$H)
stopifnot(candidates[[winner]]$K == 1L)
for (i in seq_len(nrow(geometry$tiles))) {
  point <- geometry$tiles[i, ]
  scores <- vapply(candidates, function(candidate) {
    bands_score(X, N$bands_load_labels(candidate, store), point$C,
                 point$lambda)$J
  }, 0)
  stopifnot(abs(point$J - min(scores)) < 1e-10)
}

# Both branches start from the same extra allocation, not different random
# Refinement remains finite for valid density scales whose product overflows
# or underflows; this changes neither the data nor the requested search window.
extreme_control <- control
extreme_control$refine_rounds <- 1L
extreme_control$max_C_points <- 41L
X2 <- cbind(X[, 1], sin(seq_len(nrow(X))))
for (scale in c(1e-90, 1e90)) {
  scaled <- X2 * scale
  pool <- N$bands_baselines(scaled, store, digest::digest(scaled))
  candidate <- N$bands_candidate(scaled, rep(1:2, each = 8), store, "extreme_scale")
  pool <- N$bands_add_candidate(pool, candidate)
  refined <- N$bands_geometry(pool, extreme_control, "none")
  stopifnot(length(refined$history) > 0L,
            all(is.finite(refined$C) & refined$C > 0),
            all(is.finite(refined$tiles$J)))
}

# Both branches start from the same extra allocation, not different random
# draws hidden behind a common seed in a different position of the RNG stream.
pilot <- list(candidates = candidates, n = nrow(X))
anchors <- data.frame(C = 5, lambda = 1, key = two$key, K = 2L, id = 1L)
job1 <- list(name = "original_test", kind = "original", rep = 1L,
             seed = 123L, indices = NULL)
run1 <- N$bands_training_bank(X, job1, pilot, anchors, control, "none", store,
                              runtime = "contract-test")
warm <- run1$log[run1$log$stage == "warm_inner", ]
stopifnot(nrow(warm) == 2L, length(unique(warm$initial_key)) == 1L)
extra <- run1$log[run1$log$stage == "extra_target", ]
stopifnot(nrow(extra) == 2L,
          identical(sort(extra$branch), sort(c("free", "bounded"))),
          length(unique(extra$initial_key)) == 1L,
          length(unique(extra$seed)) == 1L,
          extra$C[extra$branch == "free"] == Inf,
          extra$C[extra$branch == "bounded"] == anchors$C)

# Structural symmetry: with identical training data and identical randomization,
# changing the original/bootstrap role does not change the estimator itself.
# The actual public procedure assigns different seeds to these roles.
job2 <- job1
job2$name <- "bootstrap_test"
job2$kind <- "bootstrap"
job2$indices <- seq_len(nrow(X))
run2 <- N$bands_training_bank(X, job2, pilot, anchors, control, "none", store,
                              runtime = "contract-test")
stopifnot(identical(run1$inner_indices, run2$inner_indices),
          identical(names(run1$candidates), names(run2$candidates)))
fields <- c("anchor", "stage", "branch", "shot", "C", "lambda", "seed",
            "initial_key", "candidate_key", "valid", "stop")
stopifnot(identical(run1$log[fields], run2$log[fields]))
stopifnot(length(run1$inner_indices) == nrow(X))

# Real jobs have independent randomization keys and retain supplied outer draws.
indices <- list(seq_len(nrow(X)), rev(seq_len(nrow(X))))
jobs <- N$bands_jobs(nrow(X), 2L, 2L, 55L, indices)
stopifnot(length(unique(vapply(jobs, `[[`, 0L, "seed"))) == 4L,
          identical(lapply(jobs[3:4], `[[`, "indices"), indices))

# A bank checkpoint can be reused exactly but cannot silently accept changed
# controls. The original archive is not overwritten by this failure.
again <- N$bands_training_bank(X, job1, pilot, anchors, control, "none", store,
                               runtime = "contract-test")
stopifnot(identical(run1, again))
changed <- control
changed$extra_starts <- changed$extra_starts + 1L
bad <- try(N$bands_training_bank(X, job1, pilot, anchors, changed, "none", store,
                                 runtime = "contract-test"), silent = TRUE)
stopifnot(inherits(bad, "try-error"))
