library(CEClust)
N <- asNamespace("CEClust")
E <- if (exists("bands_cleanup", .GlobalEnv, inherits = FALSE)) .GlobalEnv else N
split_labels <- get("bands_cleanup_split", E)
reject <- get("bands_cleanup_reject", E)
cleanup <- get("bands_cleanup", E)
set.seed(731)
X <- matrix(c(rnorm(20, -3, .2), rnorm(20, 3, .2)), ncol = 1)
z <- rep(1:2, each = 20)
for (mode in c("median", "kmeans")) {
  s <- split_labels(X, z, 1:2, mode, 27)
  stopifnot(max(s) == 4L, all(vapply(split(z, s), function(x) length(unique(x)) == 1L, TRUE)),
    identical(s, split_labels(X, z, 1:2, mode, 27)))
  stopifnot(max(split_labels(matrix(1, 5, 2), rep(1L, 5), 1L, mode, 27)) == 1L)
}
# Require improvement AND (higher entropy OR saturation). Exercise all branches.
r <- list(J=1,entropy=1)
stopifnot(reject(r,list(J=.9,entropy=1.1,saturation=0),1e-8,0)["remove"],
 reject(r,list(J=.9,entropy=.9,saturation=.1),1e-8,0)["remove"],
 !reject(r,list(J=2,entropy=2,saturation=.1),1e-8,0)["remove"],
 !reject(r,list(J=.9,entropy=.9,saturation=0),1e-8,0)["remove"],
 !reject(r,list(J=1-1e-10,entropy=2,saturation=.1),1e-8,0)["remove"])
ctrl <- bands_control(orders = 1:2, starts = 1, extra_starts = 0, max_anchors = 2,
  C_points = 3, max_C_points = 5, refine_rounds = 0, projection_rounds = 0,
  max_projection_queries = 2, max_order = 4, max_iter = 30,
  algo_min = .5, rsi_min = .5)
fit <- cec_bands(X, L = 2, B = 2, workers = 1, seed = 71, control = ctrl,
  directory = tempfile("cleanup-source-"))
before <- digest::digest(fit, algo = "sha256")
archive <- tempfile("cleanup-test-"); set.seed(101); rng <- .Random.seed
out <- cleanup(fit, archive, max_rounds = 1, trials = 2, workers = 2, seed = 73)
stopifnot(identical(rng,.Random.seed))
paths <- list.files(file.path(archive, "jobs"), full.names = TRUE)
times <- file.info(paths)$mtime
again <- cleanup(fit, archive, max_rounds = 1, trials = 2, workers = 2, seed = 73)
stopifnot(identical(before, digest::digest(fit, algo = "sha256")),
  identical(out, again), identical(times, file.info(paths)$mtime),
  identical(out$diagnostics_recomputed, FALSE),
  nrow(out$remaining)==sum(!out$partitions$removed))
# Reuse checks both the source identity/settings and the exact split initialisation.
copy <- cleanup(fit,tempfile("cleanup-reuse-"),max_rounds=1,trials=2,workers=1,seed=73,
 reuse_directory=archive)
stopifnot(all(copy$descents$reused),identical(copy$partitions,out$partitions))
failed <- try(cleanup(fit,fit$store$path,max_rounds=1,trials=2),silent=TRUE)
stopifnot(inherits(failed,"try-error"))
empty <- fit;empty$tiles$keep[] <- FALSE
empty_result <- cleanup(empty,tempfile("cleanup-empty-"),max_rounds=1,trials=2)
stopifnot(nrow(empty_result$partitions)==0,nrow(empty_result$remaining)==0)
cat("Cleanup removal contracts passed.\n")
