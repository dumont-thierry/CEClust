test_that('soft EM preserves the frozen 9007 kernels and MAP allocations', {
  fixtures <- readRDS(test_path('fixtures','em_9007.rds'))
  for(z in fixtures) {
    observed <- do.call(emmap_fit_cpp,z$args)
    expect_identical(observed$labels,z$result$labels)
    expect_identical(observed$stop,z$result$stop)
    expect_identical(observed$converged,z$result$converged)
    expect_equal(observed$loglik_mean_trace,z$result$loglik_mean_trace,tolerance=1e-11)
    expect_equal(observed$means,z$result$means,tolerance=1e-11)
    expect_equal(observed$cov,z$result$cov,tolerance=1e-11)
    expect_equal(observed$weights,z$result$weights,tolerance=1e-11)
    expect_equal(observed$minimum_logdet_constraint_margin,
      z$result$minimum_logdet_constraint_margin,tolerance=1e-11)
  }
})

test_that('verified training MAP reuse preserves original and bootstrap roles', {
  set.seed(900820)
  X <- cbind(rnorm(60),rnorm(60));X[1:30,1] <- X[1:30,1]+2
  control <- CECoptimization_control(list(seed=900821,Rmax=4,Ninit=2,max_iter=20,Cmax=20),
    nrow(X),1,2,20,4)
  bank <- CECfit_em_bank(X,control)
  for(p in bank$candidates) {
    expect_identical(p$local_labels,p$labels)
    expect_identical(p$local_stats,p$stats)
    if(!is.null(p$params))expect_identical(p$labels,ca_labels(emmap_predict_cpp(X,
      p$params$means,p$params$cov,p$params$weights)))
  }
  index <- sample.int(nrow(X),replace=TRUE)
  boot <- CECfit_em_bank(X[index,,drop=FALSE],control,Xref=X)
  for(p in boot$candidates)if(!is.null(p$params)) {
    expect_identical(p$stats,ca_stats(X,p$labels))
    expect_identical(p$local_stats,ca_stats(X[index,,drop=FALSE],p$local_labels))
    expect_identical(p$labels,ca_labels(emmap_predict_cpp(X,p$params$means,
      p$params$cov,p$params$weights)))
  }
})

test_that('adaptive profile reuse keeps different bootstrap statistics separate', {
  X <- matrix(c(-3,-2.8,-2.2,-1.9,-1.2,-1,1,1.3,1.8,2.1,2.7,3),ncol=1)
  labels <- list(rep(1L,12),rep(1:2,each=6))
  candidates <- lapply(labels,function(z)list(labels=z,stats=ca_stats(X,z),local_stats=ca_stats(X,z)))
  other <- candidates
  other[[1]]$local_stats$logdet <- other[[1]]$local_stats$logdet+.8
  other[[2]]$local_stats$logdet <- other[[2]]$local_stats$logdet-.4
  banks <- list(list(candidates=candidates),list(candidates=other))
  pool <- CECgrid_collect(banks)
  for(penalty in c('none','bic')) {
    score <- CECadaptive_scorer(pool$parts,banks,pool$maps,X,list(Cmax=20),penalty)
    for(C in c(.02,.7,10))for(lambda in c(.3,1,2)) {
      observed <- score(C,lambda)
      truth <- CECscore_candidates(pool$parts,C,lambda,1,penalty)
      expect_identical(observed$id,truth$best)
      expect_equal(observed$J,truth$profile$J,tolerance=1e-12)
      for(k in seq_along(banks)) {
        local <- CECscore_candidates(banks[[k]]$candidates,C,lambda,1,penalty,'local_stats')
        expect_identical(observed$winners[k],pool$maps[[k]][local$best])
      }
    }
  }
})
