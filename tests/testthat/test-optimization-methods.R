solver_data <- function() {
  set.seed(680)
  cbind(c(rnorm(16,-2,.5),rnorm(16,2,.5)),rnorm(32,sd=.7))
}
solver_control <- function()list(seed=71,Rmax=3,Ninit=2,max_iter=30,max_batches=2)

test_that('all four optimizer and penalty choices return the stated hard criterion', {
  X <- solver_data();n <- nrow(X);d <- ncol(X);lambda <- .73
  for(method in c('cem','em_bank'))for(penalty in c('none','bic')) {
    f <- CECclassif(X,C=3,lambda=lambda,Nshots=2,Nloop=30,r0=3,
      optimization=method,penalty=penalty,optimizer_control=solver_control())
    expect_identical(f$optimization,method);expect_identical(f$penalty,penalty)
    expect_length(f$clusters,n);expect_equal(f$REO,length(unique(f$clusters)))
    weights <- as.numeric(table(f$clusters))/n
    conditional <- sum(vapply(seq_len(f$REO),function(k) {
      idx <- which(f$clusters==k)
      if(!is.null(attr(f$params$Sigma[[k]],'CEC_singular_support')))return(-length(idx)*log(3)/n)
      expect_gte(determinant(f$params$Sigma[[k]],logarithm=TRUE)$modulus[1],-d*log(2*pi)-2*log(3)-1e-8)
      -sum(mvtnorm::dmvnorm(X[idx,,drop=FALSE],f$params$m[k,],f$params$Sigma[[k]],log=TRUE))/n
    },0))
    entropy <- -lambda*sum(weights*log(weights))+conditional
    dimension_difference <- (f$REO-1L)*(d+1)*(d+2)/2
    expected_penalty <- if(penalty=='bic')dimension_difference*log(n)/(2*n) else 0
    expect_equal(f$entropy,entropy,tolerance=1e-8)
    expect_equal(f$penalty_value,expected_penalty,tolerance=1e-12)
    expect_equal(f$objective,entropy+expected_penalty,tolerance=1e-8)
    expect_equal(lambda*f$Hphi,f$objective,tolerance=1e-12)
    if(method=='cem')for(trace in f$solver$descent_traces) {
      expect_true(all(diff(trace$J)<0))
      expect_equal(trace$J,trace$CEC+trace$penalty_relative,tolerance=1e-12)
    }
    if(method=='em_bank')for(s in f$bank$starts)if(s$valid)expect_true(all(diff(s$trace)>=0))
  }
})

test_that('EM candidate generation is independent of penalty and hard scoring respects smaller ceilings', {
  X <- solver_data();control <- solver_control();control$Cmax <- 5
  a <- CECclassif(X,C=.4,lambda=.8,optimization='em_bank',penalty='none',optimizer_control=control)
  b <- CECclassif(X,C=.4,lambda=.8,optimization='em_bank',penalty='bic',optimizer_control=control)
  expect_identical(a$bank,b$bank)
  expect_lte(b$REO,a$REO)
  for(p in a$bank$candidates)if(!is.null(p$params))for(k in seq_along(p$params$weights))
    expect_gte(as.numeric(determinant(p$params$cov[,,k],logarithm=TRUE)$modulus),
      -ncol(X)*log(2*pi)-2*log(control$Cmax)-1e-8)
  expect_error(CECclassif(X,C=6,optimization='em_bank',optimizer_control=control),'Cmax')
})

test_that('explicit seeds leave caller RNG unchanged and reproduce the full bank', {
  X <- solver_data();set.seed(933);old <- .Random.seed
  a <- CECclassif(X,optimization='em_bank',optimizer_control=solver_control())
  expect_identical(.Random.seed,old)
  b <- CECclassif(X,optimization='em_bank',optimizer_control=solver_control())
  expect_identical(a$clusters,b$clusters);expect_identical(a$bank,b$bank)
})

test_that('singular data and small samples remain explicit rather than crashing native code', {
  for(X in list(matrix(1,6,2),cbind(1:6,2*(1:6)),matrix(c(0,0),1,2))) {
    z <- CECclassif(X,C=2,optimization='em_bank',optimizer_control=list(seed=1,Rmax=2,Ninit=1,max_iter=3))
    expect_true(is.finite(z$objective));expect_equal(z$saturation,1)
    expect_equal(length(z$clusters),nrow(X))
  }
  expect_error(emmap_fit_cpp(matrix(1,2,1),c(0L,1L),1),'labels')
  expect_error(emmap_fit_cpp(matrix(1,2,1),c(1L,1L),-1),'input')
  expect_error(CECclassif(data.frame(a=factor(c('a','b'))),optimization='em_bank'),'numeric')
  expect_error(CECclassif(matrix(NA_real_,2,2),optimization='em_bank'),'finite')
  expect_error(CECclassif(1:10,optimization='em_bank',optimizer_control=list(typo=3)),'Unknown')
  expect_error(CECclassif(1:10,optimization='em_bank',optimizer_control=list(3)),'named')
})

test_that('exhaustive stability uses original-bootstrap pairs and no self comparisons', {
  a <- list(c(1,1,2,2),c(1,2,1,2),c(1,2,2,1))
  b <- list(a[[1]],rep(1,4))
  z <- CECstability(a,b)
  expect_equal(z$Stab_algo,.5);expect_equal(z$Stab_boot,7/12);expect_equal(z$RSI,7/6)
  expect_equal(z$original_pairs,3);expect_equal(z$cross_pairs,6)
  expect_true(is.na(CECstability(a[1],b)$Stab_algo))
  expect_true(is.na(CECstability(a,list())$RSI))
  expect_error(CECstability(a,list(1:3)),'equal length')
})

test_that('grids preserve independent winners and recompute all diagnostics when penalty changes', {
  X <- solver_data()
  a <- CECfitBoundGrid(Z=X,C_grid=c(.2,3),lambda_grid=c(.7,1.8),k0=3,B=2,
    optimization='em_bank',penalty='none',optimizer_control=solver_control(),n_cores=1)
  b <- CECrescoreGrid(a,penalty='bic')
  expect_identical(a$banks,b$banks);expect_true(b$meta$reused_banks)
  expect_true(all(b$summary$REO<=a$summary$REO))
  for(g in list(a,b))for(i in seq_len(nrow(g$summary))) {
    labs <- lapply(g$winners[i,],function(id)g$parts[[id]]$labels)
    z <- CECstability(labs[1:3],labs[4:5])
    expect_equal(unname(unlist(g$summary[i,c('stab_algo','stab_boot','rsi')])),
      c(z$Stab_algo,z$Stab_boot,z$RSI),tolerance=1e-12)
    expect_equal(CECextractPartition(g,g$summary$C[i],g$summary$lambda[i]),g$parts[[g$selected[i]]]$labels)
  }
  expect_error(CECrescoreGrid(a,C_grid=30),'Cmax')
  damaged <- a;damaged$Z[1,1] <- 100
  expect_error(CECrescoreGrid(damaged),'changed')
  fixed <- CECextractBestPartitions(a$runs_by_C[[1]]$lambda_diag)
  expect_identical(CECrepairBestPartitionTrajectory(fixed,X)$best,fixed$best)
})

test_that('targeted RSI never fabricates stability for unmeasured or trivial cells', {
  X <- solver_data()
  g <- CECfitBoundGrid(Z=X,C_grid=c(.1,3),lambda_grid=c(.7,20),k0=2,B=2,
    optimization='em_bank',optimizer_control=solver_control(),rsi_policy='nontrivial_unsaturated')
  s <- g$summary
  skipped <- s$sat>0 | s$REO==1
  expect_true(any(skipped));expect_true(all(is.na(s$rsi[skipped])))
  expect_false(any(s$keep[skipped]));expect_false(any(s$measured[skipped]))
  for(C in g$C_grid)expect_identical(CECselectStableLambdas(g,C=C)$summary$keep,s$keep[s$C==C])
  shown <- cec_grid_apply_thresholds(g,rsi_threshold=0,stab_algo_threshold=0,sat_threshold=1)
  expect_false(any(shown$keep[skipped]))
  compact <- CECcompactGrid(g)
  expect_identical(compact$summary,g$summary)
  expect_identical(CECextractPartition(compact,3,.7),CECextractPartition(g,3,.7))
})

test_that('all explicit grid solver combinations propagate to classification and plotting', {
  X <- solver_data()
  for(method in c('cem','em_bank'))for(penalty in c('none','bic')) {
    g <- CECfitBoundGrid(Z=X,C_grid=c(.2,3),lambda_grid=c(.7,1.3),k0=2,B=1,
      optimization=method,penalty=penalty,optimizer_control=list(seed=4,Rmax=2,Ninit=1,max_iter=10,max_batches=1))
    expect_identical(g$optimization,method);expect_identical(g$penalty,penalty)
    for(run in g$runs_by_C)for(f in run$best_parts$best) {
      expect_identical(f$fit$params$optimization,method);expect_identical(f$fit$params$penalty,penalty)
    }
    pdf_file <- tempfile(fileext='.pdf');grDevices::pdf(pdf_file)
    expect_error(CECplotGrid(g),NA);expect_error(CECplotPartition(g,C=3,lambda=.7),NA)
    expect_error(CECplotPath(g,C=3,lambda=.7),NA);grDevices::dev.off()
    expect_gt(file.info(pdf_file)$size,1000)
  }
})

test_that('lambda entry point and presets also honor the requested optimizer', {
  X <- solver_data()[,1]
  d <- CECfitLambdaGrid(X,lambda_grid=c(.7,1.3),C=3,k0=2,B=1,
    optimization='em_bank',penalty='bic',optimizer_control=list(seed=3,Rmax=2,Ninit=1,max_iter=10),verbose=FALSE)
  expect_identical(d$optimization,'em_bank');expect_identical(d$penalty,'bic')
  expect_length(CECextractBestPartitions(d)$best,2)
  g <- CECfitPreset('iris',C_grid=1,lambda_grid=1,k0=2,B=0,n_cores=1,
    optimization='em_bank',optimizer_control=list(seed=2,Rmax=2,Ninit=1,max_iter=5),verbose=FALSE)
  expect_identical(g$optimization,'em_bank');expect_true(all(is.na(g$summary$rsi)))
})

test_that('archived CEM diagnostics now default to exhaustive comparisons', {
  set.seed(4);X <- matrix(rnorm(32),ncol=2)
  d <- CECfitLambdaGrid(X,C=1,lambda_grid=1,k0=3,B=2,r0=2,Nshots_fresh=1,
    Nshots_warm=1,Nloop=5,n_cores=1,checkpoint_dir=FALSE,verbose=FALSE)
  z <- CECstability(d$details[[1]]$partitions0,d$details[[1]]$partitionsB)
  expect_identical(d$stability_approach,'exhaustive')
  expect_equal(d$summary$stab0_mean,z$Stab_algo)
  expect_equal(d$summary$stabB_mean,z$Stab_boot)
})

test_that('legacy progress restores its context after errors without changing fits', {
  context <- CEClust:::.CEC_grid_progress
  saved <- context$state
  outer <- new.env(); outer$enabled <- FALSE
  inner <- new.env(); inner$enabled <- FALSE
  CEClust:::with_cec_grid_global_batch_progress({
    expect_identical(context$state,outer)
    expect_error(CEClust:::with_cec_grid_global_batch_progress(stop('probe'),inner),'probe')
    expect_identical(context$state,outer)
  },outer)
  expect_identical(context$state,saved)
  X <- solver_data()
  args <- list(Z=X,C_grid=.4,lambda_grid=1,k0=2,B=1,r0=2,familyType='gaussVector',
    Nshots_fresh=1,Nshots_warm=1,Nloop=3,n_cores=1,algo_seed=9,verbose=FALSE)
  silent <- do.call(CECfitBoundGrid,c(args,list(show_progress=FALSE)))
  output <- capture.output(shown <- do.call(CECfitBoundGrid,c(args,list(show_progress=TRUE,global_progress=TRUE))))
  expect_identical(shown$summary,silent$summary)
  expect_identical(context$state,saved)
  expect_true(any(grepl('elapsed',output)))
})
