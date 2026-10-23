reference_fixture <- function() {
  X <- CECwith_seed(18,matrix(c(rnorm(12,-3,.2),rnorm(12,3,.2)),ncol=1))
  labels <- rep(1:2,each=12L)
  control <- CECoptimization_control(list(Rmax=2,Ninit=1,max_iter=15),24,10,1,15)
  cfg <- list(X=X,boots=list(rep(c(1:10,13:14),2),24:1),seed=31L,
    settings=CECadaptive_cem_settings(control),penalty='none',signature='fixture',cache=NULL,resume=FALSE)
  z <- list(original_labels=rep(list(rep(1L,24)),3),
    diagnostics=CECstability(rep(list(rep(1L,24)),3),rep(list(rep(1L,24)),2)))
  list(cfg=cfg,z=z,reference=list(labels=labels),cell=data.frame(id=7L,C=10,lambda=1))
}

test_that('all three trigger conditions are checked before any final filtering', {
  f <- reference_fixture();calls <- 0L
  runner <- function(jobs,worker) {
    calls <<- calls+length(jobs)
    lapply(jobs,function(j)list(part=list(labels=f$reference$labels)))
  }
  a <- CECreference_stability(f$z,f$reference,f$cfg,f$cell,'match',.85,TRUE,runner)
  expect_true(a$reference_check$applied);expect_equal(calls,2L)
  expect_equal(a$reference_check$independent_stability,1)
  expect_true(a$reference_check$independent_stability_passed)
  expect_equal(a$reference_check$independent_rsi,1)
  expect_true(a$reference_check$independent_rsi_passed)
  expect_equal(a$diagnostics$Stab_algo,1);expect_equal(a$diagnostics$RSI,1)
  expect_identical(a$diagnostics_raw,f$z$diagnostics)
  expect_equal(a$diagnostics$original_pairs,0L);expect_equal(a$diagnostics$cross_pairs,2L)
  expect_identical(a$original_labels,f$z$original_labels)
  again <- CECreference_stability(a,f$reference,f$cfg,f$cell,'match',.85,TRUE,runner)
  expect_identical(again,a);expect_equal(calls,2L)
  # Equality with the similarity threshold does not trigger.
  b <- CECreference_stability(f$z,f$reference,f$cfg,f$cell,'match',.5,TRUE,runner)
  expect_false(b$reference_check$applied);expect_equal(calls,2L)
  # The one-group candidate is better when lambda is sufficiently large.
  cell <- f$cell;cell$lambda <- 10
  b <- CECreference_stability(f$z,f$reference,f$cfg,cell,'match',.85,TRUE,runner)
  expect_false(b$reference_check$original_not_better);expect_false(b$reference_check$applied)
  expect_equal(calls,2L)
  b <- CECreference_stability(f$z,f$reference,f$cfg,f$cell,'match',.85,FALSE,runner)
  expect_false(b$reference_check$applied);expect_identical(b$diagnostics,f$z$diagnostics)
  cfg <- f$cfg;cfg$boots <- list()
  expect_false(CECreference_stability(f$z,f$reference,cfg,f$cell,'match',.85,TRUE,runner)$reference_check$applied)
  # Lack of repetitions never establishes algorithmic stability.
  z <- f$z;z$diagnostics$Stab_algo <- NA_real_
  expect_false(CECreference_stability(z,f$reference,f$cfg,f$cell,'match',.85,TRUE,runner)$reference_check$applied)
})

test_that('independent algorithmic stability must strictly exceed its threshold', {
  f <- reference_fixture()
  runner <- function(...)stop('An unstable independent campaign must not launch reference bootstraps')
  for(value in c(.714342105263158,.471236842105263,.8,.85,NA_real_,NaN,Inf)) {
    z <- f$z;z$diagnostics$Stab_algo <- value
    a <- CECreference_stability(z,f$reference,f$cfg,f$cell,'match',.85,TRUE,runner)
    expect_false(a$reference_check$applied)
    expect_false(a$reference_check$independent_stability_passed)
    expect_identical(a$diagnostics,z$diagnostics)
    expect_null(a$reference_bootstrap)
  }
  z <- f$z;z$diagnostics$Stab_algo <- .8500001
  a <- CECreference_stability(z,f$reference,f$cfg,f$cell,'match',.85,TRUE,
    function(jobs,worker)lapply(jobs,function(j)list(part=list(labels=f$reference$labels))))
  expect_true(a$reference_check$applied)
  expect_true(a$reference_check$independent_stability_passed)
  # A previous fixed value of one cannot replace the raw stability in the gate.
  a$diagnostics_raw$Stab_algo <- .8
  b <- CECreference_stability(a,f$reference,f$cfg,f$cell,'match',.85,TRUE,runner)
  expect_false(b$reference_check$applied)
  expect_equal(b$reference_check$independent_stability,.8)
  expect_identical(b$diagnostics,a$diagnostics_raw)
  expect_null(b$reference_bootstrap)
})

test_that('one reference descent uses indexed labels and frozen projection without multistart', {
  f <- reference_fixture();cache <- tempfile('reference_cache');dir.create(cache)
  f$cfg$cache <- cache;f$cfg$resume <- TRUE
  local_mocked_bindings(q_fit=function(...)stop('Multistart must not be called'),
    ca_subdivide=function(...)stop('Subdivision must not be called'),.package='CEClust')
  for(penalty in c('none','bic')) {
    f$cfg$penalty <- penalty
    z <- CECreference_stability(f$z,f$reference,f$cfg,f$cell,'match',.85)
    expect_true(z$reference_check$applied)
    expect_length(z$reference_bootstrap,2L)
    for(b in 1:2) {
      r <- z$reference_bootstrap[[b]]
      expect_identical(r$initial_labels,ca_labels(f$reference$labels[f$cfg$boots[[b]]]))
      expect_identical(r$bootstrap_indices,f$cfg$boots[[b]])
      expect_equal(r$descent_count,1L)
      pred <- predict_frozen(f$cfg$X,r$fit$params,'CEC',r$projection_seed)
      expect_identical(r$part$labels,ca_labels(pred$labels))
      expect_lte(r$fit$part$J,r$fit$initial_J+1e-10)
    }
    expect_identical(z$diagnostics$RSI,mean(vapply(z$reference_bootstrap,
      function(r)ca_sim(f$reference$labels,r$part$labels),0)))
    resumed <- CECreference_stability(f$z,f$reference,f$cfg,f$cell,'match',.85)
    expect_identical(z,resumed)
  }
  expect_length(list.files(cache,pattern='[.]rds$'),4L)
})

test_that('a passing ordinary RSI may be recalculated before regional traversal advances', {
  f <- reference_fixture();s <- list(mesh=ceca_mesh(exp(0:2),c(1,2)))
  s$parts <- setNames(rep(list(list(labels=f$reference$labels,K=2L,saturation=0)),6),1:6)
  # This example's original labels have length 24, independent of mesh labels.
  s$parts <- lapply(s$parts,function(p)list(labels=f$reference$labels,K=2L,saturation=0))
  calls <- integer()
  measure <- function(i) {
    calls <<- c(calls,i)
    z <- f$z;z$diagnostics$Stab_algo <- .9;z$diagnostics$RSI <- .97
    CECreference_stability(z,f$reference,f$cfg,f$cell,'match',.85,TRUE,
      function(jobs,worker)lapply(jobs,function(j)list(part=list(labels=f$reference$labels))))
  }
  a <- CECadaptive_component_diagnostics(s,rep(TRUE,6),measure,.85,.95)
  expect_length(calls,1L);expect_equal(a$plan$table$status,'accepted')
  expect_true(all(a$values$RSI==1))
})

test_that('regional traversal continues after an independently unstable point', {
  f <- reference_fixture();s <- list(mesh=ceca_mesh(exp(0:2),c(1,2)))
  s$parts <- setNames(rep(list(list(labels=f$reference$labels,K=2L,saturation=0)),6),1:6)
  calls <- integer()
  measure <- function(i) {
    calls <<- c(calls,i)
    z <- f$z;z$diagnostics$Stab_algo <- .8;z$diagnostics$RSI <- 1
    CECreference_stability(z,f$reference,f$cfg,f$cell,'match',.85,TRUE,
      function(...)stop('Reference bootstrap must not run'))
  }
  a <- CECadaptive_component_diagnostics(s,rep(TRUE,6),measure,.85,.95)
  expect_length(calls,6L)
  expect_false(any(a$plan$table$status=='accepted'))
})

test_that('independent RSI must strictly exceed the configured threshold', {
  f <- reference_fixture()
  runner <- function(...)stop('A failing raw RSI must not launch conditional bootstraps')
  for(value in c(.0811,.2,.949999,.95,NA_real_,NaN,Inf)) {
    z <- f$z;z$diagnostics$RSI <- value
    a <- CECreference_stability(z,f$reference,f$cfg,f$cell,'match',.85,TRUE,runner)
    expect_false(a$reference_check$applied)
    expect_false(a$reference_check$independent_rsi_passed)
    expect_identical(a$diagnostics,z$diagnostics)
    expect_null(a$reference_bootstrap)
  }
  for(pair in list(c(.950001,.95),c(.81,.8),c(1.02,1.01))) {
    z <- f$z;z$diagnostics$RSI <- pair[1]
    a <- CECreference_stability(z,f$reference,f$cfg,f$cell,'match',.85,TRUE,
      function(jobs,worker)lapply(jobs,function(j)list(part=list(labels=f$reference$labels))),rsi_min=pair[2])
    expect_true(a$reference_check$applied)
    expect_equal(a$reference_check$independent_rsi_threshold,pair[2])
  }
  a$diagnostics_raw$RSI <- .0811
  b <- CECreference_stability(a,f$reference,f$cfg,f$cell,'match',.85,TRUE,runner)
  expect_false(b$reference_check$applied)
  expect_equal(b$reference_check$independent_rsi,.0811)
  expect_identical(b$diagnostics,a$diagnostics_raw)
  expect_null(b$reference_bootstrap)
})

test_that('regional traversal does not rescue a low independent RSI', {
  f <- reference_fixture();s <- list(mesh=ceca_mesh(exp(0:2),c(1,2)))
  s$parts <- setNames(rep(list(list(labels=f$reference$labels,K=2L,saturation=0)),6),1:6)
  calls <- integer()
  measure <- function(i) {
    calls <<- c(calls,i)
    z <- f$z;z$diagnostics$RSI <- .0811
    CECreference_stability(z,f$reference,f$cfg,f$cell,'match',.85,TRUE,
      function(...)stop('A failed RSI must not be corrected'),rsi_min=.95)
  }
  a <- CECadaptive_component_diagnostics(s,rep(TRUE,6),measure,.85,.95)
  expect_length(calls,6L)
  expect_false(any(a$plan$table$status=='accepted'))
})

test_that('fixed grids forward the RSI threshold to the conditional rule', {
  f <- reference_fixture();X <- f$cfg$X
  part <- function(labels)list(labels=labels,local_labels=labels,
    stats=ca_stats(X,labels),local_stats=ca_stats(X,labels))
  banks <- c(rep(list(list(candidates=list(part(rep(1L,24))),cell_winners=c(1L,1L))),3),
    rep(list(list(candidates=list(part(f$reference$labels)),cell_winners=c(1L,1L))),2))
  for(method in c('cem','em_bank'))for(penalty in c('none','bic')) {
    for(cutoff in c(.49,.5,.95,1.1)) {
      a <- CECfit_solver_grid(X,10,c(1,10),k0=3,B=2,optimization=method,penalty=penalty,
        banks=banks,bootstrap_indices=list(1:24,24:1),rsi_threshold=cutoff,
        optimizer_control=list(Rmax=2,Ninit=1,max_iter=15))
      expect_identical(a$summary$reference_bootstrap_applied,c(cutoff<.5,FALSE))
      expect_equal(a$summary$rsi_raw[1],.5)
      expect_equal(a$measurements[['1']]$reference_check$independent_rsi_threshold,cutoff)
      if(cutoff<.5) {
        expect_gt(a$summary$rsi[1],.95)
        expect_length(a$measurements[['1']]$reference_bootstrap,2L)
      } else {
        expect_equal(a$summary$rsi[1],.5)
        expect_null(a$measurements[['1']]$reference_bootstrap)
      }
      if(cutoff>1)expect_false(any(a$summary$keep))
    }
    b <- CECfit_solver_grid(X,10,c(1,10),k0=3,B=2,optimization=method,penalty=penalty,
      banks=banks,bootstrap_indices=list(1:24,24:1),reference_bootstrap=FALSE,
      optimizer_control=list(Rmax=2,Ninit=1,max_iter=15))
    expect_identical(a$winners,b$winners);expect_identical(a$selected,b$selected)
    expect_equal(b$summary$rsi[1],.5)
  }
})

test_that('adaptive CEM, EM and rescoring retain the configured RSI gate', {
  actual <- CECreference_stability;seen <- numeric()
  local_mocked_bindings(CECreference_stability=function(...,rsi_min=.95) {
    seen <<- c(seen,rsi_min)
    actual(...,rsi_min=rsi_min)
  },.package='CEClust')
  X <- matrix(seq(-1,1,length.out=12),ncol=1)
  for(method in c('cem','em_bank'))for(strategy in c('components','cells')) {
    seen <- numeric()
    a <- CECfitBoundGrid(Z=X,C_grid=c(.2,3),lambda_grid=c(.8,1.2),
      optimization=method,grid_strategy='adaptive',k0=2,B=2,algo_seed=63,
      rsi_policy='all',rsi_threshold=.73,r0=1,Nshots_fresh=1,Nloop=3,
      optimizer_control=list(Rmax=1,Ninit=1,max_iter=3),
      adaptive_control=list(max_cells=4,max_rounds=2,max_extensions=0,diagnostic_strategy=strategy))
    expect_gt(length(seen),0L)
    expect_true(all(seen==.73))
    if(method=='em_bank') {
      seen <- numeric()
      b <- CECrescoreGrid(a,penalty='bic')
      expect_gt(length(seen),0L)
      expect_true(all(seen==.73))
      expect_true(b$meta$reused_banks)
    }
  }
})
