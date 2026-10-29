agreement_fixture <- function() {
  mesh <- ceca_mesh(c(1,2),c(1,2))
  p <- list(labels=c(1L,1L,2L,2L),K=2L,saturation=0)
  list(mesh=mesh,parts=setNames(rep(list(p),4),mesh$leaves$id))
}

test_that('disagreement rejects a stable point and continues the component', {
  s <- agreement_fixture();calls <- integer()
  measure <- function(i) {
    calls <<- c(calls,i)
    list(diagnostics=list(Stab_algo=1,Stab_boot=.98,RSI=.98),
      reference_check=list(agreement=if(length(calls)==1).232 else .85))
  }
  z <- CECadaptive_component_diagnostics(s,rep(TRUE,4),measure,.85,.95,
    require_reference_agreement=TRUE)
  expect_length(calls,2L)
  expect_equal(z$trials$accepted,c(FALSE,TRUE))
  expect_equal(z$trials$rejection_reason,c('reference_disagreement','ok'))
  expect_equal(z$plan$table$attempts,2L)
  expect_true(all(z$reference_pass))
  expect_equal(z$measurements[[as.character(calls[1])]]$diagnostics$RSI,.98)
  expect_true(all(z$reference_agreement==.85))
  summary <- CECadaptive_component_summary(ceca_diagram(s),z)
  summary$stab_algo <- summary$rsi <- 1;summary$sat <- 0
  obj <- structure(list(summary=summary,meta=list(stab_algo_threshold=.85,rsi_threshold=.95,sat_threshold=0)),
    class=c('CEC_adaptive_grid','CEC_solver_grid'))
  expect_true(all(cec_grid_apply_thresholds(obj)$keep))
  expect_false(any(cec_grid_apply_thresholds(obj,stab_algo_threshold=.86)$keep))
})

test_that('exhausted disagreeing components remain rejected despite high raw diagnostics', {
  s <- agreement_fixture()
  for(a in c(.232,.84999,NA_real_,Inf)) {
    measure <- function(i)list(diagnostics=list(Stab_algo=1,Stab_boot=1,RSI=1),
      reference_check=list(agreement=a))
    z <- CECadaptive_component_diagnostics(s,rep(TRUE,4),measure,.85,.95,
      require_reference_agreement=TRUE)
    expect_equal(z$plan$table$status,'exhausted')
    expect_equal(z$plan$table$attempts,4L)
    expect_false(any(z$reference_pass))
    summary <- CECadaptive_component_summary(ceca_diagram(s),z)
    expect_false(any(summary$keep_harmonization))
    expect_true(all(summary$stab_algo==1 & summary$rsi==1))
  }
})

test_that('reference agreement does not replace either strict raw threshold', {
  s <- agreement_fixture()
  for(d in list(list(Stab_algo=.85,RSI=1),list(Stab_algo=1,RSI=.95))) {
    d$Stab_boot <- 1
    z <- CECadaptive_component_diagnostics(s,rep(TRUE,4),function(i)
      list(diagnostics=d,reference_check=list(agreement=1)),.85,.95,require_reference_agreement=TRUE)
    expect_equal(z$plan$table$status,'exhausted')
    expect_true(all(z$trials$rejection_reason=='stability_or_rsi'))
  }
  expect_true(CECadaptive_control(list(diagnostic_strategy='cells'),0)$require_reference_agreement)
  expect_error(CECfitBoundGrid(Z=matrix(1:12,ncol=1),grid_strategy='adaptive',
    adaptive_control=list(reference_bootstrap=TRUE,require_reference_agreement=TRUE)),
    'requires reference_bootstrap=FALSE')
})

test_that('an enabled reference gate rejects disagreement regardless of score ordering', {
  s <- agreement_fixture()
  expect_true(CECadaptive_control(list(),0)$require_reference_agreement)
  expect_false(CECadaptive_control(list(),0)$reference_bootstrap)
  for(scores in list(c(0,1),c(1,0),c(1,1))) {
    measure <- function(i)list(diagnostics=list(Stab_algo=1,Stab_boot=1,RSI=1),
      reference_check=list(agreement=.232,best_original_score=scores[1],reference_score=scores[2]))
    z <- CECadaptive_component_diagnostics(s,rep(TRUE,4),measure,.85,.95,
      require_reference_agreement=TRUE)
    expect_equal(z$plan$table$status,'exhausted')
    expect_false(any(z$trials$accepted))
  }
})

test_that('CEM and EM apply the gate before accepting a witness and retain it during rescoring', {
  X <- CECwith_seed(801,matrix(rnorm(20),ncol=1));calls <- 0L
  original <- CECreference_stability
  local_mocked_bindings(CECreference_stability=function(...) {
    z <- original(...);calls <<- calls+1L
    z$reference_check$agreement <- if(calls==1).232 else 1
    z
  },.package='CEClust')
  for(optimizer in c('cem','em_bank'))for(penalty in c('none','bic')) {
    calls <- 0L
    z <- CECfitBoundGrid(Z=X,C_grid=c(1,2),lambda_grid=c(10,20),grid_strategy='adaptive',
      optimization=optimizer,penalty=penalty,k0=3,B=2,
      optimizer_control=list(Rmax=1,rinit_max=1,Ninit=1,max_iter=2,Cmax=4,max_batches=1),
      adaptive_control=list(max_cells=4,max_rounds=2,max_extensions=0))
    expect_equal(calls,2L)
    expect_true(z$meta$require_reference_agreement)
    expect_false(z$meta$reference_bootstrap)
    expect_equal(length(z$measurements),2L)
    expect_true(all(z$summary$keep))
    expect_true(all(cec_grid_apply_thresholds(z)$keep))
    expect_false(any(z$summary$reference_bootstrap_applied))
    expect_equal(z$adaptive$components$trials$reference_passed,c(FALSE,TRUE))
    if(optimizer=='em_bank') {
      calls <- 0L
      rescored <- CECrescoreGrid(z,penalty=if(penalty=='none')'bic' else 'none')
      expect_true(rescored$meta$require_reference_agreement)
      expect_false(rescored$meta$reference_bootstrap)
      expect_equal(calls,2L)
    }
  }
})

test_that('fixed numeric grids use independent rejection by default and conditional opt-in remains available', {
  X <- matrix(c(-3.2,-3.1,-3,-2.9,-2.8,2.8,2.9,3,3.1,3.2),ncol=1)
  lab <- rep(1:2,each=5)
  part <- function(labels)list(labels=labels,local_labels=labels,
    stats=ca_stats(X,labels),local_stats=ca_stats(X,labels))
  banks <- c(rep(list(list(candidates=list(part(rep(1L,10))))),3),
    rep(list(list(candidates=list(part(lab)))),2))
  args <- list(Z=X,C_grid=10,lambda_grid=1,k0=3,B=2,banks=banks,
    optimization='em_bank',bootstrap_indices=list(1:10,10:1),rsi_threshold=.49,
    optimizer_control=list(Rmax=2,Ninit=1,max_iter=5))
  a <- do.call(CECfit_solver_grid,args)
  expect_true(a$meta$require_reference_agreement)
  expect_false(a$meta$reference_bootstrap)
  expect_equal(a$summary$stab_algo,1);expect_equal(a$summary$rsi,.5)
  expect_false(a$summary$keep)
  expect_false(cec_grid_apply_thresholds(a)$keep)
  expect_equal(a$summary$fail_reason,'reference_disagreement')
  b <- do.call(CECfit_solver_grid,c(args,list(require_reference_agreement=FALSE)))
  expect_true(b$summary$keep);expect_equal(b$summary$rsi,.5)
  c <- do.call(CECfit_solver_grid,c(args,list(reference_bootstrap=TRUE)))
  expect_false(c$meta$require_reference_agreement)
  expect_true(c$summary$reference_bootstrap_applied)
  expect_gt(c$summary$rsi,.95)
  expect_identical(a$winners,b$winners);expect_identical(a$selected,c$selected)
})
