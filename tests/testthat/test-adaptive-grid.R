adaptive_data <- function()CECwith_seed(19,cbind(c(rnorm(12,-2,.4),rnorm(12,2,.4)),rnorm(24,sd=.5)))
adaptive_args <- function()list(Z=adaptive_data(),C_grid=c(.2,1),lambda_grid=c(.6,1.4),
 k0=3,B=2,algo_seed=11,optimization='em_bank',penalty='none',grid_strategy='adaptive',
 optimizer_control=list(Rmax=3,orders=c(1,3),Ninit=1,Cmax=20,max_iter=12),
 adaptive_control=list(max_cells=120,max_rounds=3,max_extensions=2))

test_that('adaptive EM cells agree with direct bank scoring and independent diagnostics', {
  a <- do.call(CECfitBoundGrid,adaptive_args())
  expect_s3_class(a,'CEC_adaptive_grid');expect_false(a$adaptive$completed)
  pool <- CECgrid_collect(a$banks);ss <- a$summary
  for(i in seq_len(nrow(ss))) {
    truth <- CECscore_candidates(pool$parts,ss$C[i],ss$lambda[i],2,'none')
    expect_equal(ss$objective[i],truth$profile$J,tolerance=1e-12)
    expect_identical(a$parts[[a$selected[i]]]$labels,pool$parts[[truth$best]]$labels)
    for(k in seq_along(a$banks)) {
      local <- CECscore_candidates(a$banks[[k]]$candidates,ss$C[i],ss$lambda[i],2,'none','local_stats')
      expect_identical(a$parts[[a$winners[i,k]]]$labels,a$banks[[k]]$candidates[[local$best]]$labels)
    }
    labs <- lapply(a$winners[i,],function(j)a$parts[[j]]$labels)
    d <- CECstability(labs[1:3],labs[4:5])
    expect_equal(unname(unlist(ss[i,c('stab_algo','stab_boot','rsi')])),
      c(d$Stab_algo,d$Stab_boot,d$RSI),tolerance=1e-12)
    expect_identical(CECextractPartition(a,exp(.7*ss$u0[i]+.3*ss$u1[i]),.7*ss$v0[i]+.3*ss$v1[i]),
      a$parts[[a$selected[i]]]$labels)
  }
  b <- CECrescoreGrid(a,penalty='bic');expect_identical(a$banks,b$banks)
  expect_s3_class(b,'CEC_adaptive_grid')
  compact <- CECcompactGrid(a)
  for(run in a$runs_by_C)expect_false(is.unsorted(run$lambda_diag$summary$lambda))
  for(C in a$C_grid)expect_false(is.unsorted(CECselectStableLambdas(a,C=C)$summary$lambda))
  expect_identical(compact$adaptive,a$adaptive)
  expect_identical(CECextractPartition(compact,ss$C[1],ss$lambda[1]),CECextractPartition(a,ss$C[1],ss$lambda[1]))
  expect_error(CECextractPartition(a,exp(max(ss$u1)+1),1),'outside')
  path <- tempfile(fileext='.pdf');grDevices::pdf(path)
  expect_warning(CECplotGrid(a),NA);expect_warning(CECplotGrid(b),NA)
  expect_warning(CECplotPath(a,C=a$C_grid[1]),NA)
  grDevices::dev.off();expect_gt(file.info(path)$size,1000)
})

test_that('frontiers include equal-order partition changes and refine both sides', {
  mesh <- ceca_mesh(c(1,2),c(1,2))
  a <- list(labels=c(1L,1L,2L,2L),K=2L,saturation=0)
  b <- list(labels=c(1L,2L,1L,2L),K=2L,saturation=0)
  f <- ceca_frontiers(mesh,list(a,b,a,b),0,TRUE)
  expect_equal(length(f$cells),4L)
  expect_equal(nrow(ceca_frontiers(mesh,list(a,b,a,b),0,FALSE)$edges),0L)
  refined <- ceca_refine(mesh,f$cells)$mesh
  expect_equal(nrow(refined$leaves),16L)
  area <- function(g)sum((g$u1-g$u0)*(g$v1-g$v0))
  expect_equal(area(mesh$leaves),area(refined$leaves))
  for(i in seq_len(nrow(refined$leaves)))for(j in seq_len(nrow(refined$leaves)))if(i<j) {
    a <- refined$leaves[i,];b <- refined$leaves[j,]
    expect_false(min(a$u1,b$u1)>max(a$u0,b$u0)+1e-12 && min(a$v1,b$v1)>max(a$v0,b$v0)+1e-12)
  }
})

test_that('right stopping uses matching rows and two unsaturated extensions', {
  state <- list(mesh=ceca_mesh(c(1,2),c(1,2)),parts=list())
  fill <- function(state) {
    g <- state$mesh$leaves
    for(i in seq_len(nrow(g)))state$parts[[as.character(g$id[i])]] <-
      list(labels=c(1L,1L,2L,2L),K=2L,saturation=if(g$lambda[i]>1.5)1 else 0)
    state
  }
  state <- fill(state);memory <- ceca_right_memory(state);cfg <- list(sat_max=0)
  expect_false(ceca_right_status(state,cfg,memory)$observed['right'])
  state$mesh <- ceca_expand(state$mesh,'right')$mesh;state <- fill(state)
  memory <- ceca_right_refresh(state,memory,TRUE)
  expect_false(ceca_right_status(state,cfg,memory)$observed['right'])
  state$mesh <- ceca_expand(state$mesh,'right')$mesh;state <- fill(state)
  memory <- ceca_right_refresh(state,memory,TRUE)
  z <- ceca_right_status(state,cfg,memory)
  expect_true(unname(z$observed['right']));expect_identical(z$right_stopping$reason,'rowwise_two_extensions_REO_stable')
  state$mesh <- ceca_refine(state$mesh,ceca_border(state$mesh)$right)$mesh;state <- fill(state)
  memory <- ceca_right_refresh(state,memory)
  expect_length(memory$panels,1L);expect_gt(memory$resets,0)
})

test_that('limits are reported and no unmeasured cell passes the filter', {
  args <- adaptive_args();args$rsi_policy <- 'nontrivial_unsaturated'
  args$adaptive_control <- list(max_cells=40,max_rounds=2,max_extensions=0)
  a <- do.call(CECfitBoundGrid,args);s <- a$summary
  omitted <- s$REO==1 | s$sat>0
  expect_true(any(omitted));expect_true(all(is.na(s$rsi[omitted])))
  expect_false(any(s$keep[omitted]));expect_false(a$adaptive$completed)
  expect_equal(sum(a$adaptive$extensions),0)
  expect_false(any(cec_grid_apply_thresholds(a,0,0,1)$keep[omitted]))
  args$adaptive_control$trivial_quadrants <- TRUE
  expect_error(do.call(CECfitBoundGrid,args),'anchors')
  args <- adaptive_args();args$optimizer_control$Cmax <- 1
  b <- do.call(CECfitBoundGrid,args)
  expect_false(b$adaptive$completed)
  expect_true(all(b$summary$C<=1+1e-12))
})

test_that('quadrant coverage never manufactures diagnostic values', {
  args <- adaptive_args();args$C_grid<-c(.01,.02);args$lambda_grid<-c(10,20)
  args$adaptive_control<-list(max_cells=50,max_rounds=2,max_extensions=0,trivial_quadrants=TRUE)
  a <- do.call(CECfitBoundGrid,args)
  waived <- a$summary$measurement_status=='quadrant_unmeasured'
  expect_true(any(waived));expect_true(all(is.na(a$summary$rsi[waived])))
  expect_false(any(a$summary$keep[waived]))
  expect_true(all(!is.na(a$summary$quadrant_witness[waived])))
  audit <- a$adaptive$screening$audit_state
  expect_equal(length(audit$audit_ids),ceiling(length(audit$candidate_ids)*.1-1e-12))
  expect_true(all(audit$audit_ids%in%audit$candidate_ids))
})

test_that('adaptive CEM supports both penalties and preserves repetitions through harmony', {
  args <- adaptive_args();args$optimization <- 'cem';args$k0 <- 2;args$B <- 1
  args$optimizer_control <- list(Rmax=2,Ninit=1,Cmax=20,max_iter=8,max_batches=1)
  args$adaptive_control <- list(max_cells=60,max_rounds=2,max_extensions=0)
  for(penalty in c('none','bic')) {
    args$penalty <- penalty;a <- do.call(CECfitBoundGrid,args)
    expect_true(all(a$summary$A_harm>=0 & a$summary$A_harm<=1))
    for(h in a$adaptive$harmonization)if(nrow(h$updates))expect_true(all(h$updates$after<h$updates$before))
    for(i in seq_len(nrow(a$summary))) {
      labs <- lapply(a$winners[i,],function(j)a$parts[[j]]$labels)
      d <- CECstability(labs[1:2],labs[3])
      expect_equal(unname(unlist(a$summary[i,c('stab_algo','stab_boot','rsi')])),
        c(d$Stab_algo,d$Stab_boot,d$RSI),tolerance=1e-12)
    }
  }
})

test_that('explicit order searches and interrupted EM starts resume reproducibly', {
  X <- adaptive_data();ctl <- CECoptimization_control(list(seed=7,orders=c(1,3),Ninit=2,max_iter=10),24,2,2,10)
  cache <- tempfile('starts');dir.create(cache)
  set.seed(4);before <- .Random.seed
  a <- CECfit_em_bank(X,ctl,checkpoint_dir=cache,resume=TRUE)
  expect_identical(.Random.seed,before);expect_equal(length(a$starts),4L)
  expect_identical(a,CECfit_em_bank(X,ctl))
  expect_equal(vapply(a$starts,`[[`,0L,'r'),c(1L,1L,3L,3L))
  last <- file.path(cache,'start_0003_002.rds')
  expect_true(file.rename(last,paste0(last,'.saved')))
  b <- CECfit_em_bank(X,ctl,checkpoint_dir=cache,resume=TRUE)
  expect_identical(a,b);expect_identical(.Random.seed,before)
  expect_error(CECoptimization_control(list(orders=c(1,25)),24,2,1,10),'between')
  expect_error(CECoptimization_control(list(Rmax=2,orders=c(1,3)),24,2,1,10),'Rmax')
  expect_error(CECadaptive_initial(X,c(1,2,5),c(1,2)),'equally spaced')
})

test_that('later CEM harmony starts locally and leaves no improving neighbour transfer', {
  args <- adaptive_args()
  args$optimization <- 'cem';args$k0 <- 2;args$B <- 1
  args$C_grid <- c(.3,2);args$lambda_grid <- c(.8,1.4)
  args$optimizer_control <- list(Rmax=2,Ninit=1,Cmax=20,max_iter=10,max_batches=1)
  args$adaptive_control <- list(max_cells=80,max_rounds=2,max_extensions=1)
  for(penalty in c('none','bic')) {
    args$penalty <- penalty
    a <- do.call(CECfitBoundGrid,args)
    h <- a$adaptive$harmonization
    expect_gt(length(h),1L)
    expect_equal(h[[1]]$seed_cells,4L)
    expect_true(any(vapply(h[-1],function(x)x$seed_cells<x$evaluated_cells,TRUE)))
    for(x in h)if(nrow(x$updates))expect_true(all(x$updates$after<x$updates$before))
    profile <- CECsolver_environment(penalty)$ca_profile
    s <- a$summary
    state <- list(mesh=a$adaptive$mesh,parts=list())
    for(i in seq_len(nrow(s))) {
      p <- a$parts[[a$selected[i]]]
      state$parts[[as.character(s$id[i])]] <-
        c(profile(p$stats,s$C[i],s$lambda[i],ncol(a$Z)),p)
    }
    settings <- q_settings();settings$Nloop <- 10L
    cfg <- list(X=as.matrix(a$Z),seed=a$meta$algo_seed,settings=settings)
    fun <- ceca_harmonize;environment(fun) <- CECsolver_environment(penalty)
    audit <- fun(state,cfg,tag='full_audit_after_local')
    expect_equal(audit$changes,0L)
  }
})

test_that('rescoring a common bank never starts another EM fit', {
  a <- do.call(CECfitBoundGrid,adaptive_args())
  local_mocked_bindings(CECfit_em_bank=function(...)stop('Unexpected EM fitting'),
    .package='CEClust')
  b <- CECrescoreGrid(a,penalty='bic')
  expect_identical(a$banks,b$banks)
  fixed <- CECrescoreGrid(a,C_grid=c(.2,1),lambda_grid=c(.6,1.4),
    grid_strategy='fixed',penalty='bic')
  expect_identical(a$banks,fixed$banks)
})
