staged_args <- function()list(Z=CECwith_seed(19,cbind(c(rnorm(12,-2,.4),rnorm(12,2,.4)),rnorm(24,sd=.5))),C_grid=c(.001,.004),lambda_grid=c(.6,1.4),
  k0=3,B=2,algo_seed=11,optimization='cem',grid_strategy='adaptive',
  optimizer_control=list(Rmax=2,Ninit=1,max_iter=8,max_batches=1),
  adaptive_control=list(max_cells=60,max_rounds=2,max_extensions=0,diagnostic_strategy='cells',reference_bootstrap=FALSE))

test_that('a saturated pilot performs exactly one fit per cell and no campaigns', {
  calls <- list();worker <- CECadaptive_cem_job
  local_mocked_bindings(CECadaptive_cem_job=function(job,cfg) {
    calls[[length(calls)+1L]] <<- job
    worker(job,cfg)
  },.package='CEClust')
  a <- do.call(CECfitBoundGrid,staged_args())
  expect_identical(a$meta$rsi_policy,'saturation_band')
  expect_true(all(a$summary$sat>.1));expect_length(a$measurements,0L)
  expect_length(calls,4L)
  expect_true(all(vapply(calls,function(j)j$phase=='pilot' && j$kind=='original',TRUE)))
  expect_true(all(is.na(a$winners)));expect_true(all(is.na(a$summary$rsi)))
  expect_false(any(a$summary$keep));expect_length(a$adaptive$pilot_fits,4L)
  expect_true(all(a$adaptive$history$measured_cells==0L))
  expect_true(all(a$executions$phase=='pilot'))
})

test_that('campaigns follow pilot completion, are selective and remain independent', {
  args <- staged_args();args$C_grid <- c(1,5)
  for(penalty in c('none','bic')) {
    args$penalty <- penalty
    a <- do.call(CECfitBoundGrid,args)
    h <- a$adaptive$history;end <- which(h$stage=='pilot_complete')
    expect_true(all(h$measured_cells[seq_len(end)]==0))
    expect_gt(length(a$measurements),0L)
    expect_true(all(vapply(a$measurements,function(z)z$saturation_at_measurement<=.1,TRUE)))
    for(z in a$measurements) {
      expect_length(z$runs,5L)
      expect_true(all(vapply(z$runs,function(r)!r$fit$warm,TRUE)))
      original <- lapply(z$runs[1:3],function(r)r$part$labels)
      bootstrap <- lapply(z$runs[4:5],function(r)r$part$labels)
      expect_identical(z$diagnostics,CECstability(original,bootstrap))
      expect_equal(z$diagnostics$original_pairs,3L);expect_equal(z$diagnostics$cross_pairs,6L)
      i <- match(z$cell_id,a$summary$id)
      expect_equal(a$summary$A_harm[i],ca_sim(z$best$labels,a$parts[[a$selected[i]]]$labels))
      for(r in z$runs) {
        expect_identical(r$seed,ca_seed(a$meta$algo_seed,'measurement',z$cell_id,r$kind,r$rep))
        if(r$kind=='bootstrap') {
          pred <- predict_frozen(as.matrix(a$Z),r$fit$params,'CEC',
            ca_seed(a$meta$algo_seed,'projection',z$cell_id,r$rep))
          expect_identical(r$part$labels,ca_labels(pred$labels))
        }
      }
    }
    unmeasured <- !as.character(a$summary$id)%in%names(a$measurements)
    expect_true(all(is.na(a$winners[unmeasured,,drop=FALSE])))
    expect_true(all(is.na(a$summary$rsi[unmeasured])))
    expect_false(any(a$summary$keep[unmeasured]))
    measured_jobs <- subset(a$executions,phase=='measurement')
    expect_false(anyDuplicated(measured_jobs[,c('cell_id','kind','rep')])>0L)
    expect_true(all(table(measured_jobs$cell_id)==5L))
  }
})

test_that('child pilots receive a parent start, never a previous measurement run', {
  args <- staged_args();args$C_grid <- c(1,5);args$lambda_grid <- c(.6,1.4)
  args$adaptive_control$max_cells <- 100L
  a <- do.call(CECfitBoundGrid,args)
  children <- Filter(function(z)z$fit$warm,a$adaptive$pilot_fits)
  expect_gt(length(children),0L)
  expect_true(all(vapply(children,function(z)any(z$fit$starts$phase=='warm'),TRUE)))
  archive <- a$adaptive$mesh$archive
  expect_gt(nrow(archive),0L)
  expect_false(any(as.character(archive$id)%in%names(a$measurements)))
})

test_that('completed numerical work resumes without fitting or changing independent winners', {
  args <- staged_args();args$C_grid <- c(1,5)
  args$checkpoint_dir <- tempfile('staged_resume');args$resume <- TRUE
  a <- do.call(CECfitBoundGrid,args)
  expect_gt(length(list.files(args$checkpoint_dir,recursive=TRUE)),4L)
  local_mocked_bindings(q_fit=function(...)stop('Unexpected repeated fitting'),
    q_descent=function(...)stop('Unexpected repeated descent'),.package='CEClust')
  b <- do.call(CECfitBoundGrid,args)
  expect_identical(a$winners,b$winners);expect_identical(a$selected,b$selected)
  expect_identical(a$measurements,b$measurements)
  expect_identical(a$adaptive$mesh,b$adaptive$mesh)
  expect_equal(a$summary,b$summary)
  expect_true(all(b$executions$checkpoint_reused))
})

test_that('twenty independent runs use 190 and 400 comparisons, including trivial cells', {
  args <- staged_args();args$C_grid <- c(1,2);args$lambda_grid <- c(10,20)
  args$k0 <- args$B <- 20L
  args$optimizer_control <- list(Rmax=1,rinit_max=1,Ninit=1,max_iter=2,max_batches=1)
  a <- do.call(CECfitBoundGrid,args)
  expect_true(all(a$summary$REO==1L));expect_equal(length(a$measurements),4L)
  for(z in a$measurements) {
    expect_equal(z$diagnostics$original_pairs,190L)
    expect_equal(z$diagnostics$cross_pairs,400L)
    expect_equal(z$diagnostics$RSI,1)
  }
  expect_equal(sum(a$executions$phase=='pilot'),4L)
  expect_equal(sum(a$executions$phase=='measurement'),160L)
})

test_that('measurement improvements reopen coverage without replacing completed campaigns', {
  # Controlled objective values isolate orchestration from numerical fitting.
  # Only cell 1 is eligible initially; its improvement makes cell 2 eligible.
  local_mocked_bindings(CECadaptive_cem_job=function(job,cfg) {
    n <- nrow(cfg$X)
    labels <- if(job$phase=='pilot' || job$rep==2L)rep(1L,n) else rep(1:2,each=n/2)
    p <- CECsolver_environment(cfg$penalty)$q_part(cfg$X,labels,job$C,job$lambda)
    p$J <- if(job$phase=='pilot')100 else -10*job$id
    p$saturation <- if(job$phase=='pilot' && job$id!=1L)1 else 0
    list(phase=job$phase,kind=job$kind,rep=job$rep,cell_id=job$id,C=job$C,
      lambda=job$lambda,seed=job$seed,part=p,fit=list(warm=!is.null(job$warm)))
  },ceca_harmonize=function(state,cfg,seed_ids=NULL,tag='') {
    if(tag=='measurement_1') {
      p <- state$parts[['1']]
      p$labels <- rep(1:2,length.out=nrow(cfg$X));p$stats <- ca_stats(cfg$X,p$labels)
      p$K <- 2L;p$J <- -11;state$parts[['1']] <- p
      state$parts[['2']]$saturation <- 0
    }
    list(state=state,time=c(user=0,system=0,elapsed=0),changes=0L,
      seed_cells=length(seed_ids),updates=data.frame())
  },.package='CEClust')
  args <- staged_args();args$C_grid <- c(1,2);args$adaptive_control$max_cells <- 4L
  a <- do.call(CECfitBoundGrid,args)
  expect_setequal(names(a$measurements),c('1','2'))
  expect_equal(a$measurements[['1']]$wave,1L);expect_equal(a$measurements[['2']]$wave,2L)
  jobs <- subset(a$executions,phase=='measurement')
  expect_equal(nrow(jobs),10L);expect_false(anyDuplicated(jobs[,c('cell_id','kind','rep')])>0L)
  z <- a$measurements[['1']]
  d <- CECstability(lapply(z$runs[1:3],function(r)r$part$labels),
    lapply(z$runs[4:5],function(r)r$part$labels))
  expect_identical(z$diagnostics,d)
  i <- match(1L,a$summary$id)
  expect_lt(a$summary$A_harm[i],1)
  expect_true(a$adaptive$final_refinement_pending)
  expect_false(a$adaptive$completed)
})

test_that('rounding of a mesh centre does not reject a valid density ceiling', {
  ctl <- CECoptimization_control(list(Cmax=.004),24,.004+2e-18,1,8)
  expect_gte(ctl$Cmax,.004+2e-18)
  expect_error(CECoptimization_control(list(Cmax=.0039),24,.004,1,8),'Cmax')
})

test_that('EM ceiling restrictions leave interior frontier refinements available', {
  mesh <- ceca_mesh(c(1,2),c(1,2))
  plan <- CECadaptive_refinable(mesh,seq_len(nrow(mesh$leaves)),2)
  expect_equal(length(plan$indices),2L);expect_equal(length(plan$blocked),2L)
  refined <- ceca_refine(mesh,plan$indices)$mesh
  expect_true(all(refined$leaves$C<=2*(1+1e-12)))
  expect_true(all(mesh$leaves$id[plan$blocked]%in%refined$leaves$id))
})

test_that('trivial quadrant exemptions avoid whole CEM campaigns and retain NA', {
  args <- staged_args();args$C_grid <- c(1,2);args$lambda_grid <- c(10,20)
  args$adaptive_control$trivial_quadrants <- TRUE
  a <- do.call(CECfitBoundGrid,args)
  waived <- a$summary$measurement_status=='quadrant_unmeasured'
  expect_true(any(waived));expect_lt(length(a$measurements),nrow(a$summary))
  expect_false(any(as.character(a$summary$id[waived])%in%names(a$measurements)))
  expect_true(all(is.na(a$winners[waived,,drop=FALSE])))
  expect_true(all(is.na(a$summary$rsi[waived])));expect_false(any(a$summary$keep[waived]))
  audit <- a$adaptive$screening$audit_state
  expect_equal(length(audit$audit_ids),ceiling(.1*length(audit$candidate_ids)-1e-12))
})
