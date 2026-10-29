# Conditional diagnostic: the original-data reference is fixed by construction.
# This is kept separate from CECstability's independent exhaustive averages.
CECresolve_reference_gate <- function(reference_bootstrap,require_reference_agreement=NULL) {
  if(is.null(require_reference_agreement))require_reference_agreement <- !reference_bootstrap
  if(!is.logical(require_reference_agreement)||length(require_reference_agreement)!=1L||is.na(require_reference_agreement))
    stop('require_reference_agreement must be TRUE, FALSE or NULL.')
  if(require_reference_agreement && reference_bootstrap)
    stop('require_reference_agreement requires reference_bootstrap=FALSE.')
  require_reference_agreement
}

CECapply_reference_gate <- function(g,required,algo_min) {
  g$reference_agreement_required <- required
  if(!'selection_reference_agreement'%in%names(g))g$selection_reference_agreement <- g$reference_agreement
  if(required)g$keep_harmonization <- is.finite(g$selection_reference_agreement) &
    g$selection_reference_agreement>=algo_min
  else if(!'keep_harmonization'%in%names(g))g$keep_harmonization <- TRUE
  g
}

CECreference_job_signature <- function(job,cfg)
  digest::digest(list('single_reference_descent_v1',cfg$signature,cfg$penalty,cfg$settings,job))

CECreference_bootstrap_job <- function(job,cfg) {
  signature <- CECreference_job_signature(job,cfg)
  file <- if(is.null(cfg$cache))NULL else file.path(cfg$cache,paste0(signature,'.rds'))
  if(!is.null(file) && cfg$resume && file.exists(file)) {
    saved <- readRDS(file)
    if(!identical(saved$signature,signature))stop('Reference-bootstrap checkpoint signature mismatch.')
    return(saved$result)
  }
  result <- CECwith_seed(job$seed, {
    indices <- cfg$boots[[job$rep]];Xb <- cfg$X[indices,,drop=FALSE]
    initial <- ca_labels(job$reference_labels[indices])
    backend <- CECprepare_backend_data(as.data.frame(Xb),familyType='gaussVector')
    solver <- CECsolver_environment(cfg$penalty)
    # Exactly one descent: no q_fit, random starts, or subdivision stage.
    fit <- solver$q_descent(Xb,job$C,job$lambda,initial,cfg$settings,backend)
    if(!isTRUE(fit$valid))stop('Reference-bootstrap descent failed: ',fit$error)
    params <- CECoptParam_fast_from_clusters(backend,fit$part$labels,job$lambda,job$C,Inf)
    params$phi <- NULL;fit$params <- params;fit$warm <- TRUE
    prediction <- predict_frozen(cfg$X,params,'CEC',job$projection_seed)
    list(phase='reference_bootstrap',kind='bootstrap',rep=job$rep,cell_id=job$id,
      C=job$C,lambda=job$lambda,seed=job$seed,projection_seed=job$projection_seed,
      bootstrap_indices=indices,initial_labels=initial,descent_count=1L,fit=fit,
      part=solver$q_part(cfg$X,prediction$labels,job$C,job$lambda))
  })
  if(!is.null(file))CECwrite_rds_atomic(list(signature=signature,result=result),file)
  result
}

CECreference_stability <- function(z,reference,cfg,cell,metric,algo_min,enabled=TRUE,runner=NULL,rsi_min=.95) {
  previous_runs <- z$reference_bootstrap;previous_signatures <- z$reference_check$bootstrap_signatures
  raw <- if(is.null(z$diagnostics_raw))z$diagnostics else z$diagnostics_raw
  z$diagnostics <- raw;z$reference_bootstrap <- NULL
  original <- z$original_parts
  if(is.null(original) && !is.null(z$original_labels))original <- lapply(z$original_labels,function(labels)list(labels=labels))
  if(is.null(original))original <- lapply(Filter(function(r)r$kind=='original',z$runs),`[[`,'part')
  solver <- CECsolver_environment(cfg$penalty)
  score <- function(p) {
    stats <- if(is.null(p$stats))ca_stats(cfg$X,p$labels) else p$stats
    c(list(labels=p$labels),solver$ca_profile(stats,cell$C,cell$lambda,ncol(cfg$X)))
  }
  scored <- lapply(original,score)
  best <- if(length(scored))scored[[which.min(vapply(scored,`[[`,0,'J'))]] else NULL
  target <- score(reference)
  agreement <- if(is.null(best))NA_real_ else partition_similarity(best$labels,target$labels,method=metric)
  worse_or_equal <- !is.null(best) && !ca_better(best$J,target$J,cfg$settings$tolerance)
  independent_stability_passed <- is.finite(raw$Stab_algo) && raw$Stab_algo>algo_min
  independent_rsi_passed <- is.finite(raw$RSI) && raw$RSI>rsi_min
  applied <- enabled && length(cfg$boots)>0L && independent_stability_passed && independent_rsi_passed && is.finite(raw$Stab_boot) &&
    is.finite(agreement) && agreement<algo_min && worse_or_equal
  z$diagnostics_raw <- raw
  z$reference_check <- list(applied=applied,enabled=enabled,agreement=agreement,threshold=algo_min,
    independent_stability=raw$Stab_algo,independent_stability_passed=independent_stability_passed,
    independent_rsi=raw$RSI,independent_rsi_threshold=rsi_min,independent_rsi_passed=independent_rsi_passed,
    best_original_score=if(is.null(best))NA_real_ else best$J,reference_score=target$J,
    original_not_better=worse_or_equal,tolerance=cfg$settings$tolerance,
    best_original_index=if(is.null(best))NA_integer_ else which.min(vapply(scored,`[[`,0,'J')),
    best_original_labels=if(applied)best$labels else NULL,
    reference_labels=if(applied)target$labels else NULL,reference_hash=digest::digest(target$labels))
  z$original_parts <- NULL
  if(!applied)return(z)
  jobs <- lapply(seq_along(cfg$boots),function(b)list(phase='reference_bootstrap',kind='bootstrap',
    rep=b,id=cell$id,C=cell$C,lambda=cell$lambda,reference_labels=target$labels,
    seed=ca_seed(cfg$seed,'reference_bootstrap',cell$id,b,z$reference_check$reference_hash),
    projection_seed=ca_seed(cfg$seed,'reference_projection',cell$id,b,z$reference_check$reference_hash)))
  signatures <- vapply(jobs,CECreference_job_signature,'',cfg=cfg)
  z$reference_check$bootstrap_signatures <- signatures
  runs <- if(length(previous_runs)==length(jobs) && identical(previous_signatures,signatures))previous_runs else
    if(is.null(runner))lapply(jobs,CECreference_bootstrap_job,cfg=cfg) else runner(jobs,CECreference_bootstrap_job)
  agreement <- vapply(runs,function(r)partition_similarity(target$labels,r$part$labels,method=metric),0)
  z$reference_bootstrap <- runs
  z$diagnostics <- list(Stab_algo=1,Stab_boot=mean(agreement),RSI=mean(agreement),
    original_pairs=0L,cross_pairs=length(runs),definition='conditional_reference_bootstrap',
    algorithm_agreements=numeric(),bootstrap_agreements=agreement,reference_fixed=TRUE)
  z
}

CECreference_parallel <- function(jobs,cfg,workers) {
  if(workers<=1L)return(lapply(jobs,CECreference_bootstrap_job,cfg=cfg))
  old <- Sys.getenv('R_LIBS',unset=NA_character_)
  Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
  cl <- tryCatch(parallel::makePSOCKcluster(min(workers,length(jobs))),finally={
    if(is.na(old))Sys.unsetenv('R_LIBS') else Sys.setenv(R_LIBS=old)
  })
  on.exit(parallel::stopCluster(cl),add=TRUE)
  setup <- function(paths,version) {
    .libPaths(paths);ns <- loadNamespace('CEClust')
    if(as.character(getNamespaceVersion(ns))!=version)stop('Worker loaded a different CEClust version.')
    NULL
  }
  environment(setup) <- baseenv()
  parallel::clusterCall(cl,setup,.libPaths(),as.character(utils::packageVersion('CEClust')))
  parallel::parLapplyLB(cl,jobs,CECreference_bootstrap_job,cfg=cfg)
}

CECattach_reference_summary <- function(g,measurements) {
  g$stab_algo_raw <- g$stab_boot_raw <- g$rsi_raw <- g$reference_agreement <- NA_real_
  g$reference_bootstrap_applied <- FALSE;g$diagnostic_definition <- NA_character_
  for(i in seq_len(nrow(g))) {
    direct <- measurements[[as.character(g$id[i])]]
    if(!is.null(direct)) {
      raw <- if(is.null(direct$diagnostics_raw))direct$diagnostics else direct$diagnostics_raw
      g$stab_algo_raw[i] <- raw$Stab_algo;g$stab_boot_raw[i] <- raw$Stab_boot;g$rsi_raw[i] <- raw$RSI
      g$reference_agreement[i] <- if(is.null(direct$reference_check$agreement))NA_real_ else direct$reference_check$agreement
    }
    source <- if('diagnostic_witness'%in%names(g) && !is.na(g$diagnostic_witness[i]))g$diagnostic_witness[i] else g$id[i]
    z <- measurements[[as.character(source)]]
    if(!is.null(z)) {
      g$reference_bootstrap_applied[i] <- isTRUE(z$reference_check$applied)
      g$diagnostic_definition[i] <- z$diagnostics$definition
    }
  }
  g
}
