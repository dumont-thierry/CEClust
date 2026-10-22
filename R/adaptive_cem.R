# Staged adaptive CEM. Pilot fits never use a bootstrap sample. Measurement
# executions have their own seeds and never receive a pilot or neighbour start.
CECadaptive_cem_settings <- function(control) {
  s <- q_settings()
  s$rinit <- control$Rmax; s$shots <- control$Ninit; s$Nloop <- control$max_iter
  s$rinit_max <- max(control$rinit_max,s$rinit)
  s$min_child <- control$min_child; s$confirmations <- control$confirmations
  s$max_batches <- control$max_batches; s$tau <- control$saturation_stop
  s$tolerance <- control$tolerance
  s
}

# Jobs are checkpointed individually, including partially completed campaigns.
CECadaptive_cem_job <- function(job,cfg) {
  signature <- digest::digest(list(cfg$signature,job))
  file <- if(is.null(cfg$cache))NULL else file.path(cfg$cache,paste0(signature,'.rds'))
  if(!is.null(file) && cfg$resume && file.exists(file)) {
    saved <- readRDS(file)
    if(!identical(saved$signature,signature))stop('CEM cell checkpoint signature mismatch.')
    return(saved$result)
  }
  result <- CECwith_seed(job$seed, {
    Xfit <- if(job$kind=='bootstrap')cfg$X[cfg$boots[[job$rep]],,drop=FALSE] else cfg$X
    warm <- if(job$phase=='pilot')job$warm else NULL
    f <- CECsolver_environment(cfg$penalty)$q_fit(Xfit,job$C,job$lambda,
      job$seed,cfg$settings,warm=warm)
    if(f$status=='fit_failure')stop('A CEM execution failed; no stability is assigned.')
    p <- f$part
    if(job$kind=='bootstrap') {
      prediction <- predict_frozen(cfg$X,f$params,'CEC',
        ca_seed(cfg$seed,'projection',job$id,job$rep))
      p <- CECsolver_environment(cfg$penalty)$q_part(cfg$X,prediction$labels,job$C,job$lambda)
    }
    list(phase=job$phase,kind=job$kind,rep=job$rep,cell_id=job$id,
      C=job$C,lambda=job$lambda,seed=job$seed,fit=f,part=p)
  })
  if(!is.null(file))CECwrite_rds_atomic(list(signature=signature,result=result),file)
  result
}

CECadaptive_cem_eligible <- function(g,policy,maximum) {
  switch(policy,all=rep(TRUE,nrow(g)),unsaturated=g$saturation==0,
    nontrivial_unsaturated=g$REO>1 & g$saturation==0,
    saturation_band=g$saturation<=maximum)
}

CECfit_adaptive_cem <- function(X,initial,o,k0,B,r0,Nshots_fresh,Nloop,algo_seed,n_cores,
  penalty,optimizer_control,rsi_policy,rsi_saturation_max,stab_algo_threshold,rsi_threshold,
  sat_threshold,partition_metric,bootstrap_indices,checkpoint_dir,resume,force_recompute,
  save_results,output_dir,verbose) {
  started <- proc.time()[3]
  k0 <- CECinteger(k0,'k0'); B <- CECinteger(B,'B',0L)
  algo_seed <- CECinteger(algo_seed,'algo_seed',0L); n_cores <- CECinteger(n_cores,'n_cores')
  for(x in list(rsi_saturation_max,stab_algo_threshold,rsi_threshold,sat_threshold))
    if(!is.numeric(x)||length(x)!=1L||!is.finite(x)||x<0)stop('Diagnostic thresholds must be finite nonnegative numbers.')
  if(rsi_saturation_max>1 || sat_threshold>1 || stab_algo_threshold>1)
    stop('Saturation and agreement thresholds must not exceed one.')
  control <- CECoptimization_control(optimizer_control,nrow(X),max(initial$C),Nshots_fresh,Nloop,r0)
  if(!is.null(control$seed))algo_seed <- control$seed
  boots <- if(is.null(bootstrap_indices))CECwith_seed(ca_seed(algo_seed,'bootstrap_indices'),
    replicate(B,sample.int(nrow(X),nrow(X),replace=TRUE),simplify=FALSE)) else bootstrap_indices
  if(!is.list(boots)||length(boots)!=B||any(vapply(boots,function(x)length(x)!=nrow(X)||anyNA(x)||
      any(x!=floor(x))||any(x<1|x>nrow(X)),TRUE)))stop('Invalid bootstrap_indices.')
  signature <- digest::digest(list(protocol='staged_cem_v1',version=as.character(utils::packageVersion('CEClust')),
    X=X,boots=boots,L=k0,B=B,control=control,seed=algo_seed,initial=initial,adaptive=o,
    penalty=penalty,policy=rsi_policy,measurement=rsi_saturation_max,
    filters=c(stab_algo_threshold,rsi_threshold,sat_threshold),metric=partition_metric))
  cache <- NULL
  if(!isFALSE(checkpoint_dir) && !is.null(checkpoint_dir)) {
    if(!is.character(checkpoint_dir)||length(checkpoint_dir)!=1L)stop('checkpoint_dir must be FALSE or a directory path.')
    cache <- file.path(checkpoint_dir,paste0('staged_cem_',signature))
    dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  }
  cfg <- list(X=X,boots=boots,seed=algo_seed,settings=CECadaptive_cem_settings(control),
    penalty=penalty,signature=signature,cache=cache,resume=isTRUE(resume)&&!isTRUE(force_recompute))
  cluster <- NULL
  if(n_cores>1L) {
    old <- Sys.getenv('R_LIBS',unset=NA_character_)
    Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
    cluster <- tryCatch(parallel::makePSOCKcluster(n_cores),finally={
      if(is.na(old))Sys.unsetenv('R_LIBS') else Sys.setenv(R_LIBS=old)
    })
    on.exit(parallel::stopCluster(cluster),add=TRUE)
    setup <- function(paths,version) {
      .libPaths(paths);ns <- loadNamespace('CEClust')
      if(as.character(getNamespaceVersion(ns))!=version)stop('Worker loaded a different CEClust version.')
      NULL
    }
    environment(setup) <- baseenv()
    parallel::clusterCall(cluster,setup,.libPaths(),as.character(utils::packageVersion('CEClust')))
  }
  executions <- list();job_elapsed <- 0;agreement_elapsed <- 0;harmony_elapsed <- 0
  run_jobs <- function(jobs) {
    if(!length(jobs))return(list())
    for(job in jobs) {
      key <- digest::digest(list(cfg$signature,job))
      reused <- !is.null(cache) && cfg$resume && file.exists(file.path(cache,paste0(key,'.rds')))
      executions[[length(executions)+1L]] <<- data.frame(phase=job$phase,kind=job$kind,
        cell_id=job$id,rep=job$rep,seed=job$seed,C=job$C,lambda=job$lambda,
        warm=job$phase=='pilot' && !is.null(job$warm),checkpoint_reused=reused)
    }
    tm <- proc.time()[3]
    result <- if(is.null(cluster))lapply(jobs,CECadaptive_cem_job,cfg=cfg) else
      parallel::parLapplyLB(cluster,jobs,CECadaptive_cem_job,cfg=cfg)
    job_elapsed <<- job_elapsed+unname(proc.time()[3]-tm)
    result
  }
  state <- list(mesh=ceca_mesh(initial$C,initial$lambda),parts=list())
  history <- restrictions <- harmonies <- refinelog <- list()
  pilots <- measurements <- metadata <- list();audit <- NULL
  right <- NULL;round <- 1L;phase <- 'pilot';measurement_wave <- 0L
  counts <- setNames(integer(4),c('left','right','lower','upper'))
  blocked <- setNames(rep(FALSE,4),names(counts))
  settings <- list(sat_max=sat_threshold,sat_max2=o$frontier_saturation,max_cells=o$max_cells,
    refinement=o[c('max_rounds','max_depth','patience','position_tolerance','distinguish_partitions')])
  screen_cfg <- list(seed=algo_seed,workers=1L,
    sat_max3=if(rsi_policy=='all')1 else if(rsi_policy=='saturation_band')rsi_saturation_max else 0,
    RSI_min=rsi_threshold,screening=list(trivial_quadrants=o$trivial_quadrants,
      anchor_rule='rsi_and_algo',audit_fraction=o$audit_fraction,
      stability_min=max(stab_algo_threshold,o$anchor_stability),audit_scheme='shared_persistent'))
  borders <- function(append=FALSE) {
    right <<- ceca_right_refresh(state,right,append)
    ceca_right_status(state,settings,right)
  }
  snapshot <- function(stage,status=borders()) {
    g <- state$mesh$leaves
    history[[length(history)+1L]] <<- data.frame(step=length(history)+1L,stage=stage,phase=phase,
      round=round,cells=nrow(g),pilot_executions=length(pilots),measured_cells=length(measurements),
      Cmin=min(g$C),Cmax=max(g$C),lambda_min=min(g$lambda),lambda_max=max(g$lambda),
      left=unname(status$observed['left']),right=unname(status$observed['right']),
      lower=unname(status$observed['lower']),upper=unname(status$observed['upper']),
      right_reason=status$right_stopping$reason,elapsed=unname(proc.time()[3]-started))
    if(isTRUE(verbose))message(stage,': ',nrow(g),' cells; ',length(measurements),' measured.')
  }
  harmonize <- function(ids=NULL,tag) {
    if(!o$harmonize)return(invisible(NULL))
    tm <- proc.time()[3]
    # Replay completed propagation without rerunning its numerical descents.
    key <- digest::digest(list(signature,state,ids,tag))
    file <- if(is.null(cache))NULL else file.path(cache,paste0('harmony_',key,'.rds'))
    if(!is.null(file) && cfg$resume && file.exists(file))h <- readRDS(file) else {
      fun <- ceca_harmonize;environment(fun) <- CECsolver_environment(penalty)
      h <- CECwith_seed(ca_seed(algo_seed,'adaptive_harmony',tag),fun(state,cfg,seed_ids=ids,tag=tag))
      if(!is.null(file))CECwrite_rds_atomic(h,file)
    }
    h$evaluated_cells <- nrow(state$mesh$leaves);h$initial_cells <- ids
    state <<- h$state;harmonies[[length(harmonies)+1L]] <<- h[names(h)!='state']
    harmony_elapsed <<- harmony_elapsed+unname(proc.time()[3]-tm)
  }
  evaluate <- function(first=FALSE) {
    g <- state$mesh$leaves;fresh <- which(!as.character(g$id)%in%names(state$parts))
    jobs <- lapply(fresh,function(i) {
      parent <- state$parts[[as.character(g$parent[i])]]
      list(phase='pilot',kind='original',rep=1L,id=g$id[i],C=g$C[i],lambda=g$lambda[i],
        warm=if(is.null(parent))NULL else parent$labels,seed=ca_seed(algo_seed,'pilot',g$id[i]))
    })
    fits <- run_jobs(jobs)
    for(j in seq_along(fresh)) {
      id <- as.character(g$id[fresh[j]])
      pilots[[id]] <<- fits[[j]];state$parts[[id]] <<- fits[[j]]$part
    }
    if(length(fresh))harmonize(if(first)NULL else g$id[fresh],paste0(phase,'_',length(pilots)))
  }
  extend <- function() {
    repeat {
      m <- borders();allowed <- m$needs_extension & counts<o$max_extensions & !blocked
      if(!any(allowed))break
      side <- c('upper','left','lower','right')[which(allowed[c('upper','left','lower','right')])[1]]
      plan <- ceca_expand(state$mesh,side)
      reason <- if(any(!is.finite(plan$mesh$leaves$C)) || any(plan$mesh$leaves$C<=0) ||
        any(!is.finite(plan$mesh$leaves$lambda)) || any(plan$mesh$leaves$lambda<=0))'numeric_limit' else
          if(nrow(plan$mesh$leaves)>o$max_cells)'cell_budget' else ''
      if(nzchar(reason)) {
        restrictions[[length(restrictions)+1L]] <<- list(side=side,reason=reason,round=round)
        blocked[side] <<- TRUE;next
      }
      counts[side] <<- counts[side]+1L;state$mesh <<- plan$mesh;evaluate()
      snapshot(paste0('extension_',side),borders(side=='right'))
    }
  }
  evaluate(TRUE);snapshot('initial');extend();memory <- ceca_refinement_memory(state)
  scale <- c(diff(log(initial$C))[1],diff(initial$lambda)[1])
  repeat {
    pick <- ceca_select_refinement(state,memory,settings)
    if(pick$status!='refine'){stop_reason <- pick$status;break}
    if(round>=o$max_rounds){stop_reason <- 'round_budget';break}
    before <- state;ids <- state$mesh$leaves$id[pick$indices]
    state$mesh <- ceca_refine(state$mesh,pick$indices)$mesh;round <- round+1L;evaluate()
    update <- ceca_refinement_update(before,state,memory,ids,settings,scale)
    memory <- update$memory;refinelog[[length(refinelog)+1L]] <- cbind(round=round,update$diagnostics)
    snapshot('refinement');extend();memory <- ceca_wake_changed(state,memory,o$frontier_saturation)
  }
  pilot <- state;snapshot('pilot_complete')
  pilot_elapsed <- unname(proc.time()[3]-started)
  if(!is.null(cache))CECwrite_rds_atomic(list(state=pilot,memory=memory,stop_reason=stop_reason),file.path(cache,'pilot_complete.rds'))
  phase <- 'measurement';diagnostic_started <- proc.time()[3]
  repeat {
    plan <- ceca_screen_plan(state,screen_cfg,metadata,audit);audit <- plan$audit_state
    g <- ceca_diagram(state);eligible <- CECadaptive_cem_eligible(g,rsi_policy,rsi_saturation_max)
    todo <- plan$next_indices;todo <- todo[eligible[todo]]
    if(!length(todo))todo <- head(which(eligible & !as.character(g$id)%in%names(metadata) &
      plan$table$status=='pending' & !plan$table$awaiting_quadrant_audit),screen_cfg$workers)
    if(!length(todo))break
    measurement_wave <- measurement_wave+1L
    jobs <- unlist(lapply(todo,function(i)lapply(seq_len(k0+B),function(k) {
      kind <- if(k<=k0)'original' else 'bootstrap';rep <- if(k<=k0)k else k-k0
      list(phase='measurement',kind=kind,rep=rep,id=g$id[i],C=g$C[i],lambda=g$lambda[i],
        seed=ca_seed(algo_seed,'measurement',g$id[i],kind,rep))
    })),recursive=FALSE)
    fits <- run_jobs(jobs);changed <- integer()
    for(j in seq_along(todo)) {
      i <- todo[j];id <- as.character(g$id[i]);runs <- fits[(j-1L)*(k0+B)+seq_len(k0+B)]
      parts <- lapply(runs,`[[`,'part');labels <- lapply(parts,`[[`,'labels')
      tm <- proc.time()[3]
      d <- CECstability(labels[seq_len(k0)],labels[k0+seq_len(B)],partition_metric)
      agreement_elapsed <- agreement_elapsed+unname(proc.time()[3]-tm)
      best <- parts[[which.min(vapply(parts,`[[`,0,'J'))]]
      values <- data.frame(method='exhaustif',Stab_algo=d$Stab_algo,Stab_boot=d$Stab_boot,RSI=d$RSI)
      qualifies <- is.finite(d$RSI) && d$RSI>rsi_threshold && is.finite(d$Stab_algo) &&
        d$Stab_algo>screen_cfg$screening$stability_min
      metadata[[id]] <- list(best=best,values=values,qualifies=qualifies)
      measurements[[id]] <- list(cell_id=g$id[i],C=g$C[i],lambda=g$lambda[i],
        saturation_at_measurement=g$saturation[i],wave=measurement_wave,runs=runs,best=best,diagnostics=d)
      if(ca_better(best$J,state$parts[[id]]$J,cfg$settings$tolerance)) {
        state$parts[[id]] <- best;changed <- c(changed,g$id[i])
      }
    }
    if(length(changed))harmonize(changed,paste0('measurement_',measurement_wave))
    snapshot('measurement_update');extend()
    # Already measured cells remain frozen. New saturation and domain geometry
    # determine only whether further cells need independent campaigns.
  }
  plan <- ceca_screen_plan(state,screen_cfg,metadata,audit);audit <- plan$audit_state
  g <- ceca_diagram(state);eligible <- CECadaptive_cem_eligible(g,rsi_policy,rsi_saturation_max)
  pending <- eligible & plan$table$status=='pending'
  if(any(pending))stop('Internal error: eligible CEM cells remain unmeasured without an exemption.')
  final_memory <- ceca_wake_changed(state,memory,o$frontier_saturation)
  final_pick <- ceca_select_refinement(state,final_memory,settings)
  refinement_pending <- !final_pick$status%in%c('no_frontiers','adaptive_stable')
  snapshot('measurement_complete')
  pool <- list();keys <- new.env(parent=emptyenv())
  register <- function(p) {
    key <- digest::digest(p$labels)
    if(exists(key,keys,inherits=FALSE)) {
      id <- get(key,keys);if(!identical(pool[[id]]$labels,p$labels))stop('Partition hash collision.')
    } else {
      id <- length(pool)+1L;pool[[id]] <<- list(labels=p$labels,stats=p$stats);assign(key,id,keys)
    }
    id
  }
  ids <- as.character(g$id);p <- ceca_current(state)
  final <- vapply(p,register,0L);winners <- matrix(NA_integer_,nrow(g),k0+B)
  g$A_harm <- g$stab_algo <- g$stab_boot <- g$rsi <- NA_real_
  for(i in which(ids%in%names(measurements))) {
    z <- measurements[[ids[i]]]
    winners[i,] <- vapply(z$runs,function(r)register(r$part),0L)
    g$stab_algo[i] <- z$diagnostics$Stab_algo;g$stab_boot[i] <- z$diagnostics$Stab_boot;g$rsi[i] <- z$diagnostics$RSI
    g$A_harm[i] <- ca_sim(z$best$labels,p[[i]]$labels)
  }
  g$sat <- g$saturation;g$objective <- vapply(p,`[[`,0,'J');g$Hphi <- g$objective/g$lambda
  g$entropy <- vapply(p,`[[`,0,'J_raw');g$penalty_value <- vapply(p,`[[`,0,'BIC_penalty_relative')
  g$n_sat_clusters <- vapply(p,function(p)sum(p$saturated),0L);g$any_bound_saturation <- g$n_sat_clusters>0L
  g$measured <- is.finite(g$stab_algo)&is.finite(g$stab_boot)
  g$keep_harmonization <- is.finite(g$A_harm) & g$A_harm>=o$harmonization_min
  g$keep_stability <- g$measured & g$stab_algo>stab_algo_threshold & g$rsi>rsi_threshold
  g$keep_stability[is.na(g$keep_stability)] <- FALSE
  g$exclude_saturation <- g$sat>sat_threshold
  g$keep <- g$keep_stability & !g$exclude_saturation & g$keep_harmonization
  g$measurement_status <- ifelse(ids%in%names(measurements),'measured',
    ifelse(plan$table$status=='extrapolated_trivial_unmeasured','quadrant_unmeasured','outside_measurement_policy'))
  g$measurement_status[g$measurement_status=='measured' & !g$measured] <- 'insufficient_repetitions'
  g$quadrant_witness <- plan$table$witness;g$partition_id <- final;g$C_source <- g$C;g$C_reused <- FALSE;g$k_hat <- g$REO
  g$fail_reason <- ifelse(g$keep,'ok',ifelse(!g$measured,'not_measured',ifelse(g$exclude_saturation,'saturation',
    ifelse(!g$keep_harmonization,'harmonization','stability_or_rsi'))))
  m <- borders()
  for(side in names(counts)[m$needs_extension & counts>=o$max_extensions & !blocked])
    restrictions[[length(restrictions)+1L]] <- list(side=side,reason='extension_budget',round=round)
  meta <- list(optimization='cem',penalty=penalty,grid_strategy='adaptive',protocol='staged_cem_v1',
    stability_approach='exhaustive',partition_metric=partition_metric,rsi_policy=rsi_policy,
    rsi_saturation_max=rsi_saturation_max,stab_algo_threshold=stab_algo_threshold,rsi_threshold=rsi_threshold,
    sat_threshold=sat_threshold,k0=k0,B=B,algo_seed=algo_seed,control=control,signature=signature,
    n_cores=n_cores,pilot_elapsed=pilot_elapsed,fitting_elapsed=job_elapsed+harmony_elapsed,
    measurement_elapsed=unname(proc.time()[3]-diagnostic_started),diagnostics_elapsed=agreement_elapsed,checkpoint_file=cache,
    data_fingerprint=digest::digest(unname(X)),comparisons_per_cell=c(original=k0*(k0-1)/2,cross=k0*B),
    familyType='gaussVector',Cquali=Inf,pipeline_version=as.character(utils::packageVersion('CEClust')))
  out <- CECgrid_compatible_output(X,g,pool,final,winners,list(),boots,meta)
  out$measurements <- measurements
  out$executions <- do.call(rbind,executions)
  out$adaptive <- list(initial_C=initial$C,initial_lambda=initial$lambda,control=o,
    mesh=state$mesh,history=do.call(rbind,history),margins=m,restrictions=restrictions,extensions=counts,
    stop_reason=stop_reason,refinement=refinelog,harmonization=harmonies,pilot=pilot,pilot_fits=pilots,
    final_refinement_status=final_pick$status,final_refinement_pending=refinement_pending,
    screening=list(coverage=plan$table,quadrants=plan$quadrants,rejected=plan$rejected,audit_state=audit),
    completed=all(m$observed) && stop_reason%in%c('no_frontiers','adaptive_stable') &&
      final_pick$status%in%c('no_frontiers','adaptive_stable'))
  out$meta$elapsed <- unname(proc.time()[3]-started)
  class(out) <- c('CEC_adaptive_grid',class(out))
  if(isTRUE(save_results)) {
    dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    out$result_file <- file.path(output_dir,paste0('CEC_adaptive_',substr(signature,1,16),'.rds'))
    saveRDS(out,out$result_file)
  }
  out
}
