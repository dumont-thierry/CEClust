CECgrid_values <- function(x,name) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x<=0) || anyDuplicated(x))
    stop(name,' must contain distinct, finite positive values.',call.=FALSE)
  sort(x)
}

CECgrid_bank_selection <- function(parts,C_grid,lambda_grid,d,penalty,stats_name='stats') {
  profile <- CECsolver_environment(penalty)$ca_profile
  ss <- lapply(parts,`[[`,stats_name)
  h <- vapply(ss,function(x)-sum(x$mass*log(x$mass)),0)
  result <- integer(length(C_grid)*length(lambda_grid))
  for(i in seq_along(C_grid)) {
    intercept <- vapply(ss,function(s)profile(s,C_grid[i],0,d)$J,0)
    scores <- intercept + outer(h,lambda_grid)
    ix <- (i-1L)*length(lambda_grid)+seq_along(lambda_grid)
    result[ix] <- max.col(-t(scores),ties.method='first')
  }
  result
}

CECgrid_bank_selection_cells <- function(parts,grid,d,penalty,stats_name='stats') {
  out <- integer(nrow(grid))
  for(C in unique(grid$C)) {
    at <- which(grid$C==C)
    out[at] <- CECgrid_bank_selection(parts,C,grid$lambda[at],d,penalty,stats_name)
  }
  out
}

CECgrid_job <- function(job,X,boots,C_grid,lambda_grid,optimization,penalty,control,
                        Nshots,Nloop,r0,seed,cells=NULL) {
  k <- job$index
  Xfit <- if(job$kind=='original')X else X[boots[[job$rep]],,drop=FALSE]
  control$seed <- ca_seed(seed,optimization,job$kind,job$rep)
  if(optimization=='em_bank') return(CECfit_em_bank(Xfit,control,Xref=X,
    checkpoint_dir=job$starts_checkpoint,resume=isTRUE(job$resume)))
  grid <- if(is.null(cells))expand.grid(lambda=lambda_grid,C=C_grid) else cells
  candidates <- vector('list',nrow(grid));fits <- vector('list',nrow(grid))
  for(i in seq_len(nrow(grid))) {
    s <- ca_seed(control$seed,grid$C[i],grid$lambda[i])
    cc <- control;cc$seed <- s;cc$Cmax <- max(cc$Cmax,grid$C[i])
    f <- CECclassif(Xfit,C=grid$C[i],lambda=grid$lambda[i],r0=r0,Nshots=Nshots,Nloop=Nloop,
      familyType='gaussVector',optimization='cem',penalty=penalty,optimizer_control=cc)
    if(is.null(f)) stop('A CEM repetition failed; no stability is assigned.',call.=FALSE)
    labels <- if(job$kind=='original')f$clusters else
      CECwith_seed(s,predict_frozen(X,f$params,'CEC',s)$labels)
    candidates[[i]] <- list(local_labels=ca_labels(f$clusters),labels=ca_labels(labels),
      local_stats=ca_stats(Xfit,f$clusters),stats=ca_stats(X,labels),params=f$params,
      stop=if(is.null(f$solver))'CEM_completed' else f$solver$status)
    f$phi <- NULL;f$params$phi <- NULL;fits[[i]] <- f
  }
  list(candidates=candidates,cell_winners=seq_len(nrow(grid)),fits=fits,control=control,evaluated_cells=grid)
}

CECgrid_worker <- function(job,arguments) {
  if(!is.null(job$checkpoint) && isTRUE(job$resume) && file.exists(job$checkpoint)) {
    saved <- readRDS(job$checkpoint)
    if(!identical(saved$signature,job$signature))stop('Repetition checkpoint signature mismatch.')
    return(saved$bank)
  }
  bank <- do.call(CECgrid_job,c(list(job=job),arguments))
  if(!is.null(job$checkpoint))CECwrite_rds_atomic(list(signature=job$signature,bank=bank),job$checkpoint)
  bank
}

CECgrid_collect <- function(banks) {
  parts <- list();keys <- new.env(parent=emptyenv());maps <- vector('list',length(banks))
  for(k in seq_along(banks)) {
    maps[[k]] <- integer(length(banks[[k]]$candidates))
    for(j in seq_along(banks[[k]]$candidates)) {
      p <- banks[[k]]$candidates[[j]];key <- digest::digest(p$labels)
      if(exists(key,keys,inherits=FALSE)) {
        id <- get(key,keys)
        if(!identical(parts[[id]]$labels,p$labels)) stop('Partition hash collision.')
      } else {
        id <- length(parts)+1L;parts[[id]] <- list(labels=p$labels,stats=p$stats,origin_bank=k,origin_candidate=j)
        assign(key,id,keys)
      }
      maps[[k]][j] <- id
    }
  }
  list(parts=parts,maps=maps)
}

# Public help: man/CECstability.Rd (maintained by hand).
CECstability <- function(original,bootstrap,partition_metric=c('match','ARI')) {
  partition_metric <- match.arg(partition_metric)
  pp <- c(original,bootstrap)
  if(!is.list(original) || !is.list(bootstrap) || !length(pp) ||
     any(lengths(pp)==0L) || length(unique(lengths(pp)))!=1L || any(vapply(pp,anyNA,TRUE)))
    stop('Partitions must be nonempty lists of complete label vectors of equal length.',call.=FALSE)
  L <- length(original);B <- length(bootstrap)
  canonical <- lapply(pp,ca_labels)
  key <- vapply(canonical,digest::digest,'',algo='xxhash64')
  ids <- match(key,unique(key));representatives <- match(unique(key),key)
  for(i in seq_along(ids))if(!identical(canonical[[i]],canonical[[representatives[ids[i]]]]))
    stop('Partition hash collision in stability calculation.')
  cache <- new.env(parent=emptyenv())
  similarity <- function(i,j) {
    if(ids[i]==ids[j])return(1)
    pair <- paste(sort(c(ids[i],ids[j])),collapse=':')
    if(!exists(pair,cache,inherits=FALSE))assign(pair,
      CECcanonical_similarity(canonical[[i]],canonical[[j]],method=partition_metric),cache)
    get(pair,cache,inherits=FALSE)
  }
  a <- if(L>=2L)utils::combn(seq_len(L),2,function(ij)similarity(ij[1],ij[2])) else numeric()
  b <- if(L && B)unlist(lapply(seq_len(L),function(i)vapply(L+seq_len(B),function(j)similarity(i,j),0)),use.names=FALSE) else numeric()
  SA <- if(length(a))mean(a) else NA_real_;SB <- if(length(b))mean(b) else NA_real_
  list(Stab_algo=SA,Stab_boot=SB,RSI=if(is.finite(SA)&&SA>0)SB/SA else NA_real_,
    original_pairs=length(a),cross_pairs=length(b),definition='exhaustive_original_bootstrap',
    algorithm_agreements=a,bootstrap_agreements=b)
}

CECgrid_diagnostics <- function(winners,parts,eligible,L,B,partition_metric='match',cache=NULL) {
  out <- data.frame(Stab_algo=rep(NA_real_,nrow(winners)),Stab_boot=NA_real_,RSI=NA_real_)
  if(is.null(cache))cache <- new.env(parent=emptyenv())
  if(is.null(cache$patterns))cache$patterns <- new.env(parent=emptyenv())
  if(is.null(cache$pairs))cache$pairs <- new.env(parent=emptyenv())
  patterns <- cache$patterns;pairs <- cache$pairs
  ai <- if(L>=2L)t(utils::combn(seq_len(L),2)) else matrix(integer(),0,2)
  bi <- if(L && B)as.matrix(expand.grid(i=seq_len(L),j=L+seq_len(B))) else matrix(integer(),0,2)
  sim <- function(i,j) {
    if(i==j)return(1)
    key <- paste(sort(c(i,j)),collapse=':')
    if(!exists(key,pairs,inherits=FALSE))assign(key,
      CECcanonical_similarity(parts[[i]]$labels,parts[[j]]$labels,method=partition_metric),pairs)
    get(key,pairs)
  }
  for(i in which(eligible)) {
    ids <- winners[i,];key <- paste(ids,collapse=',')
    if(!exists(key,patterns,inherits=FALSE)) {
      SA <- if(nrow(ai))mean(vapply(seq_len(nrow(ai)),function(j)sim(ids[ai[j,1]],ids[ai[j,2]]),0)) else NA_real_
      SB <- if(nrow(bi))mean(vapply(seq_len(nrow(bi)),function(j)sim(ids[bi[j,1]],ids[bi[j,2]]),0)) else NA_real_
      assign(key,c(SA,SB,if(is.finite(SA)&&SA>0)SB/SA else NA_real_),patterns)
    }
    out[i,] <- get(key,patterns)
  }
  attr(out,'comparisons_per_cell') <- c(original=nrow(ai),cross=nrow(bi))
  attr(out,'unique_patterns') <- length(ls(patterns));out
}

CECgrid_compatible_output <- function(X,grid,parts,final,winners,banks,boots,meta) {
  C_grid <- sort(unique(grid$C));runs <- vector('list',length(C_grid))
  profile <- CECsolver_environment(meta$penalty)$ca_profile
  templates <- new.env(parent=emptyenv())
  for(ci in seq_along(C_grid)) {
    ix <- which(grid$C==C_grid[ci]);ix <- ix[order(grid$lambda[ix])]
    best <- vector('list',length(ix));ds <- grid[ix,,drop=FALSE]
    for(j in seq_along(ix)) {
      i <- ix[j];p <- parts[[final[i]]];key <- paste(ci,final[i],sep=':')
      if(!exists(key,templates,inherits=FALSE)) {
        f <- CECfit_from_partition(X,p$labels,grid$C[i],1,meta$penalty,meta$optimization)
        f$phi <- NULL;f$params$phi <- NULL;assign(key,f,templates)
      }
      f <- get(key,templates);pr <- profile(p$stats,grid$C[i],grid$lambda[i],ncol(X))
      f$params$lambda <- grid$lambda[i];f$Hphi <- pr$J/grid$lambda[i]
      f$objective <- pr$J;f$entropy <- pr$J_raw;f$penalty_value <- pr$BIC_penalty_relative
      projection <- list(clusters=p$labels,H=f$Hphi,H_class=pr$h,
        H_cond=(pr$J_raw-grid$lambda[i]*pr$h)/grid$lambda[i],REO=pr$K)
      best[[j]] <- list(lambda=grid$lambda[i],partition=p$labels,fit=f,projection=projection,
        Hphi=f$Hphi,H_projected=f$Hphi,H_raw=f$Hphi,H_class_projected=pr$h,
        H_cond_projected=projection$H_cond,REO=pr$K,k_hat=pr$K,source='candidate_pool',
        source_index=final[i],criterion='projected_H',optimization=meta$optimization,penalty=meta$penalty)
    }
    bs <- ds;bs$H_projected <- bs$Hphi;bs$H_raw <- bs$Hphi;bs$k_hat <- bs$REO
    bp <- list(best=best,summary=bs,source='all',criterion='projected_H',
      optimization=meta$optimization,penalty=meta$penalty,bank_minimum=meta$optimization=='em_bank')
    diag_summary <- data.frame(lambda=ds$lambda,stab0_mean=ds$stab_algo,stabB_mean=ds$stab_boot,
      stab_ratio_mean=ds$rsi,sat=ds$sat,measured=ds$measured,
      n_sat_clusters=ds$n_sat_clusters,any_bound_saturation=ds$any_bound_saturation)
    ld <- list(summary=diag_summary,lambda_grid=ds$lambda,partition_metric=meta$partition_metric,
      stability_approach='exhaustive',linked=FALSE,details=NULL,best_parts=bp,
      optimization=meta$optimization,penalty=meta$penalty)
    class(ld) <- c('CEC_solver_lambda','list')
    keep <- ds;keep$rsi <- ds$rsi;keep$stab_algo <- ds$stab_algo
    runs[[ci]] <- list(Z=X,meta=c(meta,list(C=C_grid[ci],C_source=C_grid[ci],C_reused=FALSE)),
      best_parts=bp,lambda_diag=ld,lambda_keep=list(summary=keep,stable_lambda=ds$lambda[ds$keep]),compact=TRUE)
  }
  out <- list(Z=as.data.frame(X),data=list(Z=as.data.frame(X)),display_data=as.data.frame(X),
    C_grid=C_grid,lambda_grid=sort(unique(grid$lambda)),summary=grid,runs_by_C=runs,
    meta=meta,banks=banks,parts=parts,winners=winners,selected=final,bootstrap_indices=boots,
    compact=TRUE,optimization=meta$optimization,penalty=meta$penalty)
  class(out) <- c('CEC_solver_grid','CEC_bound_grid','CEC_grid','list');out
}

CECfit_solver_grid <- function(Z,C_grid,lambda_grid,k0=20L,B=20L,r0=20L,
  Nshots_fresh=5L,Nloop=1000L,algo_seed=1L,n_cores=1L,
  optimization='em_bank',penalty='none',optimizer_control=list(),
  rsi_policy=c('all','unsaturated','nontrivial_unsaturated','saturation_band'),
  rsi_saturation_max=.10,stab_algo_threshold=.85,rsi_threshold=.95,sat_threshold=0,
  partition_metric=c('match','ARI'),bootstrap_indices=NULL,checkpoint_dir=FALSE,
  resume=FALSE,force_recompute=FALSE,save_results=FALSE,output_dir=tempdir(),verbose=FALSE,
  banks=NULL,cells=NULL,diagnostics=TRUE,reference_bootstrap=TRUE) {
  X <- CECnumeric_data(Z);C_grid <- CECgrid_values(C_grid,'C_grid')
  lambda_grid <- CECgrid_values(lambda_grid,'lambda_grid')
  if(!is.null(cells)) {
    if(!is.data.frame(cells) || !all(c('C','lambda')%in%names(cells)) || !nrow(cells) ||
       any(!is.finite(as.matrix(cells[,c('C','lambda')])) | as.matrix(cells[,c('C','lambda')])<=0) ||
       anyDuplicated(cells[,c('C','lambda')]))stop('Invalid evaluated cells.')
    cells <- cells[,c('C','lambda'),drop=FALSE]
    C_grid <- sort(unique(cells$C));lambda_grid <- sort(unique(cells$lambda))
  }
  k0 <- CECinteger(k0,'k0');B <- CECinteger(B,'B',0L);n_cores <- CECinteger(n_cores,'n_cores')
  algo_seed <- CECinteger(algo_seed,'algo_seed',0L)
  rsi_policy <- match.arg(rsi_policy);partition_metric <- match.arg(partition_metric)
  optimization <- match.arg(optimization,c('cem','em_bank'));penalty <- match.arg(penalty,c('none','bic'))
  if(!is.logical(reference_bootstrap)||length(reference_bootstrap)!=1L||is.na(reference_bootstrap))
    stop('reference_bootstrap must be TRUE or FALSE.')
  for(x in list(rsi_saturation_max,stab_algo_threshold,rsi_threshold,sat_threshold))
    if(!is.numeric(x)||length(x)!=1L||!is.finite(x)||x<0)stop('Diagnostic thresholds must be finite nonnegative numbers.')
  if(rsi_saturation_max>1 || sat_threshold>1 || stab_algo_threshold>1)stop('Saturation and agreement thresholds must not exceed one.')
  control <- CECoptimization_control(optimizer_control,nrow(X),max(C_grid),Nshots_fresh,Nloop,r0)
  if(!is.null(control$seed))algo_seed <- control$seed
  if(is.null(bootstrap_indices)) boots <- CECwith_seed(ca_seed(algo_seed,'bootstrap_indices'),
    replicate(B,sample.int(nrow(X),nrow(X),replace=TRUE),simplify=FALSE)) else boots <- bootstrap_indices
  if(!is.list(boots)||length(boots)!=B||any(vapply(boots,function(x)length(x)!=nrow(X)||anyNA(x)||
    any(x!=floor(x))||any(x<1|x>nrow(X)),TRUE)))stop('Invalid bootstrap_indices.',call.=FALSE)
  jobs <- lapply(seq_len(k0+B),function(k)list(index=k,kind=if(k<=k0)'original' else 'bootstrap',rep=if(k<=k0)k else k-k0))
  signature <- digest::digest(list(version=as.character(utils::packageVersion('CEClust')),X=X,boots=boots,k0=k0,B=B,
    optimization=optimization,control=control,seed=algo_seed,
    grid=if(optimization=='cem')list(C_grid,lambda_grid,penalty,cells) else NULL))
  cache <- if(isFALSE(checkpoint_dir)||is.null(checkpoint_dir))NULL else {
    if(!is.character(checkpoint_dir)||length(checkpoint_dir)!=1L)stop('checkpoint_dir must be FALSE or a directory path.')
    dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE);file.path(checkpoint_dir,paste0('banks_',signature,'.rds'))
  }
  if(is.null(banks)&&!is.null(cache)&&isTRUE(resume)&&!isTRUE(force_recompute)&&file.exists(cache)) {
    saved <- readRDS(cache)
    if(!identical(saved$signature,signature))stop('Checkpoint signature mismatch.')
    banks <- saved$banks
  }
  start <- proc.time()[3]
  if(is.null(banks)) {
    if(!is.null(cache)) {
      folder <- file.path(dirname(cache),paste0('partial_',signature))
      dir.create(folder,recursive=TRUE,showWarnings=FALSE)
      jobs <- lapply(jobs,function(job) {
        job$signature <- signature;job$resume <- isTRUE(resume)&&!isTRUE(force_recompute)
        name <- sprintf('%s_%04d',job$kind,job$rep)
        job$checkpoint <- file.path(folder,paste0(name,'.rds'))
        job$starts_checkpoint <- if(optimization=='em_bank')file.path(folder,name) else NULL
        job
      })
    }
    args <- list(X=X,boots=boots,C_grid=C_grid,lambda_grid=lambda_grid,optimization=optimization,
      penalty=penalty,control=control,Nshots=Nshots_fresh,Nloop=Nloop,r0=r0,seed=algo_seed,cells=cells)
    if(n_cores>1L) {
      old_r_libs <- Sys.getenv('R_LIBS',unset=NA_character_)
      Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
      cl <- tryCatch(parallel::makePSOCKcluster(min(n_cores,length(jobs))),finally={
        if(is.na(old_r_libs))Sys.unsetenv('R_LIBS') else Sys.setenv(R_LIBS=old_r_libs)
      })
      on.exit(parallel::stopCluster(cl),add=TRUE)
      # Do not serialize a package-namespace closure before setting worker
      # libraries: it can load an older installed CEClust during deserialization.
      setup_worker <- function(paths,expected) {
        .libPaths(paths)
        ns <- loadNamespace('CEClust')
        if(as.character(getNamespaceVersion(ns))!=expected)
          stop('Worker loaded a different CEClust version.')
        NULL
      }
      environment(setup_worker) <- baseenv()
      attributes(setup_worker) <- NULL
      parallel::clusterCall(cl,setup_worker,.libPaths(),as.character(utils::packageVersion('CEClust')))
      banks <- parallel::parLapplyLB(cl,jobs,CECgrid_worker,arguments=args)
    } else banks <- lapply(jobs,CECgrid_worker,arguments=args)
    if(!is.null(cache)) {
      tmp <- tempfile('banks_',dirname(cache));saveRDS(list(signature=signature,banks=banks),tmp)
      if(!file.copy(tmp,cache,overwrite=TRUE))stop('Cannot save checkpoint.')
      unlink(tmp)
    }
  }
  if(length(banks)!=k0+B)stop('Candidate bank count does not match repetitions.')
  fitting_elapsed <- unname(proc.time()[3]-start)
  pool <- CECgrid_collect(banks)
  grid <- if(is.null(cells))expand.grid(lambda=lambda_grid,C=C_grid) else cells;M <- nrow(grid)
  winners <- matrix(NA_integer_,M,k0+B)
  for(k in seq_along(banks)) {
    ids <- if(optimization=='em_bank')CECgrid_bank_selection_cells(banks[[k]]$candidates,grid,ncol(X),penalty,'local_stats') else banks[[k]]$cell_winners
    winners[,k] <- pool$maps[[k]][ids]
  }
  if(optimization=='em_bank')final <- CECgrid_bank_selection_cells(pool$parts,grid,ncol(X),penalty) else {
    final <- vapply(seq_len(M),function(i){ids<-unique(winners[i,]);
      ids[CECscore_candidates(pool$parts[ids],grid$C[i],grid$lambda[i],ncol(X),penalty)$best]},0L)
  }
  profile <- CECsolver_environment(penalty)$ca_profile
  profiles <- lapply(seq_len(M),function(i)profile(pool$parts[[final[i]]]$stats,grid$C[i],grid$lambda[i],ncol(X)))
  grid$REO <- vapply(profiles,`[[`,0L,'K');grid$sat <- vapply(profiles,`[[`,0,'saturation')
  eligible <- switch(rsi_policy,all=rep(TRUE,M),unsaturated=grid$sat==0,
    nontrivial_unsaturated=grid$REO>1 & grid$sat==0,saturation_band=grid$sat<=rsi_saturation_max)
  if(!isTRUE(diagnostics))eligible[] <- FALSE
  diagnostic <- CECgrid_diagnostics(winners,pool$parts,eligible,k0,B,partition_metric)
  grid$id <- seq_len(M)
  reference_cfg <- list(X=X,boots=boots,seed=algo_seed,settings=CECadaptive_cem_settings(control),
    penalty=penalty,signature=signature,cache=NULL,resume=isTRUE(resume)&&!force_recompute)
  if(!is.null(cache)) {
    reference_cfg$cache <- file.path(dirname(cache),paste0('reference_',signature))
    dir.create(reference_cfg$cache,recursive=TRUE,showWarnings=FALSE)
  }
  measurements <- list()
  for(i in which(eligible)) {
    d <- as.list(diagnostic[i,]);d$definition <- 'exhaustive_original_bootstrap'
    d$original_pairs <- k0*(k0-1)/2;d$cross_pairs <- k0*B
    z <- list(cell_id=i,C=grid$C[i],lambda=grid$lambda[i],diagnostics=d,
      original_partition_ids=winners[i,seq_len(k0)],
      original_parts=lapply(winners[i,seq_len(k0)],function(j)pool$parts[[j]]))
    runner <- function(jobs,worker)CECreference_parallel(jobs,reference_cfg,n_cores)
    z <- CECreference_stability(z,pool$parts[[final[i]]],reference_cfg,grid[i,],partition_metric,
      stab_algo_threshold,reference_bootstrap,runner,rsi_min=rsi_threshold)
    measurements[[as.character(i)]] <- z
    diagnostic[i,] <- c(z$diagnostics$Stab_algo,z$diagnostics$Stab_boot,z$diagnostics$RSI)
  }
  grid$objective <- vapply(profiles,`[[`,0,'J');grid$Hphi <- grid$objective/grid$lambda
  grid$entropy <- vapply(profiles,`[[`,0,'J_raw');grid$penalty_value <- vapply(profiles,`[[`,0,'BIC_penalty_relative')
  grid$stab_algo <- diagnostic$Stab_algo;grid$stab_boot <- diagnostic$Stab_boot;grid$rsi <- diagnostic$RSI
  grid$measured <- eligible & is.finite(grid$stab_algo) & is.finite(grid$stab_boot)
  grid <- CECattach_reference_summary(grid,measurements)
  grid$keep_stability <- grid$measured & grid$stab_algo>stab_algo_threshold & grid$rsi>rsi_threshold
  grid$keep_stability[is.na(grid$keep_stability)] <- FALSE
  grid$exclude_saturation <- grid$sat>sat_threshold;grid$keep <- grid$keep_stability & !grid$exclude_saturation
  grid$n_sat_clusters <- vapply(profiles,function(p)sum(p$saturated),0L)
  grid$any_bound_saturation <- grid$n_sat_clusters>0
  grid$partition_id <- final;grid$C_source <- grid$C;grid$C_reused <- FALSE;grid$k_hat <- grid$REO
  grid$fail_reason <- ifelse(grid$keep,'ok',ifelse(!grid$measured,'not_measured',ifelse(grid$exclude_saturation,'saturation','stability_or_rsi')))
  meta <- list(optimization=optimization,penalty=penalty,stability_approach='exhaustive',partition_metric=partition_metric,
    rsi_policy=rsi_policy,rsi_saturation_max=rsi_saturation_max,stab_algo_threshold=stab_algo_threshold,
    rsi_threshold=rsi_threshold,sat_threshold=sat_threshold,k0=k0,B=B,algo_seed=algo_seed,control=control,
    signature=digest::digest(list(signature,reference_bootstrap,rsi_policy,rsi_saturation_max,
      stab_algo_threshold,rsi_threshold,sat_threshold,partition_metric)),
    n_cores=n_cores,fitting_elapsed=fitting_elapsed,checkpoint_file=cache,
    data_fingerprint=digest::digest(unname(X)),
    comparisons_per_cell=attr(diagnostic,'comparisons_per_cell'),familyType='gaussVector',Cquali=Inf,
    reference_bootstrap=reference_bootstrap,pipeline_version=as.character(utils::packageVersion('CEClust')))
  out <- CECgrid_compatible_output(X,grid,pool$parts,final,winners,banks,boots,meta)
  out$measurements <- measurements
  if(isTRUE(save_results)) {
    dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    file <- file.path(output_dir,paste0('CEC_',optimization,'_',penalty,'_',substr(meta$signature,1,12),'.rds'))
    saveRDS(out,file);out$result_file <- file
  }
  if(isTRUE(verbose))message(M,' cells; ',length(banks),' independent repetitions; ',sum(grid$keep),' retained.')
  out
}
