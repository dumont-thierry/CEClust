# Adaptive window and rectangular leaves. No research files are sourced at runtime.
ceca_current <- function(state)unname(state$parts[as.character(state$mesh$leaves$id)])
ceca_diagram <- function(state) {
  g <- state$mesh$leaves;p <- ceca_current(state)
  g$REO <- vapply(p,`[[`,0L,'K');g$saturation <- vapply(p,`[[`,0,'saturation');g
}

CECadaptive_control <- function(control,sat_threshold) {
  defaults <- list(max_cells=12000L,max_extensions=16L,max_rounds=12L,max_depth=8L,
    patience=2L,position_tolerance=.125,distinguish_partitions=TRUE,
    frontier_saturation=sat_threshold,harmonize=TRUE,harmonization_min=.95,
    trivial_quadrants=FALSE,audit_fraction=.10,anchor_stability=.95)
  if(!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop('adaptive_control must be a named list without duplicates.',call.=FALSE)
  unknown <- setdiff(names(control),names(defaults))
  if(length(unknown))stop('Unknown adaptive_control: ',paste(unknown,collapse=', '),call.=FALSE)
  out <- utils::modifyList(defaults,control)
  for(nm in c('max_cells','max_extensions','max_rounds','max_depth','patience'))
    out[[nm]] <- CECinteger(out[[nm]],nm,if(nm=='max_extensions')0L else if(nm=='max_rounds')2L else 1L)
  for(nm in c('position_tolerance','frontier_saturation','harmonization_min','audit_fraction','anchor_stability'))
    if(!is.numeric(out[[nm]]) || length(out[[nm]])!=1L || !is.finite(out[[nm]]) || out[[nm]]<0)
      stop(nm,' must be finite and nonnegative.',call.=FALSE)
  if(out$position_tolerance<=0 || out$audit_fraction<=0 || out$anchor_stability<=0 ||
     any(unlist(out[c('frontier_saturation','harmonization_min','audit_fraction','anchor_stability')])>1) ||
     out$frontier_saturation<sat_threshold)stop('Invalid adaptive thresholds.',call.=FALSE)
  for(nm in c('harmonize','trivial_quadrants','distinguish_partitions'))
    if(!is.logical(out[[nm]]) || length(out[[nm]])!=1L || is.na(out[[nm]]))stop(nm,' must be TRUE or FALSE.')
  out
}

CECadaptive_initial <- function(X,C_grid,lambda_grid) {
  if(is.null(C_grid)) {
    s <- ca_stats(X,rep(1L,nrow(X)))
    peak <- exp(-.5*(ncol(X)*log(2*pi)+s$logdet))
    if(!is.finite(peak) || peak<=0)stop('Supply C_grid for singular data or extreme density scales.',call.=FALSE)
    C_grid <- peak*exp(seq(log(.25),log(4),length.out=5))
  }
  if(is.null(lambda_grid))lambda_grid <- seq(.5,2,by=.375)
  C_grid <- CECgrid_values(C_grid,'C_grid');lambda_grid <- CECgrid_values(lambda_grid,'lambda_grid')
  regular <- function(z)length(z)>=2 && max(abs(diff(z)-diff(z)[1]))<1e-8*max(1,abs(diff(z)[1]))
  if(!regular(log(C_grid)) || !regular(lambda_grid))
    stop('Adaptive initial C_grid must be equally spaced on the log scale, and lambda_grid on the linear scale (at least two each).',call.=FALSE)
  if(lambda_grid[1]-diff(lambda_grid)[1]/2<=0)
    stop('Initial lambda cell boundaries must be strictly positive.',call.=FALSE)
  list(C=C_grid,lambda=lambda_grid)
}

CECadaptive_scorer <- function(pool,banks,maps,X,control,penalty) {
  profile <- CECsolver_environment(penalty)$ca_profile
  cache <- new.env(parent=emptyenv());d <- ncol(X)
  entropy <- function(s)-sum(s$mass*log(s$mass))
  h <- vapply(pool,function(p)entropy(p$stats),0)
  bh <- lapply(banks,function(b)vapply(b$candidates,function(p)entropy(p$local_stats),0))
  function(C,lambda) {
    if(C>control$Cmax*(1+1e-12))stop('An adaptive cell exceeds the bank Cmax.')
    key <- sprintf('%.17g',C)
    if(!exists(key,cache,inherits=FALSE)) {
      p <- lapply(pool,function(p)profile(p$stats,C,0,d))
      intercept <- lapply(banks,function(b)vapply(b$candidates,function(p)profile(p$local_stats,C,0,d)$J,0))
      assign(key,list(a=vapply(p,`[[`,0,'J'),p=p,intercept=intercept),cache)
    }
    s <- get(key,cache);id <- which.min(s$a+lambda*h)
    winners <- vapply(seq_along(banks),function(k)maps[[k]][which.min(s$intercept[[k]]+lambda*bh[[k]])],0L)
    p <- s$p[[id]];p$J <- s$a[id]+lambda*h[id]
    c(p,list(labels=pool[[id]]$labels,stats=pool[[id]]$stats,id=id,winners=winners))
  }
}

CECfit_adaptive_grid <- function(Z,C_grid=NULL,lambda_grid=NULL,k0=20L,B=20L,r0=20L,
  Nshots_fresh=5L,Nloop=1000L,algo_seed=1L,n_cores=1L,optimization='em_bank',penalty='none',
  optimizer_control=list(),rsi_policy=c('all','unsaturated','nontrivial_unsaturated','saturation_band'),
  rsi_saturation_max=.10,stab_algo_threshold=.85,rsi_threshold=.95,sat_threshold=0,
  partition_metric=c('match','ARI'),bootstrap_indices=NULL,checkpoint_dir=FALSE,
  resume=FALSE,force_recompute=FALSE,save_results=FALSE,output_dir=tempdir(),verbose=FALSE,
  banks=NULL,adaptive_control=list()) {
  started <- proc.time()[3];X <- CECnumeric_data(Z)
  optimization <- match.arg(optimization,c('cem','em_bank'));penalty <- match.arg(penalty,c('none','bic'))
  rsi_policy <- match.arg(rsi_policy);partition_metric <- match.arg(partition_metric)
  initial <- CECadaptive_initial(X,C_grid,lambda_grid);o <- CECadaptive_control(adaptive_control,sat_threshold)
  if(length(initial$C)*length(initial$lambda)>o$max_cells)stop('Initial grid exceeds max_cells.')
  if(o$trivial_quadrants && rsi_policy=='nontrivial_unsaturated')
    stop('Trivial quadrants require a policy which measures trivial anchors.',call.=FALSE)
  if(optimization=='cem' && !is.null(optimizer_control$orders))stop('Explicit orders are for em_bank only; CEM uses Rmax and subdivision controls.')
  if(is.null(optimizer_control$Cmax))optimizer_control$Cmax <- if(optimization=='em_bank')64*max(initial$C) else max(initial$C)
  mesh <- ceca_mesh(initial$C,initial$lambda)
  args <- list(Z=X,C_grid=initial$C,lambda_grid=initial$lambda,k0=k0,B=B,r0=r0,
    Nshots_fresh=Nshots_fresh,Nloop=Nloop,algo_seed=algo_seed,n_cores=n_cores,
    optimization=optimization,penalty=penalty,optimizer_control=optimizer_control,rsi_policy=rsi_policy,
    rsi_saturation_max=rsi_saturation_max,stab_algo_threshold=stab_algo_threshold,rsi_threshold=rsi_threshold,
    sat_threshold=sat_threshold,partition_metric=partition_metric,bootstrap_indices=bootstrap_indices,
    checkpoint_dir=checkpoint_dir,resume=resume,force_recompute=force_recompute,
    banks=banks,diagnostics=FALSE,cells=mesh$leaves[,c('C','lambda')])
  base <- do.call(CECfit_solver_grid,args)
  control <- base$meta$control;algo_seed <- base$meta$algo_seed;k0 <- base$meta$k0;B <- base$meta$B
  args$bootstrap_indices <- base$bootstrap_indices;args$banks <- NULL
  settings <- list(sat_max=sat_threshold,sat_max2=o$frontier_saturation,
    max_cells=o$max_cells,refinement=o[c('max_rounds','max_depth','patience','position_tolerance','distinguish_partitions')])
  pool <- base$parts;keys <- new.env(parent=emptyenv())
  for(i in seq_along(pool))assign(digest::digest(pool[[i]]$labels),i,keys)
  register <- function(p) {
    key <- digest::digest(p$labels)
    if(exists(key,keys,inherits=FALSE)) {
      id <- get(key,keys);if(!identical(pool[[id]]$labels,p$labels))stop('Partition hash collision.')
    } else {
      id <- length(pool)+1L;pool[[id]] <<- list(labels=p$labels,stats=p$stats);assign(key,id,keys)
    }
    id
  }
  state <- list(mesh=mesh,parts=list());independent <- list();before_harmony <- list()
  history <- list();restrictions <- list();harmonies <- list();refinelog <- list()
  right <- NULL;counts <- setNames(integer(4),c('left','right','lower','upper'));round <- 1L
  blocked <- setNames(rep(FALSE,4),names(counts))
  all_banks <- base$banks;fit_time <- base$meta$fitting_elapsed
  if(optimization=='em_bank') {
    collect <- CECgrid_collect(base$banks)
    score <- CECadaptive_scorer(pool,base$banks,collect$maps,X,control,penalty)
  }
  profile <- CECsolver_environment(penalty)$ca_profile
  accept <- function(f,ids) {
    map <- vapply(f$parts,register,0L)
    for(j in seq_along(ids)) {
      key <- as.character(ids[j]);id <- map[f$selected[j]];p <- pool[[id]]
      g <- state$mesh$leaves[match(ids[j],state$mesh$leaves$id),]
      state$parts[[key]] <<- c(profile(p$stats,g$C,g$lambda,ncol(X)),p,list(id=id))
      before_harmony[[key]] <<- id
      independent[[key]] <<- map[f$winners[j,]]
    }
  }
  harmonize <- function(seed_ids=NULL) {
    if(optimization!='cem' || !o$harmonize)return(invisible(NULL))
    q <- q_settings();q$rinit <- control$Rmax;q$shots <- control$Ninit;q$Nloop <- control$max_iter
    q$tolerance <- control$tolerance
    cfg <- list(X=X,seed=algo_seed,settings=q)
    fun <- ceca_harmonize;environment(fun) <- CECsolver_environment(penalty)
    h <- CECwith_seed(ca_seed(algo_seed,'adaptive_harmony',round),
      fun(state,cfg,seed_ids=seed_ids,tag=as.character(round)))
    h$evaluated_cells <- nrow(state$mesh$leaves)
    h$initial_cells <- seed_ids
    state <<- h$state
    for(id in names(state$parts))state$parts[[id]]$id <<- register(state$parts[[id]])
    harmonies[[length(harmonies)+1L]] <<- h[names(h)!='state']
  }
  evaluate <- function(first=FALSE) {
    g <- state$mesh$leaves;new <- which(!as.character(g$id)%in%names(state$parts))
    if(!length(new))return(invisible(NULL))
    if(optimization=='em_bank')for(i in new) {
      key <- as.character(g$id[i]);p <- score(g$C[i],g$lambda[i]);state$parts[[key]] <<- p
      independent[[key]] <<- p$winners;before_harmony[[key]] <<- p$id
    } else {
      if(first)f <- base else {
        args$cells <- g[new,c('C','lambda')]
        args$optimizer_control$Cmax <- max(control$Cmax,max(args$cells$C))
        f <- do.call(CECfit_solver_grid,args);fit_time <<- fit_time+f$meta$fitting_elapsed
        for(k in seq_along(all_banks)) {
          all_banks[[k]]$candidates <<- c(all_banks[[k]]$candidates,f$banks[[k]]$candidates)
          all_banks[[k]]$fits <<- c(all_banks[[k]]$fits,f$banks[[k]]$fits)
          all_banks[[k]]$evaluated_cells <<- rbind(all_banks[[k]]$evaluated_cells,f$banks[[k]]$evaluated_cells)
          all_banks[[k]]$cell_winners <<- seq_along(all_banks[[k]]$candidates)
        }
      }
      accept(f,g$id[new]);harmonize(if(first)NULL else g$id[new])
    }
  }
  borders <- function(append=FALSE) {
    right <<- ceca_right_refresh(state,right,append)
    ceca_right_status(state,settings,right)
  }
  snapshot <- function(stage,status=borders()) {
    g <- state$mesh$leaves
    history[[length(history)+1L]] <<- data.frame(step=length(history)+1L,stage=stage,round=round,
      cells=nrow(g),Cmin=min(g$C),Cmax=max(g$C),lambda_min=min(g$lambda),lambda_max=max(g$lambda),
      left=unname(status$observed['left']),right=unname(status$observed['right']),
      lower=unname(status$observed['lower']),upper=unname(status$observed['upper']),
      right_reason=status$right_stopping$reason,elapsed=unname(proc.time()[3]-started))
    if(isTRUE(verbose))message(stage,': ',nrow(g),' cells; right: ',status$right_stopping$reason)
  }
  extend <- function() {
    repeat {
      m <- borders();allowed <- m$needs_extension & counts<o$max_extensions & !blocked
      if(!any(allowed))break
      side <- c('upper','left','lower','right')[which(allowed[c('upper','left','lower','right')])[1]]
      plan <- ceca_expand(state$mesh,side)
      reason <- if(any(!is.finite(plan$mesh$leaves$C)) || any(plan$mesh$leaves$C<=0) ||
        any(!is.finite(plan$mesh$leaves$lambda)) || any(plan$mesh$leaves$lambda<=0))'numeric_limit' else
        if(optimization=='em_bank' && any(plan$mesh$leaves$C>control$Cmax*(1+1e-12)))'EM_Cmax' else
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
    before <- state;ids <- state$mesh$leaves$id[pick$indices];plan <- ceca_refine(state$mesh,pick$indices)
    if(optimization=='em_bank' && any(plan$mesh$leaves$C>control$Cmax*(1+1e-12))) {
      stop_reason <- 'EM_Cmax_refinement';break
    }
    state$mesh <- plan$mesh;round <- round+1L;evaluate()
    update <- ceca_refinement_update(before,state,memory,ids,settings,scale)
    memory <- update$memory;refinelog[[length(refinelog)+1L]] <- cbind(round=round,update$diagnostics)
    snapshot('refinement');extend();memory <- ceca_wake_changed(state,memory,o$frontier_saturation)
  }
  g <- ceca_diagram(state);p <- ceca_current(state);ids <- as.character(g$id)
  final <- vapply(p,`[[`,0L,'id');winners <- do.call(rbind,independent[ids])
  profiles <- lapply(seq_len(nrow(g)),function(i)profile(pool[[final[i]]]$stats,g$C[i],g$lambda[i],ncol(X)))
  g$sat <- g$saturation;g$objective <- vapply(profiles,`[[`,0,'J');g$Hphi <- g$objective/g$lambda
  g$entropy <- vapply(profiles,`[[`,0,'J_raw');g$penalty_value <- vapply(profiles,`[[`,0,'BIC_penalty_relative')
  g$n_sat_clusters <- vapply(profiles,function(p)sum(p$saturated),0L);g$any_bound_saturation <- g$n_sat_clusters>0L
  g$A_harm <- vapply(seq_along(ids),function(i)ca_sim(pool[[before_harmony[[ids[i]]]]]$labels,pool[[final[i]]]$labels),0)
  g$keep_harmonization <- g$A_harm>=o$harmonization_min
  eligible <- switch(rsi_policy,all=rep(TRUE,nrow(g)),unsaturated=g$sat==0,
    nontrivial_unsaturated=g$REO>1 & g$sat==0,saturation_band=g$sat<=rsi_saturation_max)
  diagnostic_start <- proc.time()[3]
  diagnostic <- CECadaptive_diagnostics(state,g,winners,pool,eligible,k0,B,partition_metric,o,
    rsi_threshold,stab_algo_threshold,algo_seed,n_cores)
  g$stab_algo <- diagnostic$values$Stab_algo;g$stab_boot <- diagnostic$values$Stab_boot;g$rsi <- diagnostic$values$RSI
  g$measured <- is.finite(g$stab_algo)&is.finite(g$stab_boot)
  g$keep_stability <- g$measured & g$stab_algo>stab_algo_threshold & g$rsi>rsi_threshold
  g$keep_stability[is.na(g$keep_stability)] <- FALSE
  g$exclude_saturation <- g$sat>sat_threshold
  g$keep <- g$keep_stability & !g$exclude_saturation & g$keep_harmonization
  g$measurement_status <- diagnostic$status;g$quadrant_witness <- diagnostic$witness
  g$measurement_status[g$measurement_status=='measured' & !g$measured] <- 'insufficient_repetitions'
  g$partition_id <- final;g$C_source <- g$C;g$C_reused <- FALSE;g$k_hat <- g$REO
  g$fail_reason <- ifelse(g$keep,'ok',ifelse(!g$measured,'not_measured',
    ifelse(g$exclude_saturation,'saturation',ifelse(!g$keep_harmonization,'harmonization','stability_or_rsi'))))
  meta <- base$meta;meta$grid_strategy <- 'adaptive';meta$fitting_elapsed <- fit_time
  meta$diagnostics_elapsed <- unname(proc.time()[3]-diagnostic_start)
  meta$signature <- digest::digest(list(base$meta$signature,initial,o,penalty,rsi_policy,rsi_saturation_max,
    stab_algo_threshold,rsi_threshold,sat_threshold,partition_metric))
  m <- borders()
  for(side in names(counts)[m$needs_extension & counts>=o$max_extensions & !blocked])
    restrictions[[length(restrictions)+1L]] <- list(side=side,reason='extension_budget',round=round)
  out <- CECgrid_compatible_output(X,g,pool,final,winners,all_banks,base$bootstrap_indices,meta)
  out$adaptive <- list(initial_C=initial$C,initial_lambda=initial$lambda,control=o,
    mesh=state$mesh,history=do.call(rbind,history),margins=m,restrictions=restrictions,
    extensions=counts,stop_reason=stop_reason,refinement=refinelog,harmonization=harmonies,
    screening=diagnostic$screening,
    completed=all(m$observed) && stop_reason%in%c('no_frontiers','adaptive_stable'))
  out$meta$elapsed <- unname(proc.time()[3]-started)
  class(out) <- c('CEC_adaptive_grid',class(out))
  if(isTRUE(save_results)) {
    dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    out$result_file <- file.path(output_dir,paste0('CEC_adaptive_',substr(meta$signature,1,16),'.rds'))
    saveRDS(out,out$result_file)
  }
  out
}

CECadaptive_diagnostics <- function(state,g,winners,pool,eligible,L,B,metric,o,rsi_min,algo_min,seed,workers) {
  M <- nrow(g);status <- ifelse(eligible,'measured','outside_measurement_policy')
  witness <- rep(NA_integer_,M)
  if(!o$trivial_quadrants) return(list(values=CECgrid_diagnostics(winners,pool,eligible,L,B,metric),
    status=status,witness=witness,screening=NULL))
  cfg <- list(seed=seed,workers=workers,sat_max3=1,RSI_min=rsi_min,
    screening=list(trivial_quadrants=TRUE,anchor_rule='rsi_and_algo',audit_fraction=o$audit_fraction,
      stability_min=max(algo_min,o$anchor_stability),audit_scheme='shared_persistent'))
  metadata <- list();audit <- NULL;cache <- new.env(parent=emptyenv())
  values <- data.frame(Stab_algo=rep(NA_real_,M),Stab_boot=NA_real_,RSI=NA_real_)
  repeat {
    plan <- ceca_screen_plan(state,cfg,metadata,audit);audit <- plan$audit_state
    todo <- plan$next_indices;todo <- todo[eligible[todo]]
    # Ineligible cells must not block progress towards eligible anchors.
    if(!length(todo))todo <- which(eligible & !as.character(g$id)%in%names(metadata) &
      plan$table$status=='pending' & !plan$table$awaiting_quadrant_audit)
    if(!length(todo))break
    z <- CECgrid_diagnostics(winners,pool,seq_len(M)%in%todo,L,B,metric,cache)
    values[todo,] <- z[todo,]
    for(i in todo) {
      v <- data.frame(method='exhaustif',Stab_algo=z$Stab_algo[i],Stab_boot=z$Stab_boot[i],RSI=z$RSI[i])
      qualifies <- is.finite(v$RSI) && v$RSI>rsi_min && is.finite(v$Stab_algo) && v$Stab_algo>cfg$screening$stability_min
      metadata[[as.character(g$id[i])]] <- list(qualifies=qualifies,values=v,best=list(K=g$REO[i]))
    }
  }
  # A quadrant whose mandatory controls are outside the measurement policy is
  # never approved; measure its remaining eligible cells instead of extrapolating.
  remaining <- which(eligible & !as.character(g$id)%in%names(metadata) & plan$table$status!='extrapolated_trivial_unmeasured')
  if(length(remaining)) {
    z <- CECgrid_diagnostics(winners,pool,seq_len(M)%in%remaining,L,B,metric,cache)
    values[remaining,] <- z[remaining,]
  }
  skipped <- eligible & plan$table$status=='extrapolated_trivial_unmeasured'
  status[skipped] <- 'quadrant_unmeasured';witness[skipped] <- plan$table$witness[skipped]
  list(values=values,status=status,witness=witness,screening=list(coverage=plan$table,
    quadrants=plan$quadrants,rejected=plan$rejected,audit_state=audit))
}
