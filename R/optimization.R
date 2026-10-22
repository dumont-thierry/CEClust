# Public optimization choices share the same hard, density-constrained score.
CECwith_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  had <- exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get('.Random.seed', envir = .GlobalEnv)
  on.exit(if (had) assign('.Random.seed', old, envir = .GlobalEnv) else
    if (exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE))
      rm('.Random.seed', envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  force(code)
}

CECinteger <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) ||
      x < minimum || x > .Machine$integer.max) stop(name, " must be an integer >= ", minimum, call. = FALSE)
  as.integer(x)
}

CECnumeric_data <- function(Z) {
  if (is.data.frame(Z) && !all(vapply(Z, is.numeric, TRUE)))
    stop("This optimizer requires numeric Gaussian data; use unpenalized CEM for categorical or mixed data.", call. = FALSE)
  if (is.numeric(Z) && is.null(dim(Z))) Z <- matrix(Z, ncol = 1L)
  X <- as.matrix(Z)
  if (!is.numeric(X) || !nrow(X) || !ncol(X) || any(!is.finite(X)))
    stop("Z must contain finite numeric observations with no missing values.", call. = FALSE)
  storage.mode(X) <- 'double'
  X
}

CECoptimization_control <- function(control, n, C, Nshots, Nloop, r0 = NULL) {
  defaults <- list(Rmax = if (is.null(r0)) min(20L, n) else r0,
    Ninit = Nshots, Cmax = C, max_iter = Nloop, tolerance = 1e-8, soften = .02,
    seed = NULL, rinit_max = 320L, min_child = 3L, confirmations = 2L,
    max_batches = 8L, saturation_stop = .05, orders = NULL)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop("optimizer_control must be a named list without duplicates.", call. = FALSE)
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) stop("Unknown optimizer_control: ", paste(unknown, collapse=', '), call. = FALSE)
  for (nm in names(control)) defaults[nm] <- control[nm]
  for (nm in c('Rmax','Ninit','max_iter','rinit_max','min_child','confirmations','max_batches'))
    defaults[[nm]] <- CECinteger(defaults[[nm]], nm)
  defaults$Rmax <- min(n, defaults$Rmax)
  if (!is.null(defaults$orders)) {
    z <- defaults$orders
    if (!is.numeric(z) || !length(z) || any(!is.finite(z)) || any(z!=floor(z)) ||
        any(z<1 | z>n) || anyDuplicated(z))
      stop('orders must contain distinct integer orders between 1 and n.',call.=FALSE)
    defaults$orders <- sort(as.integer(z))
    if (!'Rmax' %in% names(control)) defaults$Rmax <- max(defaults$orders)
    if (max(defaults$orders)>defaults$Rmax) stop('orders must not exceed Rmax.',call.=FALSE)
  }
  if (!is.null(defaults$seed)) defaults$seed <- CECinteger(defaults$seed, 'seed', 0L)
  for (nm in c('Cmax','tolerance','soften','saturation_stop')) {
    x <- defaults[[nm]]
    if (!is.numeric(x) || length(x)!=1L || !is.finite(x)) stop(nm, ' must be finite.', call. = FALSE)
  }
  if (defaults$Cmax < C*(1-16*.Machine$double.eps) || defaults$Cmax <= 0)
    stop('Cmax must be positive and at least C.', call. = FALSE)
  # Mesh centres can differ from supplied values by a few rounding bits.
  defaults$Cmax <- max(defaults$Cmax,C)
  if (defaults$tolerance < 0 || defaults$soften <= 0 || defaults$soften >= 1 ||
      defaults$saturation_stop < 0 || defaults$saturation_stop > 1)
    stop('Invalid tolerance, soften or saturation_stop.', call. = FALSE)
  defaults
}

CECbank_initial <- function(X, r, shot) {
  n <- nrow(X)
  if (r == 1L) return(rep(1L, n))
  balanced <- function() sample(rep(seq_len(r), length.out = n))
  if (shot %% 5L == 0L) return(balanced())
  centers <- sample.int(n, 1L)
  dist2 <- rep(Inf, n)
  for (j in 2:r) {
    D <- sweep(X, 2, X[tail(centers, 1L), ], '-')
    dist2 <- pmin(dist2, rowSums(D*D)); dist2[centers] <- 0
    remaining <- setdiff(seq_len(n),centers)
    next_id <- if (is.finite(sum(dist2)) && sum(dist2) > 0)
      sample.int(n, 1L, prob = dist2) else remaining[sample.int(length(remaining),1L)]
    centers <- c(centers, next_id)
  }
  z <- tryCatch(suppressWarnings(stats::kmeans(X, centers = X[centers,,drop=FALSE],
    iter.max=50L, algorithm='Lloyd')$cluster), error=function(e) balanced())
  if (length(unique(z)) != r) z <- balanced()
  ca_labels(z)
}

CECfit_em_bank <- function(X, control, Xref = X, checkpoint_dir=NULL, resume=FALSE) {
  CECwith_seed(control$seed, {
    signature <- digest::digest(list(X,Xref,control,as.character(utils::packageVersion('CEClust'))))
    if(!is.null(checkpoint_dir))dir.create(checkpoint_dir,recursive=TRUE,showWarnings=FALSE)
    candidates <- list(); starts <- list(); keys <- new.env(parent=emptyenv())
    # A hard one-group candidate remains defined for singular datasets.
    one <- list(local_labels=rep(1L,nrow(X)), labels=rep(1L,nrow(Xref)),
      local_stats=ca_stats(X,rep(1L,nrow(X))), stats=ca_stats(Xref,rep(1L,nrow(Xref))),
      params=NULL, loglik_trace=numeric(), stop='hard_trivial_baseline', converged=NA)
    candidates[[1L]] <- one
    assign(digest::digest(list(one$local_labels,one$labels)),1L,keys)
    orders <- if(is.null(control$orders))seq_len(control$Rmax) else control$orders
    for (r in orders) for (shot in seq_len(control$Ninit)) {
      file <- if(is.null(checkpoint_dir))NULL else file.path(checkpoint_dir,sprintf('start_%04d_%03d.rds',r,shot))
      cached <- !is.null(file) && isTRUE(resume) && file.exists(file)
      if(cached) {
        saved <- readRDS(file)
        if(!identical(saved$signature,signature))stop('EM start checkpoint signature mismatch.')
        init <- saved$init;fit <- saved$fit
        if(!is.null(saved$rng))assign('.Random.seed',saved$rng,envir=.GlobalEnv) else
          if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
      } else {
      init <- CECbank_initial(X, r, shot)
      fit <- tryCatch({
        e <- emmap_fit_cpp(X,init,control$Cmax,control$max_iter,control$tolerance,control$soften)
        if (any(diff(e$loglik_mean_trace)<0) || e$minimum_logdet_constraint_margin < -1e-8)
          stop('EM monotonicity or density-constraint check failed.')
        keep <- unique(e$labels)
        params <- list(means=e$means[keep,,drop=FALSE],cov=e$cov[,,keep,drop=FALSE],
          weights=as.numeric(e$weights[keep])/sum(e$weights[keep]))
        local <- ca_labels(e$labels)
        if (!identical(local,as.integer(emmap_predict_cpp(X,params$means,params$cov,params$weights))))
          stop('MAP pruning changed the training allocation.')
        labels <- ca_labels(emmap_predict_cpp(Xref,params$means,params$cov,params$weights))
        list(local_labels=local,labels=labels,local_stats=ca_stats(X,local),stats=ca_stats(Xref,labels),
          params=params,loglik_trace=e$loglik_mean_trace,stop=e$stop,converged=e$converged,
          minimum_constraint_margin=e$minimum_logdet_constraint_margin,soft_K=e$soft_K)
      }, error=function(e)e)
      if(!is.null(file))CECwrite_rds_atomic(list(signature=signature,init=init,fit=fit,
        rng=if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))get('.Random.seed',envir=.GlobalEnv) else NULL),file)
      }
      ok <- !inherits(fit,'error')
      id <- NA_integer_
      if (ok) {
        key <- digest::digest(list(fit$local_labels,fit$labels))
        if (exists(key,keys,inherits=FALSE)) id <- get(key,keys) else {
          id <- length(candidates)+1L; candidates[[id]] <- fit; assign(key,id,keys)
        }
        # Replace the baseline with the actual EM fit when it gives one group.
        if (is.null(candidates[[id]]$params)) candidates[[id]] <- fit
      }
      starts[[length(starts)+1L]] <- list(r=r,shot=shot,initial_labels=init,candidate=id,valid=ok,
        error=if(ok)'' else conditionMessage(fit),
        stop=if(ok)fit$stop else 'fit_failure',trace=if(ok)fit$loglik_trace else numeric())
    }
    structure(list(candidates=candidates,starts=starts,control=control,
      valid_starts=sum(vapply(starts,`[[`,TRUE,'valid')), n=nrow(X),d=ncol(X)),class='CEC_candidate_bank')
  })
}

CECscore_candidates <- function(parts, C, lambda, d, penalty, stats_name='stats') {
  profile <- CECsolver_environment(penalty)$ca_profile
  scores <- lapply(parts,function(x)profile(x[[stats_name]],C,lambda,d))
  best <- which.min(vapply(scores,`[[`,0,'J'))
  list(best=best,profile=scores[[best]],scores=vapply(scores,`[[`,0,'J'))
}

CECfit_from_partition <- function(X, labels, C, lambda, penalty, optimization,
                                  solver=NULL, bank=NULL, selected_candidate=NULL) {
  labels <- ca_labels(labels)
  backend <- CECprepare_backend_data(as.data.frame(X),familyType='gaussVector')
  params <- CECoptParam_fast_from_clusters(backend,labels,lambda,C,Inf)
  profile <- CECsolver_environment(penalty)$ca_profile(ca_stats(X,labels),C,lambda,ncol(X))
  phi <- CECclusters_to_phi(labels,max(labels));params$phi <- phi
  params$optimization <- optimization;params$penalty <- penalty
  out <- list(phi=phi,params=params,Hphi=profile$J/lambda,REO=max(labels),clusters=labels,
    objective=profile$J,entropy=profile$J_raw,penalty_value=profile$BIC_penalty_relative,
    saturation=profile$saturation,optimization=optimization,penalty=penalty,
    criterion_scale='H_lambda + relative BIC; Hphi = objective/lambda',
    solver=solver,bank=bank,selected_candidate=selected_candidate)
  class(out) <- c('CEC_fit','list')
  out
}

# Public help: man/CECclassif.Rd (maintained by hand).
CECclassif <- function(Z,lambda=1,C=1,Cquali=Inf,r0=NULL,Nshots=100,Nloop=1000,
  familyType='gaussAndDiscreteVector',sizeMaxOutlier=0,autoRegroupOutliers=FALSE,
  displayRemainingTime=FALSE,focus=NULL,backend_data=NULL,n_cores=1L,
  shot_strategy=c('current','independent'),optimization=c('cem','em_bank'),
  penalty=c('none','bic'),optimizer_control=list()) {
  optimization <- match.arg(optimization);penalty <- match.arg(penalty)
  shot_strategy <- match.arg(shot_strategy)
  if (!is.numeric(lambda) || length(lambda)!=1L || !is.finite(lambda) || lambda<=0 ||
      !is.numeric(C) || length(C)!=1L || !is.finite(C) || C<=0)
    stop('C and lambda must be finite and strictly positive.',call.=FALSE)
  Nshots <- CECinteger(Nshots,'Nshots');Nloop <- CECinteger(Nloop,'Nloop')
  if(!is.null(r0))r0 <- CECinteger(r0,'r0')
  n_cores <- CECinteger(n_cores,'n_cores')
  if (!is.list(optimizer_control) || (length(optimizer_control) &&
      (is.null(names(optimizer_control)) || any(!nzchar(names(optimizer_control))) || anyDuplicated(names(optimizer_control)))))
    stop('optimizer_control must be a named list without duplicates.',call.=FALSE)
  if (optimization=='cem' && penalty=='none' && !length(setdiff(names(optimizer_control),'seed'))) {
    seed <- optimizer_control$seed
    if (!is.null(seed)) seed <- CECinteger(seed,'seed',0L)
    out <- CECwith_seed(seed, CECclassif_legacy(Z,lambda,C,Cquali,r0,Nshots,Nloop,
      familyType,sizeMaxOutlier,autoRegroupOutliers,displayRemainingTime,focus,backend_data,n_cores,shot_strategy))
    if (!is.null(out)) {
      out$optimization <- optimization;out$penalty <- penalty
      out$objective <- lambda*out$Hphi;out$entropy <- out$objective;out$penalty_value <- 0
      out$params$optimization <- optimization;out$params$penalty <- penalty
    }
    return(out)
  }
  if (!familyType %in% c('gaussVector','gaussUniv','gaussAndDiscreteVector'))
    stop('This optimizer/penalty combination requires a Gaussian family.',call.=FALSE)
  if(n_cores!=1L)
    stop('For these engines, parallelize independent grid repetitions with CECfitBoundGrid(n_cores=...). A single CECclassif call uses n_cores=1.',call.=FALSE)
  X <- CECnumeric_data(Z)
  if (!is.null(focus) || sizeMaxOutlier!=0 || autoRegroupOutliers || Cquali!=Inf)
    stop('focus, outlier regrouping and qualitative floors are not supported by this optimizer.',call.=FALSE)
  control <- CECoptimization_control(optimizer_control,nrow(X),C,Nshots,Nloop,r0)
  if (optimization=='em_bank') {
    bank <- CECfit_em_bank(X,control)
    selection <- CECscore_candidates(bank$candidates,C,lambda,ncol(X),penalty)
    p <- bank$candidates[[selection$best]]
    return(CECfit_from_partition(X,p$local_labels,C,lambda,penalty,optimization,
      solver=list(stop=p$stop,converged=p$converged),bank=bank,selected_candidate=selection$best))
  }
  settings <- q_settings()
  settings$rinit <- control$Rmax;settings$shots <- control$Ninit;settings$Nloop <- control$max_iter
  settings$rinit_max <- max(control$rinit_max,settings$rinit)
  settings$min_child <- control$min_child;settings$confirmations <- control$confirmations
  settings$max_batches <- control$max_batches;settings$tau <- control$saturation_stop
  settings$tolerance <- control$tolerance
  seed <- if(is.null(control$seed))sample.int(.Machine$integer.max,1L) else control$seed
  fit <- CECwith_seed(seed,CECsolver_environment(penalty)$q_fit(X,C,lambda,seed,settings))
  if (fit$status=='fit_failure') stop('All CEM starts failed; inspect the input and constraints.',call.=FALSE)
  CECfit_from_partition(X,fit$part$labels,C,lambda,penalty,optimization,solver=fit)
}
