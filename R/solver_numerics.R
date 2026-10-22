# Numerical routines consolidated from v052/v075; no runtime sourcing.
q_bic_beta <- function(n,d)(d+1)*(d+2)*log(n)/(4*n)
ca_ns <- function(x)getFromNamespace(x,'CEClust')
ca_seed <- function(...)strtoi(substr(digest::digest(list(...),algo='xxhash64'),1,7),16L)
ca_labels <- function(x)as.integer(match(x,unique(x)))
ca_sim <- function(a,b) {
  if(identical(ca_labels(a),ca_labels(b)))return(1)
  tab <- table(a,b);k<-max(dim(tab));m<-matrix(0,k,k)
  m[seq_len(nrow(tab)),seq_len(ncol(tab))]<-tab
  assignment<-as.integer(clue::solve_LSAP(m,maximum=TRUE))
  sum(m[cbind(seq_len(k),assignment)])/length(a)
}
ca_stats <- function(X,labels) {
  groups<-split(seq_len(nrow(X)),labels);d<-ncol(X)
  do.call(rbind,lapply(seq_along(groups),function(k){
    a<-X[groups[[k]],,drop=FALSE]
    shifted<-sweep(a,2,a[1,],'-');centered<-sweep(shifted,2,colMeans(shifted),'-')
    S<-crossprod(centered)/nrow(a);ev<-eigen(S,symmetric=TRUE,only.values=TRUE)$values
    rank<-sum(ev>8*d*.Machine$double.eps*max(abs(ev)))
    ld<-if(rank==d)2*sum(log(diag(chol((S+t(S))/2)))) else -Inf
    data.frame(cluster=k,n=nrow(a),mass=nrow(a)/nrow(X),rank=rank,logdet=ld,
      eigmin=min(ev),eigmax=max(ev))
  }))
}
ca_profile <- function(stats,C,lambda,d) {
  logdelta<- -d*log(2*pi)-2*log(C)
  singular<-stats$rank<d
  ls<-pmax(0,(logdelta-stats$logdet)/d)
  conditional<-d/2*log(2*pi)+pmax(stats$logdet,logdelta)/2+d/2*exp(-ls)
  conditional[singular]<- -log(C)
  h<- -sum(stats$mass*log(stats$mass))
  saturated<-singular | stats$logdet<=logdelta+1e-10*max(1,abs(logdelta))
  list(J=lambda*h+sum(stats$mass*conditional),h=h,K=nrow(stats),
    saturation=sum(stats$mass[saturated]),saturated=saturated)
}
ca_better <- function(new,old,tolerance=1e-8) {
  is.finite(new) & (!is.finite(old) | new<old-tolerance*pmax(1,abs(old)))
}
ca_subdivide <- function(X,labels,target,min_child=6L) {
  pieces<-split(seq_len(nrow(X)),ca_labels(labels));sizes<-lengths(pieces)
  capacity<-pmax(1L,sizes%/%min_child)
  if(target<length(pieces) || target>sum(capacity))stop('Subdivision target cannot be constructed.')
  quota<-rep(1L,length(pieces))
  # Split every sufficiently large parent once before assigning extra children.
  eligible<-which(capacity>=2L)
  eligible<-eligible[order(-sizes[eligible],eligible)]
  take<-head(eligible,target-sum(quota));quota[take]<-2L
  while(sum(quota)<target) {
    score<-sizes/quota;score[quota>=capacity]<- -Inf
    k<-which.max(score);quota[k]<-quota[k]+1L
  }
  result<-integer(nrow(X));offset<-0L;fallback<-0L
  for(k in seq_along(pieces)) {
    ids<-pieces[[k]];q<-quota[k];A<-X[ids,,drop=FALSE]
    if(q==1L)child<-rep(1L,length(ids)) else {
      child<-NULL
      for(attempt in 1:10) {
        km<-tryCatch(suppressWarnings(kmeans(A,centers=q,iter.max=50L,nstart=1L,algorithm='Lloyd')),error=function(e)NULL)
        if(!is.null(km)&&length(unique(km$cluster))==q&&min(tabulate(km$cluster,q))>=min_child) {
          child<-km$cluster;break
        }
      }
      if(is.null(child)) {
        ordering<-order(drop(A%*%rnorm(ncol(X))),runif(nrow(A)))
        child<-integer(nrow(A));child[ordering]<-floor((seq_len(nrow(A))-1L)*q/nrow(A))+1L
        fallback<-fallback+1L
      }
    }
    result[ids]<-offset+child;offset<-offset+q
  }
  stopifnot(length(unique(result))==target,
    all(vapply(split(labels,result),function(x)length(unique(x))==1L,TRUE)))
  list(labels=ca_labels(result),fallback=fallback,allocation=data.frame(parent=seq_along(pieces),n=sizes,children=quota))
}

q_time <- function() c(unname(proc.time()[1:2]),as.numeric(Sys.time()))
q_duration <- function(t) setNames(q_time()-t,c('user','system','elapsed'))
q_part <- function(X,labels,C,lambda,stats=NULL) {
  labels <- ca_labels(labels)
  if(is.null(stats))stats <- ca_stats(X,labels)
  c(list(labels=labels,stats=stats),ca_profile(stats,C,lambda,ncol(X)))
}
q_settings <- function()list(rinit=10L,rinit_max=320L,shots=5L,Nloop=100L,
  min_child=3L,tau=.05,confirmations=2L,max_batches=8L,tolerance=1e-8)

q_descent <- function(X,C,lambda,labels,settings,backend) {
  tm <- q_time();warnings <- character()
  out <- tryCatch(withCallingHandlers({
    traces<-list();deletions<-list();fit_elapsed<-0;proposal_elapsed<-0;deletion_elapsed<-0
    evaluate<-function(lab) {
      te<-q_time();p<-q_part(X,lab,C,lambda)
      params<-ca_ns('CECoptParam_fast_from_clusters')(backend,p$labels,lambda,C,Inf)
      logdens<-ca_ns('CECcompute_logdens_matrix_fast')(backend,params,lambda)
      from_density<-lambda*p$h-mean(logdens[cbind(seq_len(nrow(X)),p$labels)])
      stopifnot(abs(from_density-p$J_raw)<1e-6)
      fit_elapsed<<-fit_elapsed+q_duration(te)['elapsed']
      list(part=p,params=params,logdens=logdens)
    }
    current<-evaluate(labels);initial<-current$part
    append_trace<-function(move,loss=NA_real_) {
      p<-current$part;traces[[length(traces)+1L]]<<-data.frame(step=length(traces),move=move,
        K=p$K,CEC=p$J_raw,penalty_relative=p$BIC_penalty_relative,J=p$J,
        full_BIC_penalty=p$BIC_penalty,full_J=p$J_BIC_full,frozen_deletion_loss=loss)
    }
    append_trace('initial');beta<-q_bic_beta(nrow(X),ncol(X));reason<-'active'
    repeat {
      cem_stop<-'iteration_limit'
      for(iteration in seq_len(settings$Nloop)) {
        tp<-q_time()
        proposal<-ca_ns('CECproposal_logdens')(current$logdens,current$params)
        scores<-sweep(proposal,2,lambda*log(current$params$nu),'+')
        next_labels<-ca_labels(max.col(scores,ties.method='first'))
        proposal_elapsed<-proposal_elapsed+q_duration(tp)['elapsed']
        if(identical(next_labels,current$part$labels)){cem_stop<-'unchanged_partition';break}
        candidate<-evaluate(next_labels)
        stopifnot(candidate$part$K<=current$part$K)
        if(!ca_better(candidate$part$J,current$part$J,settings$tolerance)) {
          cem_stop<-'no_penalized_CEM_improvement';break
        }
        current<-candidate;append_trace('CEM')
      }
      if(cem_stop=='iteration_limit'){reason<-'CEM_iteration_budget';break}
      td<-q_time();tol<-settings$tolerance*max(1,abs(current$part$J))
      deletion<-q_deletion_plan(current$logdens,current$params$nu,current$part$labels,lambda,beta,tol)
      deletion_elapsed<-deletion_elapsed+q_duration(td)['elapsed']
      deletions[[length(deletions)+1L]]<-list(K=current$part$K,costs=deletion$costs,
        selected=deletion$selected,loss=deletion$loss,beta=beta,CEM_stop=cem_stop)
      if(is.null(deletion$candidate)){reason<-'no_improving_frozen_deletion';break}
      before<-current$part;candidate<-evaluate(deletion$candidate)
      stopifnot(candidate$part$K==before$K-1L,
        candidate$part$J-before$J<=deletion$loss-beta+1e-6,
        ca_better(candidate$part$J,before$J,settings$tolerance))
      current<-candidate;append_trace('delete_component',deletion$loss)
    }
    trace<-do.call(rbind,traces)
    stopifnot(all(diff(trace$J)<0),current$part$J<=initial$J+1e-10)
    list(part=current$part,valid=TRUE,error='',initial_J=initial$J,
      steps=nrow(trace)-1L,stop=reason,descent_trace=trace,deletion_plans=deletions,
      deletion_count=sum(trace$move=='delete_component'),
      stage_timings=c(refit=fit_elapsed,CEM_proposal=proposal_elapsed,deletion=deletion_elapsed),
      max_increase=if(nrow(trace)>1L)max(diff(trace$J)) else 0)
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e)list(valid=FALSE,error=conditionMessage(e)))
  out$time <- q_duration(tm);out$warnings <- paste(unique(warnings),collapse=' | ');out
}

# One complete algorithm run. Warm starts are allowed only in the pilot.
# The independent measurement jobs always pass warm=NULL.
q_fit <- function(X,C,lambda,seed,settings=q_settings(),warm=NULL,backend=NULL) {
  tm <- q_time();set.seed(seed);X <- as.matrix(X);prep <- q_time()
  if(is.null(backend))backend <- ca_ns('CECprepare_backend_data')(as.data.frame(X),familyType='gaussVector')
  prep_time <- q_duration(prep)
  best <- q_part(X,rep(1L,nrow(X)),C,lambda)
  trace <- list();starts <- list();valid_total <- 0L;descent_traces<-list();deletion_plans<-list()
  run <- function(labels,phase,batch,shot) {
    ans <- q_descent(X,C,lambda,labels,settings,backend)
    starts[[length(starts)+1L]] <<- data.frame(phase=phase,batch=batch,shot=shot,
      initial_K=length(unique(labels)),valid=ans$valid,
      final_K=if(ans$valid)ans$part$K else NA_integer_,
      J=if(ans$valid)ans$part$J else NA_real_,
      CEC=if(ans$valid)ans$part$J_raw else NA_real_,
      BIC_penalty_relative=if(ans$valid)ans$part$BIC_penalty_relative else NA_real_,
      deletions=if(ans$valid)ans$deletion_count else NA_integer_,
      descent_stop=if(ans$valid)ans$stop else 'failure',
      max_J_increase=if(ans$valid)ans$max_increase else NA_real_,
      refit_elapsed=if(ans$valid)unname(ans$stage_timings['refit']) else NA_real_,
      CEM_proposal_elapsed=if(ans$valid)unname(ans$stage_timings['CEM_proposal']) else NA_real_,
      deletion_elapsed=if(ans$valid)unname(ans$stage_timings['deletion']) else NA_real_,
      user=ans$time['user'],system=ans$time['system'],elapsed=ans$time['elapsed'],
      error=ans$error,warnings=ans$warnings,row.names=NULL)
    if(ans$valid) {
      descent_traces[[length(starts)]]<<-ans$descent_trace
      deletion_plans[[length(starts)]]<<-ans$deletion_plans
      valid_total <<- valid_total+1L
      if(ca_better(ans$part$J,best$J,settings$tolerance))best <<- ans$part
    }
    ans$valid
  }
  if(!is.null(warm)) {
    p <- q_part(X,warm,C,lambda)
    if(ca_better(p$J,best$J,settings$tolerance))best <- p
    run(warm,'warm',0L,0L)
  }
  initial_r <- min(settings$rinit,max(1L,nrow(X)%/%settings$min_child))
  subdivision_time <- c(user=0,system=0,elapsed=0)
  for(s in seq_len(settings$shots)) {
    ts <- q_time();z <- ca_subdivide(X,rep(1L,nrow(X)),initial_r,settings$min_child)
    subdivision_time <- subdivision_time+q_duration(ts)
    run(z$labels,'initial',0L,s)
  }
  attempted <- initial_r;stale <- 0L;batches <- 0L;status <- 'active'
  repeat {
    if(valid_total==0L){status<-'fit_failure';break}
    if(best$saturation>settings$tau){status<-'saturated';break}
    if(stale>=settings$confirmations){status<-'plateau';break}
    if(batches>=settings$max_batches){status<-'budget_limit';break}
    capacity <- sum(pmax(1L,as.integer(table(best$labels))%/%settings$min_child))
    limit <- min(settings$rinit_max,capacity)
    if(limit<=best$K){status<-'subdivision_exhausted';break}
    target <- min(limit,max(2L*best$K,2L*attempted))
    mode <- if(target>attempted)'grow' else 'confirm'
    batches <- batches+1L;old <- best;ok <- 0L
    for(s in seq_len(settings$shots)) {
      ts <- q_time();z <- ca_subdivide(X,old$labels,target,settings$min_child)
      subdivision_time <- subdivision_time+q_duration(ts)
      ok <- ok+run(z$labels,mode,batches,s)
    }
    improved <- ca_better(best$J,old$J,settings$tolerance)
    stale <- if(improved)0L else stale+as.integer(ok>0L)
    attempted <- target
    trace[[batches]] <- data.frame(batch=batches,target=target,mode=mode,
      improved=improved,J=best$J,K=best$K,saturation=best$saturation,stale=stale,valid=ok)
    if(!ok){status<-'fit_failure';break}
  }
  pt <- q_time()
  params <- ca_ns('CECoptParam_fast_from_clusters')(backend,best$labels,lambda,C,Inf)
  params$phi <- NULL
  parameter_time <- q_duration(pt)
  list(part=best,params=params,status=status,batches=batches,stale=stale,
    attempted=attempted,seed=seed,warm=!is.null(warm),
    starts=do.call(rbind,starts),trace=if(length(trace))do.call(rbind,trace) else data.frame(),
    descent_traces=descent_traces,deletion_plans=deletion_plans,
    time=q_duration(tm),backend_time=prep_time,subdivision_time=subdivision_time,parameter_time=parameter_time)
}

q_deletion_plan<-function(raw_logdens,weights,labels,lambda,beta,tolerance) {
 n<-length(labels);K<-length(weights)
 if(K<=1L)return(list(candidate=NULL,costs=numeric(),selected=NA_integer_,loss=Inf))
 scores<-sweep(raw_logdens,2,lambda*log(weights),'+')
 at<-cbind(seq_len(n),labels);own<-scores[at]
 stopifnot(all(is.finite(own)))
 scores[at]<- -Inf
 alternative<-max.col(scores,ties.method='first')
 best_other<-scores[cbind(seq_len(n),alternative)]
 losses<-(own-best_other)/n
 costs<-as.numeric(rowsum(losses,labels,reorder=TRUE))
 stopifnot(length(costs)==K)
 chosen<-which.min(costs)
 if(!is.finite(costs[chosen]) || costs[chosen]>=beta-tolerance)
  return(list(candidate=NULL,costs=costs,selected=NA_integer_,loss=Inf))
 candidate<-labels;affected<-labels==chosen;candidate[affected]<-alternative[affected]
 list(candidate=ca_labels(candidate),costs=costs,selected=chosen,loss=costs[chosen])
}

# Prediction only: frozen learned parameters, no density floor, no refitting.
gaussian_prediction_logdens<-function(X,mu,S) {
  X<-as.matrix(X);S<-as.matrix(S);d<-ncol(X)
  stopifnot(all(is.finite(X)),length(mu)==d,all(is.finite(mu)),all(dim(S)==c(d,d)),all(is.finite(S)))
  support<-attr(S,'CEC_singular_support',exact=TRUE)
  if(is.null(support)) {
    R<-tryCatch(chol((S+t(S))/2),error=function(e)NULL)
    if(!is.null(R)) {
      centered<-sweep(X,2,mu)
      v<-backsolve(R,t(centered),transpose=TRUE)
      return(-d/2*log(2*pi)-sum(log(diag(R)))-colSums(v*v)/2)
    }
    eig<-eigen((S+t(S))/2,symmetric=TRUE)
    tol<-8*d*.Machine$double.eps*max(abs(eig$values))
    if(min(eig$values)< -tol)stop('Invalid covariance: negative eigenvalue.')
    null<-eig$values<=tol
    if(!any(null))stop('Full-rank covariance factorization failed.')
    support<-list(null_basis=eig$vectors[,null,drop=FALSE],support_tolerance=sqrt(tol))
  }
  residual<-sqrt(rowSums((sweep(X,2,mu) %*% support$null_basis)^2))
  scale<-pmax(apply(abs(X),1,max),max(abs(mu)))
  tolerance<-support$support_tolerance+64*d*.Machine$double.eps*scale
  # User's prediction convention: infinity on the affine support, zero outside.
  ifelse(residual<=tolerance,Inf,-Inf)
}

classify_log_scores<-function(scores,u) {
  scores<-as.matrix(scores);n<-nrow(scores);G<-ncol(scores)
  stopifnot(G>=1L,length(u)==n,all(u>=0 & u<1),!anyNA(scores))
  isinf<-is.infinite(scores)&scores>0
  numinf<-rowSums(isinf)
  allzero<-rowSums(is.finite(scores))==0L & numinf==0L
  labels<-max.col(scores,ties.method='first')
  for(i in which(numinf>0L)) {
    candidates<-which(isinf[i,])
    labels[i]<-candidates[1L+floor(u[i]*length(candidates))]
  }
  labels[allzero]<-1L+floor(u[allzero]*G)
  list(labels=as.integer(labels),single_infinite=sum(numinf==1L),
    multiple_infinite=sum(numinf>1L),all_zero=sum(allzero))
}

predict_frozen<-function(X,params,method,seed) {
  before<-serialize(params,NULL)
  if(method=='CEC') {
    pro<-params$nu;mu<-t(params$m);covs<-params$Sigma;lambda<-params$lambda
  } else {
    pro<-params$pro;mu<-params$means;covs<-params$covs;lambda<-1
  }
  stopifnot(length(pro)==ncol(mu),length(pro)==length(covs),all(is.finite(pro)),
    all(pro>=0),sum(pro)>0,is.finite(lambda),lambda>0)
  active<-which(pro>0);n<-nrow(X)
  logdens<-vapply(active,function(k)gaussian_prediction_logdens(X,mu[,k],covs[[k]]),numeric(n))
  logdens<-matrix(logdens,n,length(active))
  scores<-sweep(logdens/lambda,2,log(pro[active]),'+')
  set.seed(as.integer(seed));u<-runif(n)
  result<-classify_log_scores(scores,u)
  result$labels<-active[result$labels]
  result$seed<-seed;result$K<-length(unique(result$labels));result$lambda<-lambda
  result$active_clusters<-active
  stopifnot(identical(before,serialize(params,NULL)),!anyNA(result$labels))
  result
}


CECsolver_environment <- function(penalty) {
  ns <- environment(CECsolver_environment)
  env <- new.env(parent = ns)
  for (nm in c("q_part", "q_descent", "q_fit")) {
    f <- get(nm, envir = ns)
    environment(f) <- env
    assign(nm, f, envir = env)
  }
  env$q_bic_beta <- function(n, d) {
    if (penalty == "bic") (d + 1) * (d + 2) * log(n) / (4 * n) else 0
  }
  env$ca_profile <- function(stats, C, lambda, d) {
    p <- ca_profile(stats, C, lambda, d)
    n <- sum(stats$n)
    p$J_raw <- p$J
    p$parameter_count <- p$K * (d + 1) * (d + 2) / 2 - 1
    p$BIC_penalty <- if (penalty == "bic") p$parameter_count * log(n) / (2*n) else 0
    p$BIC_penalty_relative <- env$q_bic_beta(n, d) * (p$K - 1)
    p$J_BIC_full <- p$J_raw + p$BIC_penalty
    p$J <- p$J_raw + p$BIC_penalty_relative
    p
  }
  env
}
