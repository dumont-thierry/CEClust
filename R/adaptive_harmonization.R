ceca_harmonize <- function(state,cfg,seed_ids=NULL,tag='') {
  tm<-q_time();g<-state$mesh$leaves;adj<-state$mesh$adj;n<-nrow(g);X<-cfg$X
  parts<-ceca_current(state);backend<-ca_ns('CECprepare_backend_data')(as.data.frame(X),familyType='gaussVector')
  if(is.null(seed_ids))active<-seq_len(n) else {
    changed<-match(seed_ids,g$id);changed<-changed[!is.na(changed)]
    active<-sort(unique(c(changed,unlist(adj[changed],use.names=FALSE))))
  }
  queue<-active;queued<-seq_len(n)%in%active;cursor<-1L;visits<-integer(n)
  comparisons<-0L;logs<-list();failures<-0L;descent_time<-0
  while(cursor<=length(queue)) {
    i<-queue[cursor];cursor<-cursor+1L;queued[i]<-FALSE;visits[i]<-visits[i]+1L
    candidates<-c(i,adj[[i]])
    costs<-vapply(candidates,function(j)ca_profile(parts[[j]]$stats,g$C[i],g$lambda[i],ncol(X))$J,0)
    comparisons<-comparisons+length(candidates)-1L
    j<-candidates[which.min(costs)];before<-parts[[i]]$J
    if(!ca_better(min(costs),before,cfg$settings$tolerance))next
    selected<-q_part(X,parts[[j]]$labels,g$C[i],g$lambda[i],parts[[j]]$stats)
    set.seed(ca_seed(cfg$seed,'harmony',tag,g$id[i],length(logs)+1L))
    f<-q_descent(X,g$C[i],g$lambda[i],selected$labels,cfg$settings,backend)
    descent_time<-descent_time+f$time['elapsed'];failures<-failures+!f$valid
    if(f$valid && ca_better(f$part$J,selected$J,cfg$settings$tolerance))selected<-f$part
    stopifnot(selected$J<=min(costs)+1e-7,selected$J<before)
    parts[[i]]<-selected
    logs[[length(logs)+1L]]<-data.frame(cell=g$id[i],source=g$id[j],before=before,
      transferred=min(costs),after=selected$J,fit_valid=f$valid,
      fit_elapsed=f$time['elapsed'],fit_user=f$time['user'],fit_system=f$time['system'],
      error=f$error,row.names=NULL)
    stopifnot(length(logs)<=100L*n)
    add<-adj[[i]][!queued[adj[[i]]]]
    queue<-c(queue,add);queued[add]<-TRUE
  }
  for(i in seq_len(n))state$parts[[as.character(g$id[i])]]<-parts[[i]]
  list(state=state,time=q_duration(tm),comparisons=comparisons,visits=sum(visits),
    visited_cells=sum(visits>0),seed_cells=length(active),changes=length(logs),fit_failures=failures,
    descent_elapsed=descent_time,updates=if(length(logs))do.call(rbind,logs) else data.frame())
}
