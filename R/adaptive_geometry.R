# Imported from Moteurs/CEClust_leaf_patch_v060/mesh.R
# Active rectangular leaves: C uses log coordinates, lambda linear coordinates.
# Direct neighbours share an edge segment of positive length, never just a corner.
ceca_connect <- function(mesh) {
  g <- mesh$leaves;n <- nrow(g);tol <- 1e-10;adj <- vector('list',n)
  for(i in seq_len(n)) {
    v <- (abs(g$u0-g$u1[i])<tol | abs(g$u1-g$u0[i])<tol) &
      pmin(g$v1,g$v1[i])-pmax(g$v0,g$v0[i])>tol
    h <- (abs(g$v0-g$v1[i])<tol | abs(g$v1-g$v0[i])<tol) &
      pmin(g$u1,g$u1[i])-pmax(g$u0,g$u0[i])>tol
    adj[[i]] <- which(v|h)
  }
  edges <- do.call(rbind,lapply(seq_len(n),function(i) {
    j<-adj[[i]][adj[[i]]>i];if(length(j))cbind(i=i,j=j) else NULL
  }))
  if(is.null(edges))edges <- matrix(integer(),0,2,dimnames=list(NULL,c('i','j')))
  g$C <- exp((g$u0+g$u1)/2);g$lambda <- (g$v0+g$v1)/2
  mesh$leaves <- g;mesh$adj <- adj;mesh$edges <- edges;mesh
}
ceca_mesh <- function(C,lambda) {
  borders <- function(x)c(x[1]-(x[2]-x[1])/2,(head(x,-1)+tail(x,-1))/2,
    tail(x,1)+(tail(x,1)-tail(x,2)[1])/2)
  u <- borders(log(C));v <- borders(lambda)
  stopifnot(min(v)>0,length(C)>=2,length(lambda)>=2)
  p <- expand.grid(i=seq_along(C),j=seq_along(lambda));n<-nrow(p)
  g <- data.frame(id=seq_len(n),parent=NA_integer_,depth=0L,
    u0=u[p$i],u1=u[p$i+1],v0=v[p$j],v1=v[p$j+1])
  ceca_connect(list(leaves=g,archive=g[FALSE,],next_id=n+1L))
}
ceca_refine <- function(mesh,ids) {
  g<-mesh$leaves;children<-list();next_id<-mesh$next_id
  for(i in ids) {
    p<-g[i,];u<-c(p$u0,(p$u0+p$u1)/2,p$u1);v<-c(p$v0,(p$v0+p$v1)/2,p$v1)
    for(a in 1:2)for(b in 1:2) {
      children[[length(children)+1L]]<-data.frame(id=next_id,parent=p$id,depth=p$depth+1L,
        u0=u[a],u1=u[a+1],v0=v[b],v1=v[b+1]);next_id<-next_id+1L
    }
  }
  child<-do.call(rbind,children);cols<-names(child)
  out<-ceca_connect(list(leaves=rbind(g[-ids,cols,drop=FALSE],child),
    archive=rbind(mesh$archive[,cols,drop=FALSE],g[ids,cols,drop=FALSE]),next_id=next_id))
  stopifnot(nrow(out$leaves)==nrow(g)+3L*length(ids))
  list(mesh=out,new_ids=child$id,removed_ids=g$id[ids])
}
ceca_border <- function(mesh) {
  g<-mesh$leaves;tol<-1e-10
  list(left=which(abs(g$u0-min(g$u0))<tol),right=which(abs(g$u1-max(g$u1))<tol),
    lower=which(abs(g$v0-min(g$v0))<tol),upper=which(abs(g$v1-max(g$v1))<tol))
}
ceca_margins <- function(mesh,parts,tau) {
  s<-vapply(parts,`[[`,0,'saturation');K<-vapply(parts,`[[`,0L,'K');b<-ceca_border(mesh)
  terminal<-K==1L | s>tau
  okay<-c(left=all(terminal[b$left]),right=all(terminal[b$right]),
    lower=all(s[b$lower]>tau),upper=all(K[b$upper]==1L))
  list(observed=okay,needs_extension=!okay,border=b)
}
ceca_expand <- function(mesh,side) {
  g<-mesh$leaves;ids<-ceca_border(mesh)[[side]];added<-list();nid<-mesh$next_id
  for(i in ids) {
    p<-g[i,];t<-p;t$parent<-p$id;t$id<-nid;nid<-nid+1L;t$depth<-0L
    if(side=='left'){t$u1<-p$u0;t$u0<-p$u0-log(2)}
    if(side=='right'){t$u0<-p$u1;t$u1<-p$u1+log(2)}
    if(side=='lower'){t$v1<-p$v0;t$v0<-p$v0/2}
    if(side=='upper'){t$v0<-p$v1;t$v1<-p$v1*2}
    added[[length(added)+1L]]<-t
  }
  added<-do.call(rbind,added)
  out<-ceca_connect(list(leaves=rbind(g,added),archive=mesh$archive,next_id=nid))
  list(mesh=out,new_ids=added$id,removed_ids=integer())
}

# Imported from Moteurs/CEClust_leaf_bic_patch_v075/refinement.R
ceca_partition_key<-function(p)digest::digest(ca_labels(p$labels),algo='xxhash64')
ceca_adaptive_options<-function(cfg) {
  defaults<-list(patience=2L,position_tolerance=.125,max_rounds=12L,max_depth=8L,
                 distinguish_partitions=TRUE)
  out<-modifyList(defaults,if(is.null(cfg$refinement))list() else cfg$refinement)
  stopifnot(out$patience>=1L,out$position_tolerance>0,out$max_rounds>=2L,out$max_depth>=1L)
  out
}

# Edge-sharing neighbours only. Both endpoints of every selected frontier are refined.
ceca_frontiers<-function(mesh,parts,tau,distinguish_partitions=TRUE) {
  sat<-vapply(parts,`[[`,0,'saturation');K<-vapply(parts,`[[`,0L,'K')
  e<-mesh$edges;inside<-sat<=tau
  key<-if(distinguish_partitions)vapply(parts,ceca_partition_key,'') else as.character(K)
  hit<-xor(inside[e[,1]],inside[e[,2]]) |
    (inside[e[,1]] & inside[e[,2]] & key[e[,1]]!=key[e[,2]])
  list(cells=sort(unique(as.integer(e[hit,,drop=FALSE]))),edges=e[hit,,drop=FALSE])
}

# Endpoints and midpoint of each raster frontier segment, in initial-grid units.
# This is a numerical displacement indicator, not a bound for the unknown true frontier.
ceca_frontier_points<-function(state,cfg,scale,indices=NULL,frontier=NULL) {
  g<-state$mesh$leaves;o<-ceca_adaptive_options(cfg)
  f<-if(is.null(frontier))ceca_frontiers(state$mesh,ceca_current(state),cfg$sat_max2,o$distinguish_partitions) else frontier
  e<-f$edges
  if(!is.null(indices) && nrow(e))e<-e[e[,1]%in%indices & e[,2]%in%indices,,drop=FALSE]
  if(!nrow(e))return(matrix(numeric(),0,2))
  points<-lapply(seq_len(nrow(e)),function(k) {
    a<-g[e[k,1],];b<-g[e[k,2],]
    if(abs(a$u1-b$u0)<1e-10 || abs(b$u1-a$u0)<1e-10) {
      u<-if(abs(a$u1-b$u0)<1e-10)a$u1 else b$u1
      lo<-max(a$v0,b$v0);hi<-min(a$v1,b$v1)
      cbind(rep(u,3),c(lo,(lo+hi)/2,hi))
    } else {
      v<-if(abs(a$v1-b$v0)<1e-10)a$v1 else b$v1
      lo<-max(a$u0,b$u0);hi<-min(a$u1,b$u1)
      cbind(c(lo,(lo+hi)/2,hi),rep(v,3))
    }
  })
  unique(sweep(do.call(rbind,points),2,scale,'/'))
}
ceca_frontier_distance<-function(a,b) {
  if(!nrow(a) && !nrow(b))return(0)
  if(!nrow(a) || !nrow(b))return(Inf)
  directed<-function(x,y)max(vapply(seq_len(nrow(x)),function(i)
    sqrt(min(rowSums(sweep(y,2,x[i,],'-')^2))),0))
  max(directed(a,b),directed(b,a))
}

ceca_refinement_memory<-function(state,previous=NULL) {
  ids<-as.character(state$mesh$leaves$id);parts<-ceca_current(state)
  out<-data.frame(id=state$mesh$leaves$id,quiet=0L,key=vapply(parts,ceca_partition_key,''),
                  saturation=vapply(parts,`[[`,0,'saturation'),stringsAsFactors=FALSE)
  if(!is.null(previous)) {
    j<-match(ids,as.character(previous$id));ok<-!is.na(j)
    out$quiet[ok]<-previous$quiet[j[ok]]
  }
  out
}
ceca_wake_changed<-function(state,memory,tau) {
  now<-ceca_refinement_memory(state,memory);j<-match(now$id,memory$id)
  changed<-which(is.na(j) | now$key!=memory$key[j] |
                   (now$saturation<=tau)!=(memory$saturation[j]<=tau))
  wake<-unique(c(changed,unlist(state$mesh$adj[changed],use.names=FALSE)))
  now$quiet[wake]<-0L
  now
}
ceca_refinement_update<-function(before,after,memory,refined_ids,cfg,scale) {
  o<-ceca_adaptive_options(cfg);out<-ceca_wake_changed(after,memory,cfg$sat_max2)
  oldg<-before$mesh$leaves;newg<-after$mesh$leaves
  oldp<-ceca_current(before);newp<-ceca_current(after);logs<-list()
  oldfront<-ceca_frontiers(before$mesh,oldp,cfg$sat_max2,o$distinguish_partitions)
  newfront<-ceca_frontiers(after$mesh,newp,cfg$sat_max2,o$distinguish_partitions)
  old_allkeys<-vapply(oldp,ceca_partition_key,'');new_allkeys<-vapply(newp,ceca_partition_key,'')
  old_inside<-vapply(oldp,function(p)p$saturation<=cfg$sat_max2,TRUE)
  new_inside<-vapply(newp,function(p)p$saturation<=cfg$sat_max2,TRUE)
  for(id in refined_ids) {
    i<-match(id,oldg$id);local_old<-sort(unique(c(i,before$mesh$adj[[i]])))
    local_ids<-oldg$id[local_old]
    local_new<-which(newg$id%in%local_ids | newg$parent%in%local_ids)
    oldkeys<-old_allkeys[local_old[old_inside[local_old]]]
    newkeys<-new_allkeys[local_new[new_inside[local_new]]]
    novel<-length(setdiff(newkeys,oldkeys))
    shift<-ceca_frontier_distance(ceca_frontier_points(before,cfg,scale,local_old,oldfront),
                               ceca_frontier_points(after,cfg,scale,local_new,newfront))
    quiet<-novel==0L && shift<=o$position_tolerance+1e-12
    count<-if(quiet)memory$quiet[match(id,memory$id)]+1L else 0L
    children<-which(newg$parent==id)
    out$quiet[children]<-count
    logs[[length(logs)+1L]]<-data.frame(parent=id,novel_local_winners=novel,
      displacement_initial_grid_units=shift,quiet=count,sleeping=count>=o$patience)
  }
  list(memory=out,diagnostics=if(length(logs))do.call(rbind,logs) else data.frame())
}
ceca_select_refinement<-function(state,memory,cfg) {
  o<-ceca_adaptive_options(cfg);g<-state$mesh$leaves
  f<-ceca_frontiers(state$mesh,ceca_current(state),cfg$sat_max2,o$distinguish_partitions)
  e<-f$edges
  if(!nrow(e))return(list(indices=integer(),edges=e,status='no_frontiers'))
  quiet<-memory$quiet[match(g$id,memory$id)]
  active<-quiet[e[,1]]<o$patience | quiet[e[,2]]<o$patience
  if(!any(active))return(list(indices=integer(),edges=e[FALSE,,drop=FALSE],status='adaptive_stable'))
  e<-e[active,,drop=FALSE]
  if(any(g$depth[as.integer(e)]>=o$max_depth))
    return(list(indices=integer(),edges=e,status='depth_budget'))
  ids<-sort(unique(as.integer(e)))
  if(nrow(g)+3L*length(ids)>cfg$max_cells)
    return(list(indices=integer(),edges=e,status='cell_budget'))
  list(indices=ids,edges=e,status='refine')
}

# Imported from Moteurs/CEClust_leaf_bic_patch_v075/right_stop.R
# Three actually optimized columns, with exactly matched lambda cells.
# Historical partitions are refreshed after harmonization; no interpolation.
ceca_right_panel <- function(state) {
  g <- ceca_diagram(state); p <- g[ceca_border(state$mesh)$right,]
  p <- p[order(p$lambda),c('id','C','lambda','v0','v1','REO','saturation')]
  p$key <- vapply(as.character(p$id),function(id)ceca_partition_key(state$parts[[id]]),'')
  rownames(p) <- NULL; p
}
ceca_right_memory <- function(state) list(panels=list(ceca_right_panel(state)),resets=0L,
  last_reset='initial',partition_changes=0L)

ceca_right_refresh <- function(state,memory,append=FALSE) {
  if(is.null(memory))return(ceca_right_memory(state))
  current <- ceca_right_panel(state); g <- ceca_diagram(state)
  reset <- function(reason) {
    z <- ceca_right_memory(state);z$resets <- memory$resets+1L;z$last_reset <- reason;z
  }
  for(i in seq_along(memory$panels)) {
    p <- memory$panels[[i]]; at <- match(p$id,g$id)
    if(anyNA(at))return(reset('comparison_cells_refined'))
    p$REO <- g$REO[at];p$saturation <- g$saturation[at]
    keys <- vapply(as.character(p$id),function(id)ceca_partition_key(state$parts[[id]]),'')
    memory$partition_changes <- memory$partition_changes+sum(keys!=p$key)
    p$key <- keys;memory$panels[[i]] <- p
  }
  if(append)memory$panels <- tail(c(memory$panels,list(current)),3L)
  else if(!identical(tail(memory$panels,1)[[1]]$id,current$id))
    return(reset('lambda_domain_or_right_mesh_changed'))
  ref <- current[,c('lambda','v0','v1')]
  if(!all(vapply(memory$panels,function(p)identical(p[,c('lambda','v0','v1')],ref),TRUE)))
    return(reset('lambda_rows_do_not_match'))
  memory
}
ceca_right_status <- function(state,cfg,memory) {
  m <- ceca_margins(state$mesh,ceca_current(state),cfg$sat_max)
  panels <- memory$panels; two <- length(panels)>=3L
  unsaturated <- all(vapply(panels,function(p)all(is.finite(p$saturation)&p$saturation==0),TRUE))
  same <- two && identical(panels[[1]]$REO,panels[[2]]$REO) &&
    identical(panels[[2]]$REO,panels[[3]]$REO)
  last <- tail(panels,1)[[1]]
  old_rows <- last$REO==1L | last$saturation>cfg$sat_max
  stable_rows <- rep(FALSE,nrow(last))
  if(two)stable_rows <- panels[[1]]$saturation==0 & panels[[2]]$saturation==0 & last$saturation==0 &
    panels[[1]]$REO==panels[[2]]$REO & panels[[2]]$REO==last$REO
  row_status <- data.frame(id=last$id,lambda=last$lambda,old_rule=old_rows,
    stable=stable_rows,terminal=old_rows|stable_rows)
  practical <- !isTRUE(m$observed['right']) && all(row_status$terminal)
  m$original_observed <- m$observed
  m$right_stopping <- list(old_rule=isTRUE(m$observed['right']),
    two_extensions=two,all_unsaturated=unsaturated,REO_unchanged=same,
    practical=practical,reason=if(m$observed['right'])'trivial_or_saturated' else
      if(practical)'rowwise_two_extensions_REO_stable' else 'continue',
    available_columns=length(panels),resets=memory$resets,last_reset=memory$last_reset,
    partition_changes=memory$partition_changes,rows=row_status)
  m$observed['right'] <- m$observed['right'] || practical
  m$needs_extension <- !m$observed
  m
}

# Imported from Moteurs/CEClust_leaf_patch_v072/screening.R
ceca_screen_options<-function(cfg) {
  o<-modifyList(list(trivial_quadrants=TRUE,anchor_rule='rsi_and_algo',audit_fraction=.10,
                    stability_min=.95,audit_scheme='shared_persistent'),if(is.null(cfg$screening))list() else cfg$screening)
  stopifnot(o$anchor_rule%in%c('all_trivial','rsi_and_algo'),
    o$audit_fraction>0,o$audit_fraction<=1,o$stability_min>0,o$stability_min<=1,
    identical(o$audit_scheme,'shared_persistent'))
  o
}

ceca_update_quadrant_audits<-function(audit_state,candidate_ids,cfg) {
  o<-ceca_screen_options(cfg)
  if(is.null(audit_state))audit_state<-list(candidate_ids=integer(),audit_ids=integer(),
    seed=cfg$seed,fraction=o$audit_fraction,scheme='shared_persistent')
  stopifnot(identical(audit_state$scheme,'shared_persistent'),
    identical(audit_state$seed,cfg$seed),identical(audit_state$fraction,o$audit_fraction),
    !anyDuplicated(audit_state$candidate_ids),!anyDuplicated(audit_state$audit_ids),
    all(audit_state$audit_ids%in%audit_state$candidate_ids))
  candidates<-sort(unique(c(audit_state$candidate_ids,as.integer(candidate_ids))))
  budget<-as.integer(ceiling(o$audit_fraction*length(candidates)-1e-12))
  missing<-budget-length(audit_state$audit_ids)
  stopifnot(missing>=0L)
  if(missing>0L) {
    available<-setdiff(candidates,audit_state$audit_ids)
    priority<-vapply(available,function(id)ca_seed(cfg$seed,'quadrant_audit_shared',id),0L)
    chosen<-head(available[order(priority,available)],missing)
    audit_state$audit_ids<-sort(c(audit_state$audit_ids,chosen))
  }
  audit_state$candidate_ids<-candidates
  stopifnot(length(audit_state$audit_ids)==budget)
  audit_state
}

# Coordinates are cell representatives. Geometry and the record of measurements
# are kept separate, so a measured cell remains measured inside an admitted region.
ceca_screen_plan<-function(state,cfg,metadata,audit_state=NULL) {
  g<-ceca_diagram(state);o<-ceca_screen_options(cfg)
  result<-data.frame(id=g$id,C=g$C,lambda=g$lambda,REO=g$REO,saturation=g$saturation,
    status=ifelse(g$saturation>cfg$sat_max3,'skipped_saturation','pending'),
    witness=NA_integer_,audit=FALSE,quadrant_covered=FALSE,quadrant_validated=FALSE,
    awaiting_quadrant_audit=FALSE,stringsAsFactors=FALSE)
  measured<-as.character(g$id)%in%names(metadata)
  result$status[measured]<-'measured'
  rejected<-list();regions<-list();candidate_ids<-integer()
  if(o$trivial_quadrants) {
    anchors<-which(measured & g$REO==1L & g$saturation<=cfg$sat_max3)
    anchors<-anchors[order(-g$C[anchors],g$lambda[anchors],g$id[anchors])]
    for(i in anchors) {
      if(!metadata[[as.character(g$id[i])]]$qualifies || metadata[[as.character(g$id[i])]]$best$K!=1L)next
      quad<-which(g$C<=g$C[i]*(1+1e-12) & g$lambda>=g$lambda[i]-1e-12)
      if(any(g$REO[quad]!=1L)) {
        rejected[[length(rejected)+1L]]<-data.frame(anchor=g$id[i],reason='nontrivial_cell_in_quadrant')
        next
      }
      targets<-setdiff(quad[g$saturation[quad]<=cfg$sat_max3],i)
      if(!length(targets))next
      # Keep the existing safeguard: a known failed control invalidates this region.
      known<-targets[measured[targets]]
      if(any(!vapply(known,function(j)metadata[[as.character(g$id[j])]]$qualifies &&
        metadata[[as.character(g$id[j])]]$best$K==1L,TRUE))) {
        rejected[[length(rejected)+1L]]<-data.frame(anchor=g$id[i],reason='measured_control_failed')
        next
      }
      regions[[length(regions)+1L]]<-list(index=i,quad=quad,targets=targets)
      candidate_ids<-c(candidate_ids,g$id[targets])
    }
  }
  audit_state<-ceca_update_quadrant_audits(audit_state,candidate_ids,cfg)
  result$audit<-g$id%in%intersect(audit_state$audit_ids,unique(candidate_ids))
  quadrant_rows<-list()
  for(region in regions) {
    i<-region$index;quad<-region$quad;targets<-region$targets
    result$quadrant_covered[quad]<-TRUE
    audits<-targets[result$audit[targets]]
    pending_audit<-audits[!measured[audits]]
    ready<-!length(pending_audit)
    quadrant_rows[[length(quadrant_rows)+1L]]<-data.frame(anchor=g$id[i],C=g$C[i],lambda=g$lambda[i],
      C_lower=exp(min(g$u0)),lambda_upper=max(g$v1),cells=length(quad),
      eligible_cells=length(targets),shared_controls=length(audits),pending_controls=length(pending_audit),ready=ready)
    if(!ready) {
      result$awaiting_quadrant_audit[targets[!measured[targets] & !result$audit[targets]]]<-TRUE
      result$witness[pending_audit]<-g$id[i]
      next
    }
    result$quadrant_validated[quad]<-TRUE
    omit<-targets[!measured[targets] & !result$audit[targets]]
    result$status[omit]<-'extrapolated_trivial_unmeasured'
    result$witness[omit]<-g$id[i]
  }
  # A selected control cannot be waived by another overlapping quadrant. Pending
  # ordinary cells inside a region awaiting controls do not fill spare worker slots.
  result$status[result$audit & !measured]<-'pending'
  result$awaiting_quadrant_audit[result$status!='pending' | result$audit]<-FALSE
  result$RSI_exhaustif<-NA_real_
  result$Stab_algo_exhaustif<-result$Stab_boot_exhaustif<-NA_real_
  for(i in which(measured)) {
    v<-metadata[[as.character(g$id[i])]]$values
    for(m in 'exhaustif')result[i,paste0('RSI_',m)]<-v$RSI[match(m,v$method)]
    result$Stab_algo_exhaustif[i]<-v$Stab_algo[match('exhaustif',v$method)]
    result$Stab_boot_exhaustif[i]<-v$Stab_boot[match('exhaustif',v$method)]
  }
  pending<-which(result$status=='pending' & !result$awaiting_quadrant_audit)
  # Audits first, then the southeast (Pareto) frontier of the remaining trivial cells.
  audits<-pending[result$audit[pending]]
  trivial<-pending[g$REO[pending]==1L & !result$audit[pending]]
  pareto<-trivial[vapply(trivial,function(i)!any(g$C[trivial]>=g$C[i] & g$lambda[trivial]<=g$lambda[i] &
    (g$C[trivial]>g$C[i] | g$lambda[trivial]<g$lambda[i])),TRUE)]
  other<-setdiff(pending,c(audits,pareto))
  ordered<-c(audits,pareto[order(-g$C[pareto],g$lambda[pareto])],other[order(g$REO[other]==1L,-g$C[other],g$lambda[other])])
  list(table=result,next_indices=head(ordered,cfg$workers),audit_state=audit_state,
    quadrants=if(length(quadrant_rows))do.call(rbind,quadrant_rows) else data.frame(anchor=integer(),C=numeric(),lambda=numeric(),
      C_lower=numeric(),lambda_upper=numeric(),cells=integer(),eligible_cells=integer(),shared_controls=integer(),pending_controls=integer(),ready=logical()),
    rejected=if(length(rejected))do.call(rbind,rejected) else data.frame(anchor=integer(),reason=character()))
}
