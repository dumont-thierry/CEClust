# Regional diagnostics are a witness-based approximation, not measurements at
# every cell. Geometry is frozen before this traversal starts.
CECadaptive_components <- function(state,eligible) {
  g <- state$mesh$leaves;p <- ceca_current(state);n <- nrow(g)
  keys <- vapply(p,function(z)digest::digest(ca_labels(z$labels)),"")
  membership <- integer(n);members <- paths <- list();rows <- list()
  for(first in order(g$id)) {
    if(membership[first])next
    k <- length(members)+1L;queue <- first;membership[first] <- k;at <- 1L
    while(at<=length(queue)) {
      i <- queue[at];at <- at+1L
      neighbours <- state$mesh$adj[[i]]
      neighbours <- neighbours[membership[neighbours]==0L & keys[neighbours]==keys[first]]
      for(j in neighbours) {
        if(!identical(ca_labels(p[[j]]$labels),ca_labels(p[[first]]$labels)))
          stop('Partition hash collision in component construction.')
        membership[j] <- k
      }
      queue <- c(queue,neighbours)
    }
    ii <- sort(queue);area <- (g$u1[ii]-g$u0[ii])*(g$v1[ii]-g$v0[ii])
    u <- sum(area*(g$u0[ii]+g$u1[ii])/2)/sum(area)
    v <- sum(area*(g$v0[ii]+g$v1[ii])/2)/sum(area)
    distance <- sqrt((log(g$C[ii])-u)^2+(g$lambda[ii]-v)^2)
    ordered <- order(distance,g$id[ii]);candidates <- ii[ordered][eligible[ii[ordered]]]
    members[[k]] <- ii;paths[[k]] <- candidates
    rows[[k]] <- data.frame(component=k,first_cell_id=min(g$id[ii]),cells=length(ii),
      candidates=length(candidates),REO=p[[first]]$K,barycentre_logC=u,barycentre_lambda=v,
      status=if(length(candidates))'pending' else 'outside_measurement_policy',
      attempts=0L,witness_cell_id=NA_integer_)
  }
  list(membership=membership,members=members,paths=paths,table=do.call(rbind,rows))
}

CECadaptive_component_diagnostics <- function(state,eligible,measure,algo_min,rsi_min,verbose=FALSE,
  require_reference_agreement=FALSE) {
  plan <- CECadaptive_components(state,eligible);g <- state$mesh$leaves;M <- nrow(g)
  direct <- data.frame(Stab_algo=rep(NA_real_,M),Stab_boot=NA_real_,RSI=NA_real_)
  agreement <- rep(NA_real_,M);reference_pass <- rep(!require_reference_agreement,M)
  measurements <- list();trials <- list()
  for(k in seq_along(plan$paths)) {
    path <- plan$paths[[k]]
    for(rank in seq_along(path)) {
      i <- path[rank];z <- measure(i);d <- z$diagnostics
      direct[i,] <- c(d$Stab_algo,d$Stab_boot,d$RSI)
      measurements[[as.character(g$id[i])]] <- z
      diagnostic_pass <- is.finite(d$Stab_algo) && d$Stab_algo>algo_min && is.finite(d$RSI) && d$RSI>rsi_min
      a <- z$reference_check$agreement
      if(is.numeric(a) && length(a)==1L)agreement[i] <- a
      reference_pass[i] <- !require_reference_agreement || (is.finite(agreement[i]) && agreement[i]>=algo_min)
      pass <- diagnostic_pass && reference_pass[i]
      plan$table$attempts[k] <- rank
      trials[[length(trials)+1L]] <- data.frame(component=k,rank=rank,cell_id=g$id[i],
        distance_barycentre=sqrt((log(g$C[i])-plan$table$barycentre_logC[k])^2+
          (g$lambda[i]-plan$table$barycentre_lambda[k])^2),accepted=pass,
        diagnostics_passed=diagnostic_pass,reference_agreement=agreement[i],reference_passed=reference_pass[i],
        rejection_reason=if(pass)'ok' else if(!diagnostic_pass)'stability_or_rsi' else 'reference_disagreement')
      if(isTRUE(verbose))message('Component ',k,', cell ',g$id[i],': ',if(pass)'accepted' else 'next candidate')
      if(pass) {
        plan$table$status[k] <- 'accepted';plan$table$witness_cell_id[k] <- g$id[i];break
      }
    }
    if(length(path) && plan$table$status[k]=='pending')plan$table$status[k] <- 'exhausted'
  }
  values <- direct;witness <- rep(NA_integer_,M)
  origin <- ifelse(as.character(g$id)%in%names(measurements),'measured','not_measured')
  for(k in which(plan$table$status=='accepted')) {
    ii <- plan$members[[k]];id <- plan$table$witness_cell_id[k];j <- match(id,g$id)
    values[ii,] <- direct[rep(j,length(ii)),,drop=FALSE]
    agreement[ii] <- agreement[j];reference_pass[ii] <- reference_pass[j]
    witness[ii] <- id;origin[ii] <- ifelse(g$id[ii]==id,'measured','component_witness')
  }
  list(values=values,direct=direct,measurements=measurements,witness=witness,origin=origin,
    plan=plan,trials=if(length(trials))do.call(rbind,trials) else data.frame(),
    require_reference_agreement=require_reference_agreement,reference_agreement=agreement,reference_pass=reference_pass)
}

CECadaptive_component_summary <- function(g,regional) {
  g$stab_algo <- regional$values$Stab_algo;g$stab_boot <- regional$values$Stab_boot;g$rsi <- regional$values$RSI
  g$stab_algo_direct <- regional$direct$Stab_algo;g$stab_boot_direct <- regional$direct$Stab_boot
  g$rsi_direct <- regional$direct$RSI
  g$measured <- as.character(g$id)%in%names(regional$measurements)
  g$diagnostic_origin <- regional$origin;g$diagnostic_witness <- regional$witness
  g$component <- regional$plan$membership
  g$component_status <- regional$plan$table$status[g$component]
  g$diagnostics_available <- is.finite(g$stab_algo)&is.finite(g$stab_boot)
  g$measurement_status <- ifelse(g$diagnostic_origin=='component_witness','component_witness',
    ifelse(g$measured,'measured',ifelse(g$component_status=='outside_measurement_policy',
      'outside_measurement_policy','component_not_tested')))
  # The optional independent-reference gate participates in traversal and in
  # selection, including directly measured points of exhausted components.
  g$reference_agreement_required <- isTRUE(regional$require_reference_agreement)
  g$selection_reference_agreement <- regional$reference_agreement
  g$keep_harmonization <- regional$reference_pass
  g
}
