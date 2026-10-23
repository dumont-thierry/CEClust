component_fixture <- function() {
  mesh <- ceca_mesh(exp(0:2),c(1,2))
  p <- list(labels=c(1L,1L,2L,2L),K=2L,saturation=0)
  q <- list(labels=c(1L,2L,1L,2L),K=2L,saturation=0)
  parts <- rep(list(p),6);parts[[2]] <- q;parts[[5]] <- q
  list(mesh=mesh,parts=setNames(parts,mesh$leaves$id))
}

test_that('edge connectivity separates equal-order partitions and disconnected regions', {
  s <- component_fixture();p <- CECadaptive_components(s,rep(TRUE,6))
  expect_equal(nrow(p$table),3L)
  expect_equal(p$membership,c(1L,2L,3L,1L,2L,3L))
  # Relabeling is not a different partition.
  s$parts[['4']]$labels <- c(7L,7L,9L,9L)
  expect_equal(CECadaptive_components(s,rep(TRUE,6))$membership,p$membership)
})

test_that('barycentres use full component area, raw log-C and lambda, then saturation gating', {
  s <- component_fixture();s$mesh <- ceca_refine(s$mesh,1L)$mesh
  p <- s$parts[['1']];s$parts[['1']] <- NULL
  for(i in 7:10)s$parts[[as.character(i)]] <- p
  g <- s$mesh$leaves;eligible <- g$id!=4L
  z <- CECadaptive_components(s,eligible);k <- z$membership[match(4L,g$id)]
  expect_equal(z$table$barycentre_logC[k],0)
  expect_equal(z$table$barycentre_lambda[k],1.5)
  expect_false(match(4L,g$id)%in%z$paths[[k]])
  ii <- z$paths[[k]]
  ds <- sqrt((log(g$C[ii])-z$table$barycentre_logC[k])^2+(g$lambda[ii]-1.5)^2)
  expect_false(is.unsorted(ds))
  # Refining one leaf has not moved the component's geometric barycentre.
  expect_equal(z$table$cells[k],5L)
})

test_that('traversal stops on first success and records regional values separately', {
  s <- component_fixture();calls <- integer()
  measure <- function(i) {
    calls <<- c(calls,i)
    pass <- !i%in%c(1L,2L)
    list(diagnostics=list(Stab_algo=if(pass)1 else .8,Stab_boot=1,RSI=1))
  }
  z <- CECadaptive_component_diagnostics(s,rep(TRUE,6),measure,.85,.95)
  expect_equal(calls,c(1L,4L,2L,5L,3L))
  expect_equal(z$plan$table$attempts,c(2L,2L,1L))
  expect_equal(z$witness,c(4L,5L,3L,4L,5L,3L))
  expect_equal(z$origin,c('component_witness','component_witness','measured','measured','measured','component_witness'))
  expect_equal(z$direct$Stab_algo[1],.8);expect_equal(z$values$Stab_algo[1],1)
  expect_true(is.na(z$direct$Stab_algo[6]))
})

test_that('regional CEM freezes the pilot and skips entire independent campaigns', {
  args <- list(Z=CECwith_seed(9,matrix(rnorm(48),24,2)),C_grid=c(1,2),lambda_grid=c(10,20),
    k0=3,B=2,optimization='cem',grid_strategy='adaptive')
  args$optimizer_control <- list(Rmax=1,rinit_max=1,Ninit=1,max_iter=2,max_batches=1)
  args$adaptive_control <- list(max_cells=4,max_rounds=2,max_extensions=0,reference_bootstrap=FALSE)
  a <- do.call(CECfitBoundGrid,args)
  expect_identical(a$adaptive$control$diagnostic_strategy,'components')
  expect_equal(length(a$measurements),1L)
  expect_equal(sum(a$executions$phase=='measurement'),5L)
  expect_equal(sum(a$summary$measured),1L);expect_true(all(a$summary$keep))
  expect_true(all(cec_grid_apply_thresholds(a)$keep))
  expect_true(all(is.na(a$winners[!a$summary$measured,,drop=FALSE])))
  for(i in seq_len(nrow(a$summary)))expect_identical(a$parts[[a$selected[i]]]$labels,
    a$adaptive$pilot$parts[[as.character(a$summary$id[i])]]$labels)
  expect_true(all(is.na(a$summary$rsi_direct[!a$summary$measured])))
  expect_error(CECadaptive_control(list(trivial_quadrants=TRUE),0),'alternative')
})
