test_that('left extension requires saturation at every lambda, including trivial cells', {
  mesh <- ceca_mesh(c(1,2),c(.5,1,1.5))
  border <- ceca_border(mesh)
  parts <- replicate(nrow(mesh$leaves),list(K=1L,saturation=0),simplify=FALSE)
  before <- ceca_margins(mesh,parts,0)
  expect_false(unname(before$observed['left']))
  expect_true(unname(before$needs_extension['left']))
  expect_true(unname(before$observed['right']))
  expect_true(unname(before$observed['upper']))
  expect_false(unname(before$observed['lower']))
  for(i in border$left)parts[[i]]$saturation <- .1
  expect_true(unname(ceca_margins(mesh,parts,0)$observed['left']))
  expect_false(unname(ceca_margins(mesh,parts,.1)$observed['left']))
  parts[[border$left[2]]] <- list(K=3L,saturation=0)
  expect_false(unname(ceca_margins(mesh,parts,0)$observed['left']))
  parts[[border$left[2]]]$saturation <- .001
  expect_true(unname(ceca_margins(mesh,parts,0)$observed['left']))
  refined <- ceca_refine(mesh,border$left)$mesh
  parts <- replicate(nrow(refined$leaves),list(K=1L,saturation=1),simplify=FALSE)
  left <- ceca_border(refined)$left
  expect_true(unname(ceca_margins(refined,parts,0)$observed['left']))
  parts[[tail(left,1)]]$saturation <- 0
  expect_false(unname(ceca_margins(refined,parts,0)$observed['left']))
})

test_that('both adaptive optimizers extend a trivial unsaturated left edge', {
  X <- matrix(seq(-2,2,length.out=24),ncol=1)
  for(engine in c('cem','em_bank'))for(penalty in c('none','bic')) {
    args <- list(Z=X,C_grid=c(2,4),lambda_grid=c(2,3),k0=2,B=1,
      algo_seed=13,optimization=engine,penalty=penalty,grid_strategy='adaptive',
      optimizer_control=list(Rmax=1,Ninit=1,max_iter=5,max_batches=1,Cmax=16),
      adaptive_control=list(max_cells=100,max_rounds=2,max_extensions=6,
        diagnostic_strategy='components',harmonize=FALSE),checkpoint_dir=FALSE)
    a <- do.call(CECfitBoundGrid,args)
    h <- a$adaptive$history
    expect_false(h$left[1])
    expect_gt(unname(a$adaptive$extensions['left']),0)
    expect_true(any(h$stage=='extension_left'))
    expect_true(unname(a$adaptive$margins$observed['left']))
    left <- ceca_border(a$adaptive$mesh)$left
    expect_true(all(a$summary$saturation[left]>0))
  }
})
