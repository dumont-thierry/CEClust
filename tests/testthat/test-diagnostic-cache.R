test_that('shared cache preserves every exhaustive agreement and its order', {
  a <- c(1L,1L,2L,2L,3L,3L)
  b <- c(1L,2L,2L,2L,3L,3L)
  c <- c(1L,1L,1L,2L,2L,2L)
  cache <- new.env(parent=emptyenv())
  for(metric in c('match','ARI'))for(L in 0:3)for(B in 0:2) {
    if(L+B==0)next
    originals <- list(a,b,c)[seq_len(L)];boots <- list(c,a)[seq_len(B)]
    z <- CECstability_cached(originals,boots,metric,cache)
    expect_identical(z,CECstability(originals,boots,metric))
    aa <- if(L>=2)utils::combn(seq_len(L),2,function(ij)
      partition_similarity(originals[[ij[1]]],originals[[ij[2]]],method=metric)) else numeric()
    bb <- if(L && B)unlist(lapply(originals,function(x)vapply(boots,function(y)
      partition_similarity(x,y,method=metric),0)),use.names=FALSE) else numeric()
    expect_equal(z$algorithm_agreements,aa,tolerance=1e-15)
    expect_equal(z$bootstrap_agreements,bb,tolerance=1e-15)
    expect_equal(z$original_pairs,L*(L-1)/2)
    expect_equal(z$cross_pairs,L*B)
  }
})

test_that('pair keys survive reordered cells and distinguish the metric and sample size', {
  a <- c(1L,1L,2L,2L);b <- c(1L,2L,2L,2L);c <- c(1L,1L,1L,2L)
  cache <- new.env(parent=emptyenv())
  first <- CECstability_cached(list(a,b),list(c,a),'match',cache)
  evaluated <- cache$pair_evaluations
  expect_equal(evaluated,3L)
  again <- CECstability_cached(list(a,b),list(c,a),'match',cache)
  expect_identical(again,first);expect_equal(cache$pattern_hits,1L)
  second <- CECstability_cached(list(b,a),list(a,c),'match',cache)
  expect_identical(second,CECstability(list(b,a),list(a,c)))
  expect_equal(cache$pair_evaluations,evaluated)
  ari <- CECstability_cached(list(a,b),list(c,a),'ARI',cache)
  expect_identical(ari,CECstability(list(a,b),list(c,a),'ARI'))
  expect_equal(cache$pair_evaluations,2*evaluated)
  shorter <- CECstability_cached(list(a[1:3],b[1:3]),list(c[1:3]),'match',cache)
  expect_identical(shorter,CECstability(list(a[1:3],b[1:3]),list(c[1:3])))
})

test_that('cached outputs are not mutated by a caller and no random draw is consumed', {
  a <- c(1L,1L,2L,2L);b <- c(1L,2L,2L,2L)
  cache <- new.env(parent=emptyenv())
  CECwith_seed(174,{
    before <- .Random.seed
    z <- CECstability_cached(list(a,b),list(a),'match',cache)
    expected <- z;z$Stab_algo <- 0;z$algorithm_agreements[1] <- 0
    expect_identical(CECstability_cached(list(a,b),list(a),'match',cache),expected)
    expect_identical(.Random.seed,before)
  })
  key <- digest::digest(a,algo='xxhash64')
  assign(key,b,cache$partitions)
  expect_error(CECstability_cached(list(a,b),list(a),'match',cache),'hash collision')
})

test_that('a diagram with no eligible measurement has a valid empty cache', {
  expect_identical(CECstability_cache_summary(new.env(parent=emptyenv())),
    list(pair_evaluations=0L,pattern_hits=0L,patterns=0L))
  X <- matrix(seq(-1,1,length.out=20),ncol=1)
  for(method in c('cem','em_bank'))for(penalty in c('none','bic')) {
    z <- CECfitBoundGrid(Z=X,C_grid=c(.001,.002),lambda_grid=c(2,3),
      grid_strategy='adaptive',optimization=method,penalty=penalty,
      k0=2,B=2,rsi_policy='unsaturated',
      optimizer_control=list(Rmax=1,rinit_max=1,Ninit=1,max_iter=2,Cmax=2,max_batches=1),
      adaptive_control=list(max_cells=4,max_rounds=2,max_extensions=0))
    expect_length(z$measurements,0L)
    expect_identical(z$meta$agreement_cache,list(pair_evaluations=0L,pattern_hits=0L,patterns=0L))
    expect_false(any(z$summary$keep))
    expect_true(all(z$summary$sat==1))
  }
})
