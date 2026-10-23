test_that('fused Gaussian evaluation preserves determinant profiles and raw densities', {
  set.seed(63)
  fixtures <- list()
  for(d in c(1,2,6))for(K in c(1,2,8))for(shift in c(0,1e8))
    fixtures[[length(fixtures)+1L]] <- list(X=matrix(rnorm(160*d),160,d)+shift,z=sample(rep(1:K,length.out=160)))
  for(X in list(matrix(0,20,2),cbind(1:20,2*(1:20))))
    for(z in list(rep(1L,20),rep(1:10,each=2),1:20))fixtures[[length(fixtures)+1L]] <- list(X=X,z=z)
  for(f in fixtures)for(C in c(.02,1,20)) {
    X <- f$X;z <- ca_labels(f$z)
    expected <- CECsolver_environment('none')$q_part(X,z,C,.93)
    a <- cec_cpp_cem_evaluate(X,z,C,.93)
    p <- CECsolver_environment('none')$q_part(X,z,C,.93,stats=a$stats)
    expect_equal(p$J,expected$J,tolerance=1e-10)
    expect_identical(p$saturated,expected$saturated)
    backend <- CECprepare_backend_data(as.data.frame(X),familyType='gaussVector')
    param <- CECoptParam_fast_from_clusters(backend,z,.93,C,Inf)
    density <- CECcompute_logdens_matrix_fast(backend,param,.93)
    expect_identical(is.finite(a$logdens),is.finite(density))
    finite <- is.finite(density)
    expect_equal(a$logdens[finite],density[finite],tolerance=1e-7)
  }
})

test_that('compiled and reference CEM keep starts, accepted moves and winners', {
  set.seed(4)
  for(d in c(1,2,6))for(penalty in c('none','bic'))for(lambda in c(.7,1.4)) {
    X <- matrix(rnorm(64*d),64,d);X[1:32,1] <- X[1:32,1]+3
    settings <- q_settings();settings$rinit <- 4L;settings$shots <- 2L
    settings$rinit_max <- 12L;settings$max_batches <- 2L;settings$Nloop <- 40L
    solver <- CECsolver_environment(penalty)
    settings$cem_engine <- 'reference'
    old <- solver$q_fit(X,2,lambda,18L,settings)
    settings$cem_engine <- 'compiled'
    new <- solver$q_fit(X,2,lambda,18L,settings)
    expect_identical(new$part$labels,old$part$labels)
    expect_equal(new$part$J,old$part$J,tolerance=1e-10)
    expect_identical(new$status,old$status)
    expect_identical(new$starts$initial_K,old$starts$initial_K)
    expect_identical(new$starts$final_K,old$starts$final_K)
    expect_identical(lapply(new$descent_traces,function(x)x[,c('move','K')]),
      lapply(old$descent_traces,function(x)x[,c('move','K')]))
    for(i in seq_along(new$descent_traces))expect_equal(new$descent_traces[[i]]$J,
      old$descent_traces[[i]]$J,tolerance=1e-10)
    expect_true(all(new$starts$engine%in%c('compiled','reference_fallback')))
  }
})

test_that('a compiled numerical error uses the reference descent without dropping the start', {
  X <- matrix(seq_len(24),ncol=2);settings <- q_settings()
  settings$Nloop <- 4L
  backend <- CECprepare_backend_data(as.data.frame(X),familyType='gaussVector')
  solver <- CECsolver_environment('bic')
  old <- solver$q_descent_R(X,2,1,rep(1:2,each=6),settings,backend)
  expect_identical(names(old$stage_timings),c('refit','CEM_proposal','deletion'))
  expect_true(all(is.finite(old$stage_timings)))
  local_mocked_bindings(cec_cpp_cem_descent=function(...)stop('Forced numerical failure'),.package='CEClust')
  solver <- CECsolver_environment('bic')
  new <- solver$q_descent(X,2,1,rep(1:2,each=6),settings,backend)
  expect_identical(new$part,old$part);expect_identical(new$descent_trace,old$descent_trace)
  expect_identical(new$engine,'reference_fallback');expect_match(new$compiled_fallback,'Forced')
})

test_that('memoized exhaustive stability preserves every agreement and its order', {
  set.seed(27)
  labels <- c(rep(list(rep(1L,100)),5),rep(list(rep(1:2,each=50)),7),
    replicate(8,sample(1:3,100,replace=TRUE),simplify=FALSE))
  original <- labels[1:12];bootstrap <- labels[13:20]
  for(metric in c('match','ARI')) {
    a <- utils::combn(seq_along(original),2,function(ij)
      partition_similarity(original[[ij[1]]],original[[ij[2]]],method=metric))
    b <- unlist(lapply(original,function(x)vapply(bootstrap,function(y)partition_similarity(x,y,method=metric),0)),use.names=FALSE)
    got <- CECstability(original,bootstrap,metric)
    expect_identical(got$algorithm_agreements,a);expect_identical(got$bootstrap_agreements,b)
    expect_identical(got$Stab_algo,mean(a));expect_identical(got$Stab_boot,mean(b))
  }
})
