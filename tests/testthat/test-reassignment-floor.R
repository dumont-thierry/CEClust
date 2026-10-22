test_that("the proposal floor is positive below C for all finite positive C", {
  floor_log <- getFromNamespace("CECreassignment_log_floor", "CEClust")
  for (C in c(1e-100, .01, .2, 1, 10, 100, 1e100)) {
    b <- exp(floor_log(C))
    expect_gt(b, 0)
    expect_lt(b, C)
    expect_equal(b, min(1/C, C/1000), tolerance=1e-12)
  }
  expect_equal(floor_log(100, 100), floor_log(100, 1000))
})

test_that("singular and regular components share the same proposal floor", {
  raw_density <- getFromNamespace("CECcompute_logdens_matrix_fast", "CEClust")
  fit_cov <- getFromNamespace("CECfit_covariance_determinant", "CEClust")
  proposal <- getFromNamespace("CECproposal_logdens", "CEClust")
  toggle <- getFromNamespace("CECset_fast_backend", "CEClust")
  on.exit(toggle(TRUE))
  for (C in c(.2, 1, 100)) {
    b <- rep(1/(sqrt(2*pi)*sqrt(C)), 2)
    params <- list(states=1:2, nu=c(.8,.2), m=matrix(0,2,2), C=C, lambda=1,
      Sigma=list(fit_cov(diag(c(1,0)),b)$Sigma_corrected,
                 fit_cov(diag(1e6,2),b)$Sigma_corrected), familyType="gaussVector")
    X <- matrix(c(0,1e6),nrow=1)
    backend <- CECprepare_backend_data(X, "gaussVector")
    for (enabled in c(FALSE,TRUE)) {
      toggle(enabled)
      raw <- raw_density(backend, params)
      floored <- proposal(raw, params)
      expect_equal(as.numeric(exp(floored)), rep(min(1/C,C/1000),2))
      posterior <- sweep(exp(floored),2,params$nu,"*")
      posterior <- posterior/rowSums(posterior)
      expect_equal(as.numeric(posterior), params$nu)
      step <- getFromNamespace("CECoptPhi_backend", "CEClust")(
        X,params,backend_data=backend)
      expect_equal(step$clusters,1L)
      expect_equal(as.numeric(step$assigned_logdens), -Inf)
      expect_equal(raw_density(backend,params),raw)
    }
  }
})

test_that("a floor-induced uphill merge is rejected using original entropy", {
  X <- rbind(1000*rbind(c(-1,-1),c(-1,1),c(1,-1),c(1,1)),
             sweep(1000*rbind(c(-1,-1),c(1,-1),c(0,1)),2,c(1e5,0),"+"))
  labels <- c(rep(1L,4),rep(2L,3))
  make_phi <- getFromNamespace("CECclusters_to_phi", "CEClust")
  phi <- make_phi(labels,2)
  entropy <- getFromNamespace("evalCompositeEntropy", "CEClust")
  initial <- entropy(phi,X,lambda=1,C=1,familyType="gaussVector")
  merged <- entropy(make_phi(rep(1L,7),1),X,lambda=1,C=1,familyType="gaussVector")
  expect_gt(merged,initial)
  backend <- CECprepare_backend_data(X,"gaussVector")
  for (runtime in c("base","fast")) {
    CECconfigure_runtime(runtime)
    for (prepared in list(backend,list(optimized=FALSE,raw=X))) {
      result <- getFromNamespace("CECclassifOneShot", "CEClust")(
        X,lambda=1,C=1,Nloop=5,phi0=phi,familyType="gaussVector",backend_data=prepared)
      got_labels <- max.col(matrix(result$phi,nrow(X)),ties.method="first")
      expect_equal(got_labels,labels)
      expect_equal(result$Hphi,initial,tolerance=1e-10)
      expect_equal(result$descent_diagnostics$rejected_steps,1L)
      expect_equal(result$descent_diagnostics$stop_reason,"refitted_objective_increased")
      expect_equal(result$descent_diagnostics$n_assigned_below_floor,nrow(X))
      expect_true(all(diff(result$descent_diagnostics$H_trace) <= 0))
      clipped_cost <- -sum(result$params$nu*log(result$params$nu))-log(.001)
      expect_gt(result$Hphi,clipped_cost)
    }
  }
})

test_that("a regular single-cluster entropy is unchanged by the score floor", {
  X <- rbind(c(-1,0),c(1,0),c(0,-1),c(0,1),c(1e4,1e4))
  S <- cov(X)*(nrow(X)-1)/nrow(X)
  expected <- log(2*pi*exp(1)) + sum(log(diag(chol(S))))
  result <- CECclassif(X,lambda=1,C=1,r0=1,Nshots=1,Nloop=3,familyType="gaussVector")
  expect_equal(result$Hphi,expected,tolerance=1e-10)
  expect_gt(result$descent_diagnostics$n_assigned_below_floor,0)
  expect_equal(result$descent_diagnostics$acceptance_criterion,
                "original_entropy_without_lower_floor")
})
