# Registered native functions; regeneration is unnecessary for this fixed API.
bands_state_cpp <- function(X, labels, C, lambda) {
  .Call(`_CEClust_bands_state_cpp`, X, labels, C, lambda)
}

bands_cem_cpp <- function(X, labels, C, lambda, max_iter = 1000L,
                          tolerance = 1e-8, beta = 0) {
  .Call(`_CEClust_bands_cem_cpp`, X, labels, C, lambda,
        max_iter, tolerance, beta)
}
