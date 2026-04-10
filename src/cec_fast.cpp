#include <RcppArmadillo.h>
#include <limits>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

static arma::mat cec_chol_safe(const arma::mat& sigma) {
  arma::mat Sigma = sigma;
  arma::mat L;
  if (arma::chol(L, Sigma, "lower")) {
    return L;
  }

  const double base_jitter = 1e-10;
  for (int k = 0; k < 8; ++k) {
    double jitter = base_jitter * std::pow(10.0, k);
    arma::mat Sigma_j = Sigma + jitter * arma::eye(Sigma.n_rows, Sigma.n_cols);
    if (arma::chol(L, Sigma_j, "lower")) {
      return L;
    }
  }

  stop("Unable to compute a Cholesky factor for Sigma.");
  return L;
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_choose_clusters(const arma::mat& logdens, const arma::vec& nu) {
  const int n = logdens.n_rows;
  const int r = logdens.n_cols;

  arma::vec lognu(r, arma::fill::zeros);
  for (int k = 0; k < r; ++k) {
    lognu[k] = (nu[k] > 0.0) ? std::log(nu[k]) : -std::numeric_limits<double>::infinity();
  }

  IntegerVector clusters(n);
  NumericVector assigned_logdens(n);
  NumericVector assigned_scores(n);

  for (int i = 0; i < n; ++i) {
    int best_k = 0;
    double best_score = logdens(i, 0) + lognu[0];
    for (int k = 1; k < r; ++k) {
      double score = logdens(i, k) + lognu[k];
      if (score > best_score) {
        best_score = score;
        best_k = k;
      }
    }
    clusters[i] = best_k + 1;
    assigned_logdens[i] = logdens(i, best_k);
    assigned_scores[i] = best_score;
  }

  return List::create(
    _["clusters"] = clusters,
    _["assigned_logdens"] = assigned_logdens,
    _["assigned_scores"] = assigned_scores
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cec_cpp_phi_from_clusters(const IntegerVector& clusters, int r) {
  const int n = clusters.size();
  NumericVector phi(n * r);

  for (int i = 0; i < n; ++i) {
    int cl = clusters[i];
    if (cl >= 1 && cl <= r) {
      phi[i + n * (cl - 1)] = 1.0;
    }
  }

  return phi;
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_gaussian_stats(const arma::mat& X, const IntegerVector& clusters, int r) {
  const int n = X.n_rows;
  const int p = X.n_cols;

  arma::vec counts(r, arma::fill::zeros);
  arma::mat sums(r, p, arma::fill::zeros);

  for (int i = 0; i < n; ++i) {
    int cl = clusters[i] - 1;
    if (cl >= 0 && cl < r) {
      counts[cl] += 1.0;
      sums.row(cl) += X.row(i);
    }
  }

  arma::mat means(r, p, arma::fill::zeros);
  for (int k = 0; k < r; ++k) {
    if (counts[k] > 0.0) {
      means.row(k) = sums.row(k) / counts[k];
    }
  }

  List covs(r);
  for (int k = 0; k < r; ++k) {
    arma::mat cov_k(p, p, arma::fill::zeros);
    if (counts[k] > 0.0) {
      for (int i = 0; i < n; ++i) {
        if (clusters[i] == (k + 1)) {
          arma::rowvec diff = X.row(i) - means.row(k);
          cov_k += diff.t() * diff;
        }
      }
      cov_k /= counts[k];
    }
    covs[k] = cov_k;
  }

  return List::create(
    _["counts"] = counts,
    _["means"] = means,
    _["covs"] = covs
  );
}

// [[Rcpp::export]]
arma::mat cec_cpp_logdens_gaussian(const arma::mat& X, const arma::mat& means, const List& sigma_list, double lambda) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int r = means.n_rows;

  arma::mat out(n, r, arma::fill::zeros);
  const double log2pi = std::log(2.0 * M_PI);
  const double lambda_term = std::log(lambda);

  for (int k = 0; k < r; ++k) {
    arma::mat Sigma = as<arma::mat>(sigma_list[k]);
    arma::mat L = cec_chol_safe(Sigma);
    double logdet = 2.0 * arma::sum(arma::log(L.diag()));
    arma::mat centered = X.each_row() - means.row(k);
    arma::mat solved = arma::solve(arma::trimatl(L), centered.t());
    arma::rowvec quad = arma::sum(arma::square(solved), 0);
    out.col(k) = (-0.5 * (p * log2pi + logdet) - 0.5 * quad.t()) - lambda_term;
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_discrete_counts(const IntegerMatrix& X, const IntegerVector& clusters, int r, int n_levels) {
  const int n = X.nrow();
  const int p = X.ncol();

  std::vector< arma::mat > counts_vec(r, arma::mat(n_levels, p, arma::fill::zeros));
  arma::vec cluster_sizes(r, arma::fill::zeros);

  for (int i = 0; i < n; ++i) {
    int cl = clusters[i] - 1;
    if (cl < 0 || cl >= r) continue;
    cluster_sizes[cl] += 1.0;
    for (int j = 0; j < p; ++j) {
      int code = X(i, j);
      if (code >= 1 && code <= n_levels) {
        counts_vec[cl](code - 1, j) += 1.0;
      }
    }
  }

  List counts(r);
  for (int k = 0; k < r; ++k) {
    counts[k] = counts_vec[k];
  }

  return List::create(
    _["counts"] = counts,
    _["cluster_sizes"] = cluster_sizes
  );
}

// [[Rcpp::export]]
arma::mat cec_cpp_logdens_discrete(const IntegerMatrix& X, const List& prob_list) {
  const int n = X.nrow();
  const int p = X.ncol();
  const int r = prob_list.size();

  arma::mat out(n, r, arma::fill::zeros);

  for (int k = 0; k < r; ++k) {
    arma::mat probs = as<arma::mat>(prob_list[k]);
    for (int i = 0; i < n; ++i) {
      double acc = 0.0;
      for (int j = 0; j < p; ++j) {
        int code = X(i, j);
        if (code < 1 || code > probs.n_rows) {
          acc = -std::numeric_limits<double>::infinity();
          break;
        }
        double pr = probs(code - 1, j);
        if (pr <= 0.0) {
          acc = -std::numeric_limits<double>::infinity();
          break;
        }
        acc += std::log(pr);
      }
      out(i, k) = acc;
    }
  }

  return out;
}

// [[Rcpp::export]]
double cec_cpp_cluster_fingerprint(const IntegerVector& clusters) {
  const int n = clusters.size();
  if (n == 0) return 0.0;

  double acc = 0.0;
  for (int i = 0; i < n; ++i) {
    double angle = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(n);
    acc += static_cast<double>(clusters[i]) * std::cos(angle);
  }
  return acc;
}
