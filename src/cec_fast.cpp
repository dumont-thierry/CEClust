#include <RcppArmadillo.h>
#include <R_ext/Random.h>
#include <limits>
#include <cmath>
#include <vector>

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
IntegerVector cec_cpp_sample_int_replace(int r, int n) {
  IntegerVector out(n);
  int* out_ptr = out.begin();
  const double r_double = static_cast<double>(r);

  if (r == 1) {
    std::fill(out.begin(), out.end(), 1);
    for (int i = 0; i < n; ++i) {
      R_unif_index(r_double);
    }
    return out;
  }

  for (int i = 0; i < n; ++i) {
    out_ptr[i] = static_cast<int>(R_unif_index(r_double)) + 1;
  }

  return out;
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
Rcpp::List cec_cpp_choose_clusters_from_logdens(const arma::mat& logdens, const arma::vec& nu, double lambda) {
  const int n = logdens.n_rows;
  const int r = logdens.n_cols;

  arma::vec lognu(r, arma::fill::zeros);
  for (int k = 0; k < r; ++k) {
    lognu[k] = (nu[k] > 0.0) ? std::log(nu[k]) : -std::numeric_limits<double>::infinity();
  }

  IntegerVector clusters(n);
  NumericVector assigned_logdens(n);
  NumericVector assigned_scores(n);
  const double inv_lambda = 1.0 / lambda;

  for (int i = 0; i < n; ++i) {
    int best_k = 0;
    double best_score = logdens(i, 0) * inv_lambda + lognu[0];
    for (int k = 1; k < r; ++k) {
      double score = logdens(i, k) * inv_lambda + lognu[k];
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
IntegerVector cec_cpp_gauss_univ_assign(
    const NumericVector& z,
    const NumericVector& nu,
    const NumericVector& m,
    const NumericVector& s,
    double lambda) {
  const int n = z.size();
  const int r = nu.size();
  IntegerVector clusters(n);
  int* clusters_ptr = clusters.begin();

  if (r == 1) {
    std::fill(clusters.begin(), clusters.end(), 1);
    return clusters;
  }

  const double* z_ptr = z.begin();
  const double* nu_ptr = nu.begin();
  const double* m_ptr = m.begin();
  const double* s_ptr = s.begin();
  const double inv_lambda = 1.0 / lambda;
  std::vector<double> score_const(r);
  std::vector<double> score_quad(r);
  for (int k = 0; k < r; ++k) {
    score_const[k] = std::log(nu_ptr[k]) - std::log(s_ptr[k]) * inv_lambda;
    score_quad[k] = -0.5 * inv_lambda / (s_ptr[k] * s_ptr[k]);
  }

  for (int i = 0; i < n; ++i) {
    const double zi = z_ptr[i];
    int best_k = 0;
    double diff = zi - m_ptr[0];
    double best_score = score_const[0] + score_quad[0] * diff * diff;

    for (int k = 1; k < r; ++k) {
      diff = zi - m_ptr[k];
      const double score = score_const[k] + score_quad[k] * diff * diff;
      if (score > best_score) {
        best_score = score;
        best_k = k;
      }
    }

    clusters_ptr[i] = best_k + 1;
  }

  return clusters;
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_gauss_univ_cluster_sums(
    const NumericVector& z,
    const IntegerVector& clusters,
    int r) {
  const int n = z.size();
  NumericVector counts(r);
  NumericVector sum_x(r);
  NumericVector sum_x2(r);

  const double* z_ptr = z.begin();
  const int* clusters_ptr = clusters.begin();
  double* counts_ptr = counts.begin();
  double* sum_x_ptr = sum_x.begin();
  double* sum_x2_ptr = sum_x2.begin();

  for (int i = 0; i < n; ++i) {
    const int k = clusters_ptr[i] - 1;
    if (k >= 0 && k < r) {
      const double zi = z_ptr[i];
      counts_ptr[k] += 1.0;
      sum_x_ptr[k] += zi;
      sum_x2_ptr[k] += zi * zi;
    }
  }

  return Rcpp::List::create(
    Rcpp::_["counts"] = counts,
    Rcpp::_["sum_x"] = sum_x,
    Rcpp::_["sum_x2"] = sum_x2
  );
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_gauss_univ_cluster_sums_known_counts(
    const NumericVector& z,
    const IntegerVector& clusters,
    int r) {
  const int n = z.size();
  NumericVector sum_x(r);
  NumericVector sum_x2(r);

  const double* z_ptr = z.begin();
  const int* clusters_ptr = clusters.begin();
  double* sum_x_ptr = sum_x.begin();
  double* sum_x2_ptr = sum_x2.begin();

  for (int i = 0; i < n; ++i) {
    const int k = clusters_ptr[i] - 1;
    if (k >= 0 && k < r) {
      const double zi = z_ptr[i];
      sum_x_ptr[k] += zi;
      sum_x2_ptr[k] += zi * zi;
    }
  }

  return Rcpp::List::create(
    Rcpp::_["sum_x"] = sum_x,
    Rcpp::_["sum_x2"] = sum_x2
  );
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_gauss_univ_params_known_counts(
    const NumericVector& z,
    const IntegerVector& clusters,
    const NumericVector& counts,
    int r,
    double lambda,
    double C) {
  const int n = z.size();
  NumericVector sum_x(r);
  NumericVector sum_x2(r);

  const double* z_ptr = z.begin();
  const int* clusters_ptr = clusters.begin();
  double* sum_x_ptr = sum_x.begin();
  double* sum_x2_ptr = sum_x2.begin();

  for (int i = 0; i < n; ++i) {
    const int k = clusters_ptr[i] - 1;
    if (k >= 0 && k < r) {
      const double zi = z_ptr[i];
      sum_x_ptr[k] += zi;
      sum_x2_ptr[k] += zi * zi;
    }
  }

  IntegerVector states(r);
  NumericVector nu(r);
  NumericVector m(r);
  NumericVector s(r);
  NumericVector varPhi(r);
  std::vector<int> bounded;
  bounded.reserve(r);

  const double* counts_ptr = counts.begin();
  const double s_min = 1.0 / (std::sqrt(2.0 * M_PI) * C);

  for (int k = 0; k < r; ++k) {
    states[k] = k + 1;
    const double count = counts_ptr[k];
    nu[k] = count / static_cast<double>(n);
    m[k] = sum_x_ptr[k] / count;
    double var = sum_x2_ptr[k] / count - m[k] * m[k];
    if (var < 0.0) {
      var = 0.0;
    }
    varPhi[k] = var;
    double sk = std::sqrt(var);
    if (sk < s_min) {
      sk = s_min;
      bounded.push_back(k + 1);
    }
    s[k] = sk;
  }

  IntegerVector functionBoundReached(bounded.begin(), bounded.end());

  return Rcpp::List::create(
    Rcpp::_["states"] = states,
    Rcpp::_["nu"] = nu,
    Rcpp::_["m"] = m,
    Rcpp::_["s"] = s,
    Rcpp::_["varPhi"] = varPhi,
    Rcpp::_["lambda"] = lambda,
    Rcpp::_["C"] = C,
    Rcpp::_["familyType"] = "gaussUniv",
    Rcpp::_["functionBoundReached"] = functionBoundReached
  );
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_compress_integer_clusters_info(
    const IntegerVector& clusters,
    int r) {
  const int n = clusters.size();
  const int* clusters_ptr = clusters.begin();
  IntegerVector counts_int(r);
  int* counts_int_ptr = counts_int.begin();

  for (int i = 0; i < n; ++i) {
    const int cl = clusters_ptr[i];
    if (cl >= 1 && cl <= r) {
      counts_int_ptr[cl - 1] += 1;
    }
  }

  bool all_present = true;
  int r_out = 0;
  for (int k = 0; k < r; ++k) {
    if (counts_int_ptr[k] > 0) {
      ++r_out;
    } else {
      all_present = false;
    }
  }

  if (all_present) {
    NumericVector counts(r);
    double* counts_ptr = counts.begin();
    for (int k = 0; k < r; ++k) {
      counts_ptr[k] = counts_int_ptr[k];
    }
    return Rcpp::List::create(
      Rcpp::_["clusters"] = clusters,
      Rcpp::_["r"] = r,
      Rcpp::_["counts"] = counts
    );
  }

  IntegerVector map(r);
  int* map_ptr = map.begin();
  NumericVector counts(r_out);
  double* counts_ptr = counts.begin();
  int next = 0;
  for (int k = 0; k < r; ++k) {
    if (counts_int_ptr[k] > 0) {
      map_ptr[k] = next + 1;
      counts_ptr[next] = counts_int_ptr[k];
      ++next;
    }
  }

  IntegerVector compressed(n);
  int* compressed_ptr = compressed.begin();
  for (int i = 0; i < n; ++i) {
    const int cl = clusters_ptr[i];
    if (cl >= 1 && cl <= r) {
      compressed_ptr[i] = map_ptr[cl - 1];
    } else {
      compressed_ptr[i] = NA_INTEGER;
    }
  }

  return Rcpp::List::create(
    Rcpp::_["clusters"] = compressed,
    Rcpp::_["r"] = r_out,
    Rcpp::_["counts"] = counts
  );
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_compress_integer_clusters_info_fingerprint(
    const IntegerVector& clusters,
    int r,
    const NumericVector& weights) {
  const int n = clusters.size();
  const int* clusters_ptr = clusters.begin();
  const double* weights_ptr = weights.begin();
  IntegerVector counts_int(r);
  int* counts_int_ptr = counts_int.begin();

  for (int i = 0; i < n; ++i) {
    const int cl = clusters_ptr[i];
    if (cl >= 1 && cl <= r) {
      counts_int_ptr[cl - 1] += 1;
    }
  }

  bool all_present = true;
  int r_out = 0;
  for (int k = 0; k < r; ++k) {
    if (counts_int_ptr[k] > 0) {
      ++r_out;
    } else {
      all_present = false;
    }
  }

  if (all_present) {
    NumericVector counts(r);
    double* counts_ptr = counts.begin();
    for (int k = 0; k < r; ++k) {
      counts_ptr[k] = counts_int_ptr[k];
    }
    double fingerprint = 0.0;
    for (int i = 0; i < n; ++i) {
      fingerprint += static_cast<double>(clusters_ptr[i]) * weights_ptr[i];
    }
    return Rcpp::List::create(
      Rcpp::_["clusters"] = clusters,
      Rcpp::_["r"] = r,
      Rcpp::_["counts"] = counts,
      Rcpp::_["fingerprint"] = fingerprint
    );
  }

  IntegerVector map(r);
  int* map_ptr = map.begin();
  NumericVector counts(r_out);
  double* counts_ptr = counts.begin();
  int next = 0;
  for (int k = 0; k < r; ++k) {
    if (counts_int_ptr[k] > 0) {
      map_ptr[k] = next + 1;
      counts_ptr[next] = counts_int_ptr[k];
      ++next;
    }
  }

  IntegerVector compressed(n);
  int* compressed_ptr = compressed.begin();
  double fingerprint = 0.0;
  for (int i = 0; i < n; ++i) {
    const int cl = clusters_ptr[i];
    if (cl >= 1 && cl <= r) {
      const int compressed_cl = map_ptr[cl - 1];
      compressed_ptr[i] = compressed_cl;
      fingerprint += static_cast<double>(compressed_cl) * weights_ptr[i];
    } else {
      compressed_ptr[i] = NA_INTEGER;
    }
  }

  return Rcpp::List::create(
    Rcpp::_["clusters"] = compressed,
    Rcpp::_["r"] = r_out,
    Rcpp::_["counts"] = counts,
    Rcpp::_["fingerprint"] = fingerprint
  );
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_gauss_univ_compress_params_fingerprint(
    const NumericVector& z,
    const IntegerVector& clusters,
    int r,
    const NumericVector& weights,
    double lambda,
    double C) {
  const int n = clusters.size();
  const double* z_ptr = z.begin();
  const int* clusters_ptr = clusters.begin();
  const double* weights_ptr = weights.begin();
  IntegerVector counts_int(r);
  int* counts_int_ptr = counts_int.begin();
  NumericVector sum_x_raw(r);
  NumericVector sum_x2_raw(r);
  double* sum_x_raw_ptr = sum_x_raw.begin();
  double* sum_x2_raw_ptr = sum_x2_raw.begin();
  double fingerprint_raw = 0.0;

  for (int i = 0; i < n; ++i) {
    const int cl = clusters_ptr[i];
    if (cl >= 1 && cl <= r) {
      const int k = cl - 1;
      const double zi = z_ptr[i];
      counts_int_ptr[k] += 1;
      sum_x_raw_ptr[k] += zi;
      sum_x2_raw_ptr[k] += zi * zi;
      fingerprint_raw += static_cast<double>(cl) * weights_ptr[i];
    }
  }

  bool all_present = true;
  int r_out = 0;
  for (int k = 0; k < r; ++k) {
    if (counts_int_ptr[k] > 0) {
      ++r_out;
    } else {
      all_present = false;
    }
  }

  NumericVector counts(r_out);
  NumericVector sum_x(r_out);
  NumericVector sum_x2(r_out);
  double* counts_ptr = counts.begin();
  double* sum_x_ptr = sum_x.begin();
  double* sum_x2_ptr = sum_x2.begin();
  IntegerVector out_clusters = clusters;
  double fingerprint = fingerprint_raw;

  if (all_present) {
    for (int k = 0; k < r; ++k) {
      counts_ptr[k] = counts_int_ptr[k];
      sum_x_ptr[k] = sum_x_raw_ptr[k];
      sum_x2_ptr[k] = sum_x2_raw_ptr[k];
    }
  } else {
    IntegerVector map(r);
    int* map_ptr = map.begin();
    int next = 0;
    for (int k = 0; k < r; ++k) {
      if (counts_int_ptr[k] > 0) {
        map_ptr[k] = next + 1;
        counts_ptr[next] = counts_int_ptr[k];
        ++next;
      }
    }

    out_clusters = IntegerVector(n);
    int* out_clusters_ptr = out_clusters.begin();
    fingerprint = 0.0;
    for (int i = 0; i < n; ++i) {
      const int cl = clusters_ptr[i];
      if (cl >= 1 && cl <= r) {
        const int compressed_cl = map_ptr[cl - 1];
        const int k = compressed_cl - 1;
        const double zi = z_ptr[i];
        out_clusters_ptr[i] = compressed_cl;
        sum_x_ptr[k] += zi;
        sum_x2_ptr[k] += zi * zi;
        fingerprint += static_cast<double>(compressed_cl) * weights_ptr[i];
      } else {
        out_clusters_ptr[i] = NA_INTEGER;
      }
    }
  }

  IntegerVector states(r_out);
  NumericVector nu(r_out);
  NumericVector m(r_out);
  NumericVector s(r_out);
  NumericVector varPhi(r_out);
  std::vector<int> bounded;
  bounded.reserve(r_out);
  const double s_min = 1.0 / (std::sqrt(2.0 * M_PI) * C);

  for (int k = 0; k < r_out; ++k) {
    states[k] = k + 1;
    const double count = counts_ptr[k];
    nu[k] = count / static_cast<double>(n);
    m[k] = sum_x_ptr[k] / count;
    double var = sum_x2_ptr[k] / count - m[k] * m[k];
    if (var < 0.0) {
      var = 0.0;
    }
    varPhi[k] = var;
    double sk = std::sqrt(var);
    if (sk < s_min) {
      sk = s_min;
      bounded.push_back(k + 1);
    }
    s[k] = sk;
  }

  IntegerVector functionBoundReached(bounded.begin(), bounded.end());
  Rcpp::List params = Rcpp::List::create(
    Rcpp::_["states"] = states,
    Rcpp::_["nu"] = nu,
    Rcpp::_["m"] = m,
    Rcpp::_["s"] = s,
    Rcpp::_["varPhi"] = varPhi,
    Rcpp::_["lambda"] = lambda,
    Rcpp::_["C"] = C,
    Rcpp::_["familyType"] = "gaussUniv",
    Rcpp::_["functionBoundReached"] = functionBoundReached
  );

  return Rcpp::List::create(
    Rcpp::_["clusters"] = out_clusters,
    Rcpp::_["r"] = r_out,
    Rcpp::_["counts"] = counts,
    Rcpp::_["fingerprint"] = fingerprint,
    Rcpp::_["params"] = params
  );
}

// [[Rcpp::export]]
Rcpp::List cec_cpp_gauss_univ_params_from_phiM(
    const NumericVector& z,
    const NumericMatrix& phiM,
    double lambda,
    double C) {
  const int n = phiM.nrow();
  const int r = phiM.ncol();
  const double* z_ptr = z.begin();
  const double* phi_ptr = phiM.begin();
  NumericVector counts_raw(r);
  NumericVector sum_x_raw(r);
  NumericVector sum_x2_raw(r);
  double* counts_raw_ptr = counts_raw.begin();
  double* sum_x_raw_ptr = sum_x_raw.begin();
  double* sum_x2_raw_ptr = sum_x2_raw.begin();

  for (int k = 0; k < r; ++k) {
    const double* phi_col = phi_ptr + static_cast<R_xlen_t>(k) * n;
    double count = 0.0;
    double sum_x = 0.0;
    double sum_x2 = 0.0;
    for (int i = 0; i < n; ++i) {
      const double w = phi_col[i];
      if (w != 0.0) {
        const double zi = z_ptr[i];
        count += w;
        sum_x += w * zi;
        sum_x2 += w * zi * zi;
      }
    }
    counts_raw_ptr[k] = count;
    sum_x_raw_ptr[k] = sum_x;
    sum_x2_raw_ptr[k] = sum_x2;
  }

  int r_out = 0;
  for (int k = 0; k < r; ++k) {
    if (counts_raw_ptr[k] > 0.0) {
      ++r_out;
    }
  }

  IntegerVector states(r_out);
  NumericVector nu(r_out);
  NumericVector m(r_out);
  NumericVector s(r_out);
  NumericVector varPhi(r_out);
  std::vector<int> bounded;
  bounded.reserve(r_out);
  const double s_min = 1.0 / (std::sqrt(2.0 * M_PI) * C);

  int out = 0;
  for (int k = 0; k < r; ++k) {
    const double count = counts_raw_ptr[k];
    if (count <= 0.0) {
      continue;
    }
    states[out] = out + 1;
    nu[out] = count / static_cast<double>(n);
    m[out] = sum_x_raw_ptr[k] / count;
    double var = sum_x2_raw_ptr[k] / count - m[out] * m[out];
    if (var < 0.0) {
      var = 0.0;
    }
    varPhi[out] = var;
    double sk = std::sqrt(var);
    if (sk < s_min) {
      sk = s_min;
      bounded.push_back(out + 1);
    }
    s[out] = sk;
    ++out;
  }

  IntegerVector functionBoundReached(bounded.begin(), bounded.end());
  return Rcpp::List::create(
    Rcpp::_["states"] = states,
    Rcpp::_["nu"] = nu,
    Rcpp::_["m"] = m,
    Rcpp::_["s"] = s,
    Rcpp::_["varPhi"] = varPhi,
    Rcpp::_["lambda"] = lambda,
    Rcpp::_["C"] = C,
    Rcpp::_["familyType"] = "gaussUniv",
    Rcpp::_["functionBoundReached"] = functionBoundReached
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cec_cpp_gauss_univ_perturb_warm_phiM(
    const NumericVector& phi,
    int r,
    const NumericVector& nu,
    int n) {
  std::vector<int> keep;
  keep.reserve(r);
  const double* nu_ptr = nu.begin();
  for (int k = 0; k < r; ++k) {
    if (nu_ptr[k] * static_cast<double>(n) >= 10.0) {
      keep.push_back(k);
    }
  }
  if (keep.empty()) {
    int best = 0;
    double best_nu = nu_ptr[0];
    for (int k = 1; k < r; ++k) {
      if (nu_ptr[k] > best_nu) {
        best_nu = nu_ptr[k];
        best = k;
      }
    }
    keep.push_back(best);
  }

  const int r_out = static_cast<int>(keep.size());
  NumericMatrix out = no_init(n, r_out);
  NumericVector row_sums(n);
  const double* phi_ptr = phi.begin();
  double* out_ptr = out.begin();
  double* row_sums_ptr = row_sums.begin();

  for (int out_k = 0; out_k < r_out; ++out_k) {
    const int src_k = keep[out_k];
    const double* phi_col = phi_ptr + static_cast<R_xlen_t>(src_k) * n;
    double* out_col = out_ptr + static_cast<R_xlen_t>(out_k) * n;
    for (int i = 0; i < n; ++i) {
      const double value = phi_col[i] + 0.1 * R::unif_rand();
      out_col[i] = value;
      row_sums_ptr[i] += value;
    }
  }

  for (int out_k = 0; out_k < r_out; ++out_k) {
    double* out_col = out_ptr + static_cast<R_xlen_t>(out_k) * n;
    for (int i = 0; i < n; ++i) {
      const double denom = row_sums_ptr[i] > 0.0 ? row_sums_ptr[i] : 1.0;
      out_col[i] /= denom;
    }
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cec_cpp_gauss_univ_expand_phi_random(
    const NumericVector& phi,
    int n) {
  const int phi_len = phi.size();
  const int r_old = phi_len / n;
  NumericMatrix out(n, r_old + 1);

  const double* phi_ptr = phi.begin();
  double* out_ptr = out.begin();

  for (int i = 0; i < phi_len; ++i) {
    out_ptr[i] = phi_ptr[i] * R::unif_rand();
  }

  double* new_col = out_ptr + static_cast<R_xlen_t>(r_old) * n;
  for (int i = 0; i < n; ++i) {
    new_col[i] = R::unif_rand();
  }

  return out;
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
  (void)lambda;
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int r = means.n_rows;

  arma::mat out(n, r, arma::fill::zeros);
  const double log2pi = std::log(2.0 * M_PI);

  for (int k = 0; k < r; ++k) {
    arma::mat Sigma = as<arma::mat>(sigma_list[k]);
    arma::mat L = cec_chol_safe(Sigma);
    double logdet = 2.0 * arma::sum(arma::log(L.diag()));
    arma::mat centered = X.each_row() - means.row(k);
    arma::mat solved = arma::solve(arma::trimatl(L), centered.t());
    arma::rowvec quad = arma::sum(arma::square(solved), 0);
    out.col(k) = -0.5 * (p * log2pi + logdet) - 0.5 * quad.t();
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
        if (code < 1 || code > static_cast<int>(probs.n_rows)) {
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
