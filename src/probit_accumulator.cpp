#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <vector>
#include "epmnp.h"
#include "hmcmnp.h"
#include "metmnp.h"

using namespace Rcpp;

namespace {

int x_index(int i, int r, int k, int n, int m) {
  return i + n * r + n * m * k;
}

int array3_index(int r, int c, int k, int m) {
  return r + m * c + m * m * k;
}

int array4_index(int r, int c, int k, int l, int m, int p) {
  return r + m * c + m * m * k + m * m * p * l;
}

void validate_method(const std::string& E_method) {
  if (E_method != "EPMNP" && E_method != "HMC" && E_method != "MET") {
    stop("mnp_probit_accelerated supports E_method = 'EPMNP', 'HMC', or 'MET'");
  }
}

} // namespace

//' Accumulate E-step sufficient statistics for multinomial probit
//'
//' @param X \eqn{n \times m \times p} covariate array.
//' @param Y Integer response codes with base choice coded as 1.
//' @param obs_set One-based observation indices to include.
//' @param beta Current or effective coefficient matrix.
//' @param Sigma Current or effective latent utility covariance.
//' @param Precision Current latent utility precision.
//' @param E_method E-step moment method: `"EPMNP"`, `"HMC"`, or `"MET"`.
//' @param n_mc Number of retained Monte Carlo samples for HMC or MET.
//' @param return_cross_tensors Whether to return the legacy contraction tensors.
//'   Set to `false` to return posterior means for the optimized second pass.
//' @return A list of accumulated M-step sufficient statistics.
// [[Rcpp::export]]
List mnp_probit_accumulate_cpp(
    NumericVector X,
    IntegerVector Y,
    IntegerVector obs_set,
    NumericMatrix beta,
    NumericMatrix Sigma,
    NumericMatrix Precision,
    std::string E_method,
    int n_mc,
    bool return_cross_tensors = true) {

  validate_method(E_method);
  if ((E_method == "HMC" || E_method == "MET") && n_mc < 2) {
    stop("n_mc must be at least 2 for HMC or MET");
  }

  IntegerVector dims = X.attr("dim");
  if (dims.size() != 3) {
    stop("X must be a three-dimensional array");
  }

  const int n = dims[0];
  const int m = dims[1];
  const int p = dims[2];
  const int n_used = obs_set.size();

  if (n_used < 1) {
    stop("obs_set must contain at least one observation");
  }

  if (Y.size() != n) {
    stop("Y length must match the first dimension of X");
  }
  if (beta.nrow() != p || beta.ncol() != 1) {
    stop("beta must be a p x 1 matrix");
  }
  if (Sigma.nrow() != m || Sigma.ncol() != m) {
    stop("Sigma dimensions must match the second dimension of X");
  }
  if (Precision.nrow() != m || Precision.ncol() != m) {
    stop("Precision dimensions must match the second dimension of X");
  }

  NumericMatrix gls_a(p, p);
  NumericMatrix gls_b(p, 1);
  NumericMatrix E_second_moment(m, m);
  NumericVector mu_X_sum;
  NumericVector X_cross_sum;
  NumericMatrix moment_mu_by_obs;

  if (return_cross_tensors) {
    mu_X_sum = NumericVector(m * m * p);
    X_cross_sum = NumericVector(m * m * p * p);
    mu_X_sum.attr("dim") = IntegerVector::create(m, m, p);
    X_cross_sum.attr("dim") = IntegerVector::create(m, m, p, p);
  } else {
    // Columns are observations so each observation's moments are contiguous.
    moment_mu_by_obs = NumericMatrix(m, n_used);
  }

  std::vector<double> tXPrecision(p * m);

  for (int obs_pos = 0; obs_pos < n_used; ++obs_pos) {
    const int obs_idx = obs_set[obs_pos] - 1;
    if (obs_idx < 0 || obs_idx >= n) {
      stop("obs_set contains an out-of-range observation index");
    }

    const int choice_index = Y[obs_idx] - 1;
    if (choice_index < 0 || choice_index > m) {
      stop("Y contains a choice outside the valid range for X");
    }

    NumericVector moment_mu(m);
    NumericMatrix moment_sigma(m, m);

    for (int r = 0; r < m; ++r) {
      double value = 0.0;
      for (int k = 0; k < p; ++k) {
        value += X[x_index(obs_idx, r, k, n, m)] * beta(k, 0);
      }
      moment_mu[r] = value;
    }

    if (E_method == "EPMNP") {
      std::copy(Sigma.begin(), Sigma.end(), moment_sigma.begin());
      epmnp_moments_inplace(moment_mu, moment_sigma, choice_index);
    } else if (E_method == "HMC") {
      List moments = hmcmnp_moments(moment_mu, Precision, choice_index, n_mc);
      moment_mu = as<NumericVector>(moments["mu"]);
      moment_sigma = as<NumericMatrix>(moments["Sigma"]);
    } else {
      List moments = metmnp_moments(moment_mu, Sigma, choice_index, n_mc);
      moment_mu = as<NumericVector>(moments["mu"]);
      moment_sigma = as<NumericMatrix>(moments["Sigma"]);
    }

    if (!return_cross_tensors) {
      for (int r = 0; r < m; ++r) {
        moment_mu_by_obs(r, obs_pos) = moment_mu[r];
      }
    }

    for (int k = 0; k < p; ++k) {
      for (int c = 0; c < m; ++c) {
        double value = 0.0;
        for (int r = 0; r < m; ++r) {
          value += X[x_index(obs_idx, r, k, n, m)] * Precision(r, c);
        }
        tXPrecision[k + p * c] = value;
      }
    }

    for (int k = 0; k < p; ++k) {
      for (int l = 0; l < p; ++l) {
        double value = 0.0;
        for (int c = 0; c < m; ++c) {
          value += tXPrecision[k + p * c] * X[x_index(obs_idx, c, l, n, m)];
        }
        gls_a(k, l) += value;
      }

      double b_value = 0.0;
      for (int c = 0; c < m; ++c) {
        b_value += tXPrecision[k + p * c] * moment_mu[c];
      }
      gls_b(k, 0) += b_value;
    }

    for (int r = 0; r < m; ++r) {
      for (int c = 0; c < m; ++c) {
        E_second_moment(r, c) +=
          moment_sigma(r, c) + moment_mu[r] * moment_mu[c];
      }
    }

    if (return_cross_tensors) {
      for (int k = 0; k < p; ++k) {
        for (int r = 0; r < m; ++r) {
          for (int c = 0; c < m; ++c) {
            mu_X_sum[array3_index(r, c, k, m)] +=
              moment_mu[r] * X[x_index(obs_idx, c, k, n, m)];
          }
        }
      }

      for (int k = 0; k < p; ++k) {
        for (int l = 0; l < p; ++l) {
          for (int r = 0; r < m; ++r) {
            const double x_rk = X[x_index(obs_idx, r, k, n, m)];
            for (int c = 0; c < m; ++c) {
              X_cross_sum[array4_index(r, c, k, l, m, p)] +=
                x_rk * X[x_index(obs_idx, c, l, n, m)];
            }
          }
        }
      }
    }
  }

  return List::create(
    Named("gls_a") = gls_a,
    Named("gls_b") = gls_b,
    Named("E_second_moment") = E_second_moment,
    Named("mu_X_sum") = mu_X_sum,
    Named("X_cross_sum") = X_cross_sum,
    Named("moment_mu") = moment_mu_by_obs,
    Named("n_used") = n_used
  );
}

//' Complete the expected covariance statistic after estimating beta
//'
//' @param X \eqn{n \times m \times p} covariate array.
//' @param obs_set One-based observation indices used by the E-step.
//' @param moment_mu \eqn{m \times n_used} posterior mean matrix, with one
//'   observation per column.
//' @param beta_new Newly estimated \eqn{p \times 1} coefficient matrix.
//' @param E_second_moment Sum of posterior second moments from the E-step.
//' @return The averaged expected residual covariance matrix.
// [[Rcpp::export]]
NumericMatrix mnp_probit_residual_cov_cpp(
    NumericVector X,
    IntegerVector obs_set,
    NumericMatrix moment_mu,
    NumericMatrix beta_new,
    NumericMatrix E_second_moment) {

  IntegerVector dims = X.attr("dim");
  if (dims.size() != 3) stop("X must be a three-dimensional array");

  const int n = dims[0];
  const int m = dims[1];
  const int p = dims[2];
  const int n_used = obs_set.size();

  if (n_used < 1) stop("obs_set must contain at least one observation");
  if (moment_mu.nrow() != m || moment_mu.ncol() != n_used) {
    stop("moment_mu must have dimensions m x length(obs_set)");
  }
  if (beta_new.nrow() != p || beta_new.ncol() != 1) {
    stop("beta_new must be a p x 1 matrix");
  }
  if (E_second_moment.nrow() != m || E_second_moment.ncol() != m) {
    stop("E_second_moment must be an m x m matrix");
  }

  NumericMatrix covariance = clone(E_second_moment);
  std::vector<double> x_beta(m);

  for (int k = 0; k < p; ++k) {
    if (!R_finite(beta_new(k, 0))) stop("beta_new must contain only finite values");
  }
  for (int r = 0; r < m; ++r) {
    for (int c = 0; c < m; ++c) {
      if (!R_finite(E_second_moment(r, c))) {
        stop("E_second_moment must contain only finite values");
      }
    }
  }

  for (int obs_pos = 0; obs_pos < n_used; ++obs_pos) {
    const int obs_idx = obs_set[obs_pos] - 1;
    if (obs_idx < 0 || obs_idx >= n) {
      stop("obs_set contains an out-of-range observation index");
    }

    for (int r = 0; r < m; ++r) {
      if (!R_finite(moment_mu(r, obs_pos))) {
        stop("moment_mu must contain only finite values");
      }
      double value = 0.0;
      for (int k = 0; k < p; ++k) {
        value += X[x_index(obs_idx, r, k, n, m)] * beta_new(k, 0);
      }
      if (!R_finite(value)) {
        stop("X and beta_new must produce finite linear predictors");
      }
      x_beta[r] = value;
    }

    for (int r = 0; r < m; ++r) {
      for (int c = 0; c < m; ++c) {
        covariance(r, c) +=
          -moment_mu(r, obs_pos) * x_beta[c]
          -x_beta[r] * moment_mu(c, obs_pos)
          +x_beta[r] * x_beta[c];
      }
    }
  }

  for (int r = 0; r < m; ++r) {
    covariance(r, r) /= n_used;
    for (int c = r + 1; c < m; ++c) {
      const double value =
        0.5 * (covariance(r, c) + covariance(c, r)) / n_used;
      covariance(r, c) = value;
      covariance(c, r) = value;
    }
  }
  return covariance;
}
