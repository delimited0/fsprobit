#include <Rcpp.h>
#include <string>
#include <vector>
#include "epmnp.h"

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
  if (E_method != "EPMNP") {
    stop("mnp_probit_accelerated currently supports C++ accumulator E_method = 'EPMNP' only");
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
//' @param E_method E-step moment method. Currently only `"EPMNP"`.
//' @return A list of accumulated M-step sufficient statistics.
// [[Rcpp::export]]
List mnp_probit_accumulate_cpp(
    NumericVector X,
    IntegerVector Y,
    IntegerVector obs_set,
    NumericMatrix beta,
    NumericMatrix Sigma,
    NumericMatrix Precision,
    std::string E_method) {

  validate_method(E_method);

  IntegerVector dims = X.attr("dim");
  if (dims.size() != 3) {
    stop("X must be a three-dimensional array");
  }

  const int n = dims[0];
  const int m = dims[1];
  const int p = dims[2];
  const int n_used = obs_set.size();

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
  NumericVector mu_X_sum(m * m * p);
  NumericVector X_cross_sum(m * m * p * p);

  mu_X_sum.attr("dim") = IntegerVector::create(m, m, p);
  X_cross_sum.attr("dim") = IntegerVector::create(m, m, p, p);

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

    NumericVector mu_ep(m);
    NumericMatrix Sigma_ep = clone(Sigma);

    for (int r = 0; r < m; ++r) {
      double value = 0.0;
      for (int k = 0; k < p; ++k) {
        value += X[x_index(obs_idx, r, k, n, m)] * beta(k, 0);
      }
      mu_ep[r] = value;
    }

    epmnp_moments_inplace(mu_ep, Sigma_ep, choice_index);

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
        b_value += tXPrecision[k + p * c] * mu_ep[c];
      }
      gls_b(k, 0) += b_value;
    }

    for (int r = 0; r < m; ++r) {
      for (int c = 0; c < m; ++c) {
        E_second_moment(r, c) += Sigma_ep(r, c) + mu_ep[r] * mu_ep[c];
      }
    }

    for (int k = 0; k < p; ++k) {
      for (int r = 0; r < m; ++r) {
        for (int c = 0; c < m; ++c) {
          mu_X_sum[array3_index(r, c, k, m)] +=
            mu_ep[r] * X[x_index(obs_idx, c, k, n, m)];
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

  return List::create(
    Named("gls_a") = gls_a,
    Named("gls_b") = gls_b,
    Named("E_second_moment") = E_second_moment,
    Named("mu_X_sum") = mu_X_sum,
    Named("X_cross_sum") = X_cross_sum,
    Named("n_used") = n_used
  );
}
