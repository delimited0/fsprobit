#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {

const double VAR_FLOOR = 1e-12;
const double DENOM_FLOOR = 1e-12;
const double TOL = 1e-8;
const int MAX_ITER = 1000;

double mills_upper(double alpha) {
  const double log_phi = R::dnorm(alpha, 0.0, 1.0, true);
  const double log_tail = R::pnorm(alpha, 0.0, 1.0, false, true);
  return std::exp(log_phi - log_tail);
}

double mills_lower(double beta) {
  const double log_phi = R::dnorm(beta, 0.0, 1.0, true);
  const double log_cdf = R::pnorm(beta, 0.0, 1.0, true, true);
  return std::exp(log_phi - log_cdf);
}

void truncated_normal_moments(
    double cavity_mean,
    double cavity_var,
    bool lower_zero,
    double& tilted_mean,
    double& tilted_var) {

  if (!R_finite(cavity_mean) || !R_finite(cavity_var) ||
      cavity_var <= VAR_FLOOR) {
    stop("Invalid cavity moments in EPMNP update");
  }

  const double sd = std::sqrt(cavity_var);
  double var_multiplier;

  if (lower_zero) {
    const double alpha = -cavity_mean / sd;
    const double lambda = mills_upper(alpha);
    tilted_mean = cavity_mean + sd * lambda;
    var_multiplier = 1.0 + alpha * lambda - lambda * lambda;
  } else {
    const double beta = -cavity_mean / sd;
    const double lambda = mills_lower(beta);
    tilted_mean = cavity_mean - sd * lambda;
    var_multiplier = 1.0 - beta * lambda - lambda * lambda;
  }

  tilted_var = cavity_var * var_multiplier;
  if (!R_finite(tilted_mean) || !R_finite(tilted_var) ||
      tilted_var <= VAR_FLOOR) {
    tilted_var = VAR_FLOOR;
  }
}

double symmetrize_and_diff(NumericMatrix& Sigma, const NumericMatrix& Sigma_old) {
  const int m = Sigma.nrow();
  double max_diff = 0.0;

  for (int r = 0; r < m; ++r) {
    for (int c = r + 1; c < m; ++c) {
      const double value = 0.5 * (Sigma(r, c) + Sigma(c, r));
      Sigma(r, c) = value;
      Sigma(c, r) = value;
    }
  }

  for (int r = 0; r < m; ++r) {
    for (int c = 0; c < m; ++c) {
      max_diff = std::max(max_diff, std::abs(Sigma(r, c) - Sigma_old(r, c)));
    }
  }

  return max_diff;
}

} // namespace

//' Sparse EP moment approximation for multinomial probit
//'
//' @param mu Mean vector of the latent utilities.
//' @param Sigma Covariance matrix of the latent utilities.
//' @param choice_index Choice index: 0 for the base choice, otherwise 1..m for
//'   the selected non-base latent utility.
//' @return A list with approximate truncated-normal `mu` and `Sigma`.
// [[Rcpp::export]]
List epmnp(NumericVector mu, NumericMatrix Sigma, int choice_index) {
  const int m = mu.size();

  if (m < 1) {
    stop("mu must have positive length");
  }
  if (Sigma.nrow() != m || Sigma.ncol() != m) {
    stop("Sigma must be a square matrix with dimensions matching mu");
  }
  if (choice_index < 0 || choice_index > m) {
    stop("choice_index must be between 0 and length(mu)");
  }

  NumericVector mu_ep = clone(mu);
  NumericMatrix Sigma_ep = clone(Sigma);
  NumericVector tau(m);
  NumericVector eta(m);

  const bool base_choice = choice_index == 0;
  const int y = choice_index - 1;

  for (int iter = 0; iter < MAX_ITER; ++iter) {
    NumericVector mu_old = clone(mu_ep);
    NumericMatrix Sigma_old = clone(Sigma_ep);
    double max_site_change = 0.0;

    for (int j = 0; j < m; ++j) {
      double d = 0.0;
      double s2 = 0.0;
      NumericVector g(m);
      bool lower_zero = true;

      if (base_choice) {
        d = mu_ep[j];
        s2 = Sigma_ep(j, j);
        lower_zero = false;
        for (int r = 0; r < m; ++r) {
          g[r] = Sigma_ep(r, j);
        }
      } else if (j == y) {
        d = mu_ep[y];
        s2 = Sigma_ep(y, y);
        lower_zero = true;
        for (int r = 0; r < m; ++r) {
          g[r] = Sigma_ep(r, y);
        }
      } else {
        d = mu_ep[y] - mu_ep[j];
        s2 = Sigma_ep(y, y) - 2.0 * Sigma_ep(y, j) + Sigma_ep(j, j);
        lower_zero = true;
        for (int r = 0; r < m; ++r) {
          g[r] = Sigma_ep(r, y) - Sigma_ep(r, j);
        }
      }

      if (!R_finite(s2) || s2 <= VAR_FLOOR) {
        s2 = VAR_FLOOR;
      }

      const double marginal_tau = 1.0 / s2;
      const double marginal_eta = d / s2;
      const double cavity_tau = marginal_tau - tau[j];
      const double cavity_eta = marginal_eta - eta[j];

      if (!R_finite(cavity_tau) || cavity_tau <= VAR_FLOOR) {
        continue;
      }

      const double cavity_var = 1.0 / cavity_tau;
      const double cavity_mean = cavity_eta / cavity_tau;
      double tilted_mean = 0.0;
      double tilted_var = 0.0;

      truncated_normal_moments(
        cavity_mean,
        cavity_var,
        lower_zero,
        tilted_mean,
        tilted_var
      );

      const double tau_new = 1.0 / tilted_var - cavity_tau;
      const double eta_new = tilted_mean / tilted_var - cavity_eta;
      const double delta_tau = tau_new - tau[j];
      const double delta_eta = eta_new - eta[j];

      if (!R_finite(delta_tau) || !R_finite(delta_eta)) {
        continue;
      }

      tau[j] += delta_tau;
      eta[j] += delta_eta;

      const double denom = 1.0 + delta_tau * s2;
      if (!R_finite(denom) || std::abs(denom) <= DENOM_FLOOR) {
        stop("Numerically unstable EPMNP covariance update");
      }

      const double cov_scale = delta_tau / denom;
      const double mean_scale = (delta_eta - delta_tau * d) / denom;

      for (int r = 0; r < m; ++r) {
        mu_ep[r] += mean_scale * g[r];
      }

      for (int r = 0; r < m; ++r) {
        for (int c = 0; c < m; ++c) {
          Sigma_ep(r, c) -= cov_scale * g[r] * g[c];
        }
      }

      max_site_change = std::max(max_site_change, std::abs(delta_tau));
      max_site_change = std::max(max_site_change, std::abs(delta_eta));
    }

    double max_moment_change = symmetrize_and_diff(Sigma_ep, Sigma_old);
    for (int j = 0; j < m; ++j) {
      max_moment_change = std::max(max_moment_change, std::abs(mu_ep[j] - mu_old[j]));
    }

    if (std::max(max_site_change, max_moment_change) < TOL) {
      break;
    }
  }

  return List::create(
    Named("mu") = mu_ep,
    Named("Sigma") = Sigma_ep
  );
}
