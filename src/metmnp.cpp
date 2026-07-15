// Minimax exponentially tilted accept-reject sampling for Gaussian moments,
// specialized to reference-identified multinomial probit constraints.
//
// The algorithm follows Botev (2017), "The normal law under linear
// restrictions: simulation and estimation via minimax tilting", with
// reference to the MIT-licensed Python implementation by Paul Brunzema:
// https://github.com/brunzema/truncated-mvn-sampler
// Copyright (c) 2021 Paul Brunzema. See inst/COPYRIGHTS.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include "metmnp.h"

using namespace Rcpp;

namespace {

const double LOG_2PI = 1.837877066409345483560659472811;
const double OPT_TOL = 1e-8;
const int MAX_OPT_ITER = 500;
const long long MIN_PROPOSAL_LIMIT = 1000000LL;
const long long PROPOSALS_PER_SAMPLE = 10000LL;

typedef std::vector<double> Vec;
typedef std::vector<Vec> Mat;

double dot(const Vec& a, const Vec& b) {
  double out = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) out += a[i] * b[i];
  return out;
}

double norm2(const Vec& x) { return std::sqrt(dot(x, x)); }

double log_upper_tail(double x) {
  return R::pnorm(x, 0.0, 1.0, 0, 1);
}

double inverse_mills(double x, double log_tail) {
  const double log_phi = -0.5 * x * x - 0.5 * LOG_2PI;
  return std::exp(log_phi - log_tail);
}

bool solve_linear(Mat a, Vec b, Vec& x) {
  const int n = b.size();
  x.assign(n, 0.0);
  for (int k = 0; k < n; ++k) {
    int pivot = k;
    double largest = std::abs(a[k][k]);
    for (int i = k + 1; i < n; ++i) {
      if (std::abs(a[i][k]) > largest) {
        largest = std::abs(a[i][k]);
        pivot = i;
      }
    }
    if (!R_finite(largest) || largest < 1e-12) return false;
    if (pivot != k) {
      std::swap(a[pivot], a[k]);
      std::swap(b[pivot], b[k]);
    }
    for (int i = k + 1; i < n; ++i) {
      const double scale = a[i][k] / a[k][k];
      a[i][k] = 0.0;
      for (int j = k + 1; j < n; ++j) a[i][j] -= scale * a[k][j];
      b[i] -= scale * b[k];
    }
  }
  for (int i = n - 1; i >= 0; --i) {
    double value = b[i];
    for (int j = i + 1; j < n; ++j) value -= a[i][j] * x[j];
    x[i] = value / a[i][i];
    if (!R_finite(x[i])) return false;
  }
  return true;
}

struct CholeskyResult {
  Mat lower;
  Vec bounds;
  std::vector<int> permutation;
};

CholeskyResult reordered_cholesky(Mat covariance, Vec bounds) {
  const int d = bounds.size();
  Mat lower(d, Vec(d, 0.0));
  Vec expected(d, 0.0);
  std::vector<int> permutation(d);
  std::iota(permutation.begin(), permutation.end(), 0);

  for (int j = 0; j < d; ++j) {
    int best = j;
    double best_log_probability = std::numeric_limits<double>::infinity();
    for (int i = j; i < d; ++i) {
      double variance = covariance[i][i];
      double shift = 0.0;
      for (int k = 0; k < j; ++k) {
        variance -= lower[i][k] * lower[i][k];
        shift += lower[i][k] * expected[k];
      }
      if (variance <= 0.0 || !R_finite(variance)) {
        stop("Sigma must be positive definite");
      }
      const double standardized = (bounds[i] - shift) / std::sqrt(variance);
      const double probability = log_upper_tail(standardized);
      if (probability < best_log_probability) {
        best_log_probability = probability;
        best = i;
      }
    }

    if (best != j) {
      std::swap(covariance[best], covariance[j]);
      for (int i = 0; i < d; ++i) std::swap(covariance[i][best], covariance[i][j]);
      std::swap(lower[best], lower[j]);
      std::swap(bounds[best], bounds[j]);
      std::swap(permutation[best], permutation[j]);
    }

    double diagonal = covariance[j][j];
    for (int k = 0; k < j; ++k) diagonal -= lower[j][k] * lower[j][k];
    if (diagonal <= 0.0 || !R_finite(diagonal)) {
      stop("Sigma must be positive definite");
    }
    lower[j][j] = std::sqrt(diagonal);
    for (int i = j + 1; i < d; ++i) {
      double value = covariance[i][j];
      for (int k = 0; k < j; ++k) value -= lower[i][k] * lower[j][k];
      lower[i][j] = value / lower[j][j];
    }

    double shift = 0.0;
    for (int k = 0; k < j; ++k) shift += lower[j][k] * expected[k];
    const double a = (bounds[j] - shift) / lower[j][j];
    expected[j] = inverse_mills(a, log_upper_tail(a));
  }
  return {lower, bounds, permutation};
}

struct TiltSystem {
  Mat lower;
  Vec bounds;
  int d;

  void evaluate(const Vec& y, Vec& gradient, Mat* jacobian) const {
    const int q = d - 1;
    Vec x(d, 0.0), tilt(d, 0.0), mills(d), derivative(d);
    for (int i = 0; i < q; ++i) {
      x[i] = y[i];
      tilt[i] = y[q + i];
    }
    for (int i = 0; i < d; ++i) {
      double conditional = 0.0;
      for (int j = 0; j < i; ++j) conditional += lower[i][j] * x[j];
      const double a = bounds[i] - tilt[i] - conditional;
      const double log_probability = log_upper_tail(a);
      mills[i] = inverse_mills(a, log_probability);
      derivative[i] = a * mills[i] - mills[i] * mills[i];
      if (!R_finite(mills[i]) || !R_finite(derivative[i])) {
        stop("MET tilting equations became numerically non-finite");
      }
    }

    gradient.assign(2 * q, 0.0);
    for (int j = 0; j < q; ++j) {
      gradient[j] = -tilt[j];
      for (int i = j + 1; i < d; ++i) gradient[j] += lower[i][j] * mills[i];
      gradient[q + j] = tilt[j] - x[j] + mills[j];
    }
    if (jacobian == NULL) return;

    jacobian->assign(2 * q, Vec(2 * q, 0.0));
    for (int r = 0; r < q; ++r) {
      for (int c = 0; c < q; ++c) {
        double xx = 0.0;
        for (int i = std::max(r, c) + 1; i < d; ++i) {
          xx += lower[i][r] * derivative[i] * lower[i][c];
        }
        (*jacobian)[r][c] = xx;
        const double mx = (r == c ? -1.0 : 0.0) + derivative[r] * lower[r][c];
        (*jacobian)[q + r][c] = mx;
        (*jacobian)[c][q + r] = mx;
      }
      (*jacobian)[q + r][q + r] = 1.0 + derivative[r];
    }
  }
};

Vec dogleg_step(const Vec& gradient, const Mat& jacobian, double radius) {
  const int n = gradient.size();
  Vec rhs(n), newton;
  for (int i = 0; i < n; ++i) rhs[i] = -gradient[i];
  const bool have_newton = solve_linear(jacobian, rhs, newton);
  if (have_newton && norm2(newton) <= radius) return newton;

  Vec objective_gradient(n, 0.0);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) objective_gradient[j] += jacobian[i][j] * gradient[i];
  }
  Vec jg(n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) jg[i] += jacobian[i][j] * objective_gradient[j];
  }
  const double denominator = dot(jg, jg);
  double alpha = denominator > 0.0 ? dot(objective_gradient, objective_gradient) / denominator : 1.0;
  Vec cauchy(n);
  for (int i = 0; i < n; ++i) cauchy[i] = -alpha * objective_gradient[i];
  const double cauchy_norm = norm2(cauchy);
  if (!have_newton || cauchy_norm >= radius) {
    const double scale = radius / std::max(cauchy_norm, 1e-300);
    for (double& value : cauchy) value *= scale;
    return cauchy;
  }

  Vec direction(n);
  for (int i = 0; i < n; ++i) direction[i] = newton[i] - cauchy[i];
  const double aa = dot(direction, direction);
  const double bb = 2.0 * dot(cauchy, direction);
  const double cc = dot(cauchy, cauchy) - radius * radius;
  const double tau = (-bb + std::sqrt(std::max(0.0, bb * bb - 4.0 * aa * cc))) / (2.0 * aa);
  for (int i = 0; i < n; ++i) cauchy[i] += tau * direction[i];
  return cauchy;
}

Vec solve_tilting(const TiltSystem& system) {
  const int n = 2 * (system.d - 1);
  Vec y(n, 0.0), gradient;
  Mat jacobian;
  double radius = 1.0;
  system.evaluate(y, gradient, &jacobian);

  for (int iteration = 0; iteration < MAX_OPT_ITER; ++iteration) {
    if (norm2(gradient) <= OPT_TOL) return y;
    Vec step = dogleg_step(gradient, jacobian, radius);
    Vec candidate = y;
    for (int i = 0; i < n; ++i) candidate[i] += step[i];
    Vec candidate_gradient;
    system.evaluate(candidate, candidate_gradient, NULL);

    Vec linear = gradient;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) linear[i] += jacobian[i][j] * step[j];
    }
    const double actual = dot(gradient, gradient) - dot(candidate_gradient, candidate_gradient);
    const double predicted = dot(gradient, gradient) - dot(linear, linear);
    const double ratio = predicted > 0.0 ? actual / predicted : -1.0;
    if (ratio < 0.25) radius *= 0.25;
    else if (ratio > 0.75 && std::abs(norm2(step) - radius) < 1e-8 * std::max(1.0, radius)) {
      radius = std::min(100.0, 2.0 * radius);
    }
    if (ratio > 1e-4) {
      y.swap(candidate);
      gradient.swap(candidate_gradient);
      system.evaluate(y, gradient, &jacobian);
    }
    if (radius < 1e-12) break;
  }
  stop("MET minimax tilting optimization did not converge");
}

double truncated_normal_lower(double lower) {
  if (lower > 0.66) {
    const double c = 0.5 * lower * lower;
    while (true) {
      const double x = c + R::rexp(1.0);
      if (R::runif(0.0, 1.0) * R::runif(0.0, 1.0) * x <= c) {
        return std::sqrt(2.0 * x);
      }
    }
  }
  while (true) {
    const double x = R::rnorm(0.0, 1.0);
    if (x >= lower) return x;
  }
}

Mat transform_covariance(const NumericMatrix& sigma, int choice_index) {
  const int d = sigma.nrow();
  Mat out(d, Vec(d, 0.0));
  if (choice_index == 0) {
    for (int i = 0; i < d; ++i) for (int j = 0; j < d; ++j) out[i][j] = sigma(i, j);
    return out;
  }
  const int y = choice_index - 1;
  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < d; ++j) {
      out[i][j] = sigma(y, y);
      if (i != y) out[i][j] -= sigma(i, y);
      if (j != y) out[i][j] -= sigma(y, j);
      if (i != y && j != y) out[i][j] += sigma(i, j);
    }
  }
  return out;
}

Vec transform_mean(const NumericVector& mu, int choice_index) {
  const int d = mu.size();
  Vec out(d);
  if (choice_index == 0) {
    for (int i = 0; i < d; ++i) out[i] = -mu[i];
  } else {
    const int y = choice_index - 1;
    for (int i = 0; i < d; ++i) out[i] = (i == y) ? mu[y] : mu[y] - mu[i];
  }
  return out;
}

Vec inverse_transform(const Vec& w, int choice_index) {
  const int d = w.size();
  Vec z(d);
  if (choice_index == 0) {
    for (int i = 0; i < d; ++i) z[i] = -w[i];
  } else {
    const int y = choice_index - 1;
    for (int i = 0; i < d; ++i) z[i] = (i == y) ? w[y] : w[y] - w[i];
  }
  return z;
}

void validate_inputs(const NumericVector& mu, const NumericMatrix& sigma,
                     int choice_index, int n_mc) {
  const int d = mu.size();
  if (d < 1) stop("mu must have positive length");
  if (sigma.nrow() != d || sigma.ncol() != d) stop("Sigma dimensions must match mu");
  if (choice_index < 0 || choice_index > d) stop("choice_index must be between 0 and length(mu)");
  if (n_mc < 2) stop("n_mc must be at least 2");
  for (int i = 0; i < d; ++i) {
    if (!R_finite(mu[i])) stop("mu must contain only finite values");
    for (int j = 0; j < d; ++j) {
      if (!R_finite(sigma(i, j))) stop("Sigma must contain only finite values");
      const double scale = std::max(1.0, std::max(std::abs(sigma(i, j)), std::abs(sigma(j, i))));
      if (std::abs(sigma(i, j) - sigma(j, i)) > 1e-10 * scale) stop("Sigma must be symmetric");
    }
  }
}

} // namespace

List metmnp_moments(const NumericVector& mu, const NumericMatrix& Sigma,
                    int choice_index, int n_mc) {
  validate_inputs(mu, Sigma, choice_index, n_mc);
  const int d = mu.size();
  const Vec transformed_mu = transform_mean(mu, choice_index);
  Mat transformed_sigma = transform_covariance(Sigma, choice_index);
  Vec centered_bounds(d);
  for (int i = 0; i < d; ++i) centered_bounds[i] = -transformed_mu[i];
  CholeskyResult factor = reordered_cholesky(transformed_sigma, centered_bounds);

  Vec diagonal(d);
  Mat scaled_lower(d, Vec(d, 0.0));
  Vec scaled_bounds(d);
  for (int i = 0; i < d; ++i) {
    diagonal[i] = factor.lower[i][i];
    scaled_bounds[i] = factor.bounds[i] / diagonal[i];
    for (int j = 0; j < i; ++j) scaled_lower[i][j] = factor.lower[i][j] / diagonal[i];
  }

  Vec tilt(d, 0.0), saddle_x(d, 0.0);
  if (d > 1) {
    TiltSystem system{scaled_lower, scaled_bounds, d};
    Vec solution = solve_tilting(system);
    for (int i = 0; i < d - 1; ++i) {
      saddle_x[i] = solution[i];
      tilt[i] = solution[d - 1 + i];
    }
    // If the unconstrained stationary point is outside the truncation region,
    // psi_star remains a valid (but less tight) likelihood-ratio upper bound.
    // This preserves exactness at the cost of a potentially lower acceptance
    // rate, matching the simulation references for Botev's method.
  }

  double psi_star = 0.0;
  for (int i = 0; i < d; ++i) {
    double conditional = 0.0;
    for (int j = 0; j < i; ++j) conditional += scaled_lower[i][j] * saddle_x[j];
    const double a = scaled_bounds[i] - tilt[i] - conditional;
    psi_star += log_upper_tail(a) + 0.5 * tilt[i] * tilt[i] - saddle_x[i] * tilt[i];
  }

  Vec mean(d, 0.0);
  NumericMatrix covariance_sum(d, d);
  int accepted = 0;
  long long proposals = 0;
  const long long proposal_limit = std::max(MIN_PROPOSAL_LIMIT,
    PROPOSALS_PER_SAMPLE * static_cast<long long>(n_mc));
  while (accepted < n_mc) {
    if (++proposals > proposal_limit) stop("MET acceptance rate is too low to obtain n_mc samples");
    Vec standardized(d, 0.0);
    double log_ratio = 0.0;
    for (int i = 0; i < d; ++i) {
      double conditional = 0.0;
      for (int j = 0; j < i; ++j) conditional += scaled_lower[i][j] * standardized[j];
      const double lower = scaled_bounds[i] - tilt[i] - conditional;
      standardized[i] = tilt[i] + truncated_normal_lower(lower);
      log_ratio += log_upper_tail(lower) + 0.5 * tilt[i] * tilt[i] - tilt[i] * standardized[i];
    }
    double rejection_log = psi_star - log_ratio;
    if (rejection_log < -1e-7) stop("MET likelihood-ratio bound was violated numerically");
    rejection_log = std::max(0.0, rejection_log);
    if (R::rexp(1.0) <= rejection_log) continue;

    Vec centered_permuted(d, 0.0), w(d);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j <= i; ++j) centered_permuted[i] += factor.lower[i][j] * standardized[j];
      const int original = factor.permutation[i];
      w[original] = centered_permuted[i] + transformed_mu[original];
    }
    Vec sample = inverse_transform(w, choice_index);
    ++accepted;
    Vec delta(d), delta_after(d);
    for (int i = 0; i < d; ++i) {
      delta[i] = sample[i] - mean[i];
      mean[i] += delta[i] / accepted;
      delta_after[i] = sample[i] - mean[i];
    }
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) covariance_sum(i, j) += delta[i] * delta_after[j];
    }
  }

  NumericVector mean_out(d);
  NumericMatrix covariance(d, d);
  for (int i = 0; i < d; ++i) {
    mean_out[i] = mean[i];
    for (int j = 0; j < d; ++j) covariance(i, j) = covariance_sum(i, j) / (n_mc - 1.0);
  }
  return List::create(Named("mu") = mean_out, Named("Sigma") = covariance);
}

//' Minimax-tilted moments for reference-identified multinomial probit
//'
//' @param mu Mean vector of the relative latent utilities.
//' @param Sigma Covariance matrix of the relative latent utilities.
//' @param choice_index Choice index: 0 for the base choice, otherwise 1..m.
//' @param n_mc Number of accepted independent Monte Carlo samples.
//' @return A list containing Monte Carlo estimates `mu` and `Sigma`.
// [[Rcpp::export]]
List metmnp(NumericVector mu, NumericMatrix Sigma, int choice_index, int n_mc) {
  RNGScope scope;
  return metmnp_moments(mu, Sigma, choice_index, n_mc);
}
