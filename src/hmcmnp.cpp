// Exact HMC for linearly truncated Gaussians, specialized to the
// reference-identified multinomial probit constraints.
//
// The trajectory and reflection implementation follows Pakman and Paninski
// (2014), with reference to the MIT-licensed implementations at
// https://github.com/aripakman/hmc-tmg and
// https://github.com/erik-a-bensen/tmg_hmc.
// Copyright (c) 2020 Ari Pakman; copyright (c) 2025 Erik A. Bensen.
// Referenced implementation portions are used under the MIT License; see
// inst/COPYRIGHTS.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "hmcmnp.h"

using namespace Rcpp;

namespace {

const double PI = 3.141592653589793238462643383279502884;
const double TRAJECTORY_TIME = PI / 2.0;
const double HIT_EPS = 1e-10;
const double FEASIBILITY_TOL = 1e-8;
const int BURN_IN = 30;
const int MAX_RETRIES = 100;
const int MAX_BOUNCES = 10000;

double dot(const std::vector<double>& a, const std::vector<double>& b) {
  double value = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) value += a[i] * b[i];
  return value;
}

NumericMatrix cholesky_upper(const NumericMatrix& precision) {
  const int m = precision.nrow();
  NumericMatrix upper(m, m);

  for (int i = 0; i < m; ++i) {
    for (int j = i; j < m; ++j) {
      double value = precision(i, j);
      for (int k = 0; k < i; ++k) value -= upper(k, i) * upper(k, j);

      if (i == j) {
        if (!R_finite(value) || value <= 0.0) {
          stop("Precision must be positive definite");
        }
        upper(i, i) = std::sqrt(value);
      } else {
        upper(i, j) = value / upper(i, i);
      }
    }
  }
  return upper;
}

std::vector<double> upper_multiply(
    const NumericMatrix& upper,
    const std::vector<double>& x) {
  const int m = upper.nrow();
  std::vector<double> result(m, 0.0);
  for (int i = 0; i < m; ++i) {
    for (int j = i; j < m; ++j) result[i] += upper(i, j) * x[j];
  }
  return result;
}

std::vector<double> upper_solve(
    const NumericMatrix& upper,
    const std::vector<double>& b) {
  const int m = upper.nrow();
  std::vector<double> x(m, 0.0);
  for (int i = m - 1; i >= 0; --i) {
    double value = b[i];
    for (int j = i + 1; j < m; ++j) value -= upper(i, j) * x[j];
    x[i] = value / upper(i, i);
  }
  return x;
}

std::vector<double> lower_solve_from_upper(
    const NumericMatrix& upper,
    const std::vector<double>& b) {
  const int m = upper.nrow();
  std::vector<double> x(m, 0.0);
  for (int i = 0; i < m; ++i) {
    double value = b[i];
    for (int j = 0; j < i; ++j) value -= upper(j, i) * x[j];
    x[i] = value / upper(i, i);
  }
  return x;
}

struct LinearConstraint {
  std::vector<double> normal;
  double offset;
  double norm_squared;
};

std::vector<LinearConstraint> build_constraints(
    const NumericVector& mu,
    const NumericMatrix& upper,
    int choice_index) {
  const int m = mu.size();
  std::vector<LinearConstraint> constraints;
  constraints.reserve(m);

  for (int j = 0; j < m; ++j) {
    std::vector<double> f(m, 0.0);
    if (choice_index == 0) {
      f[j] = -1.0;
    } else {
      const int y = choice_index - 1;
      f[y] = 1.0;
      if (j != y) f[j] = -1.0;
    }

    std::vector<double> normal = lower_solve_from_upper(upper, f);
    const double norm_squared = dot(normal, normal);
    double offset = 0.0;
    for (int k = 0; k < m; ++k) offset += f[k] * mu[k];
    constraints.push_back({normal, offset, norm_squared});
  }
  return constraints;
}

bool feasible(
    const std::vector<double>& position,
    const std::vector<LinearConstraint>& constraints) {
  for (const LinearConstraint& constraint : constraints) {
    if (dot(constraint.normal, position) + constraint.offset <
        -FEASIBILITY_TOL) return false;
  }
  return true;
}

double normalize_time(double value) {
  const double period = 2.0 * PI;
  value = std::fmod(value, period);
  if (value < 0.0) value += period;
  return value;
}

double next_hit_time(
    const std::vector<double>& position,
    const std::vector<double>& velocity,
    const LinearConstraint& constraint) {
  const double a = dot(constraint.normal, velocity);
  const double b = dot(constraint.normal, position);
  const double amplitude = std::hypot(a, b);
  if (amplitude <= std::abs(constraint.offset)) {
    return std::numeric_limits<double>::infinity();
  }

  const double phi = std::atan2(-a, b);
  const double angle = std::acos(
    std::max(-1.0, std::min(1.0, -constraint.offset / amplitude))
  );
  const double candidates[2] = {
    normalize_time(angle - phi),
    normalize_time(-angle - phi)
  };

  double best = std::numeric_limits<double>::infinity();
  for (double candidate : candidates) {
    if (candidate <= HIT_EPS) continue;
    const double derivative = a * std::cos(candidate) - b * std::sin(candidate);
    if (derivative < -HIT_EPS && candidate < best) best = candidate;
  }
  return best;
}

void propagate(
    std::vector<double>& position,
    std::vector<double>& velocity,
    double time) {
  const double sine = std::sin(time);
  const double cosine = std::cos(time);
  for (std::size_t i = 0; i < position.size(); ++i) {
    const double old_position = position[i];
    const double old_velocity = velocity[i];
    position[i] = old_velocity * sine + old_position * cosine;
    velocity[i] = old_velocity * cosine - old_position * sine;
  }
}

bool transition(
    std::vector<double>& state,
    const std::vector<LinearConstraint>& constraints) {
  const int m = state.size();

  for (int retry = 0; retry < MAX_RETRIES; ++retry) {
    std::vector<double> position = state;
    std::vector<double> velocity(m);
    for (int j = 0; j < m; ++j) velocity[j] = R::rnorm(0.0, 1.0);
    double remaining = TRAJECTORY_TIME;
    int bounces = 0;

    while (remaining > HIT_EPS) {
      double hit_time = std::numeric_limits<double>::infinity();
      int hit_index = -1;
      for (std::size_t j = 0; j < constraints.size(); ++j) {
        const double candidate = next_hit_time(position, velocity, constraints[j]);
        if (candidate < hit_time) {
          hit_time = candidate;
          hit_index = static_cast<int>(j);
        }
      }

      if (hit_index < 0 || hit_time >= remaining) {
        propagate(position, velocity, remaining);
        remaining = 0.0;
        break;
      }

      propagate(position, velocity, hit_time);
      remaining -= hit_time;
      const LinearConstraint& constraint = constraints[hit_index];
      const double scale = 2.0 * dot(constraint.normal, velocity) /
        constraint.norm_squared;
      for (int j = 0; j < m; ++j) {
        velocity[j] -= scale * constraint.normal[j];
      }

      if (++bounces > MAX_BOUNCES) break;
    }

    if (bounces <= MAX_BOUNCES && feasible(position, constraints)) {
      state.swap(position);
      return true;
    }
  }
  return false;
}

void validate_inputs(
    const NumericVector& mu,
    const NumericMatrix& precision,
    int choice_index,
    int n_mc) {
  const int m = mu.size();
  if (m < 1) stop("mu must have positive length");
  if (precision.nrow() != m || precision.ncol() != m) {
    stop("Precision must be a square matrix with dimensions matching mu");
  }
  if (choice_index < 0 || choice_index > m) {
    stop("choice_index must be between 0 and length(mu)");
  }
  if (n_mc < 2) stop("n_mc must be at least 2");

  for (int i = 0; i < m; ++i) {
    if (!R_finite(mu[i])) stop("mu must contain only finite values");
    for (int j = 0; j < m; ++j) {
      if (!R_finite(precision(i, j))) {
        stop("Precision must contain only finite values");
      }
      const double scale = std::max(1.0, std::max(
        std::abs(precision(i, j)), std::abs(precision(j, i))
      ));
      if (std::abs(precision(i, j) - precision(j, i)) > 1e-10 * scale) {
        stop("Precision must be symmetric");
      }
    }
  }
}

} // namespace

List hmcmnp_moments(
    const NumericVector& mu,
    const NumericMatrix& Precision,
    int choice_index,
    int n_mc) {
  validate_inputs(mu, Precision, choice_index, n_mc);
  const int m = mu.size();
  const NumericMatrix upper = cholesky_upper(Precision);
  const std::vector<LinearConstraint> constraints =
    build_constraints(mu, upper, choice_index);

  std::vector<double> initial(m, -1.0);
  if (choice_index > 0) initial[choice_index - 1] = 1.0;
  std::vector<double> centered(m);
  for (int j = 0; j < m; ++j) centered[j] = initial[j] - mu[j];
  std::vector<double> state = upper_multiply(upper, centered);
  if (!feasible(state, constraints)) stop("Failed to construct a feasible HMC initial point");

  for (int i = 0; i < BURN_IN; ++i) {
    if (!transition(state, constraints)) stop("HMC transition failed during burn-in");
  }

  std::vector<double> mean(m, 0.0);
  NumericMatrix covariance_sum(m, m);
  for (int sample_index = 0; sample_index < n_mc; ++sample_index) {
    if (!transition(state, constraints)) stop("HMC transition failed while sampling");
    std::vector<double> sample = upper_solve(upper, state);
    for (int j = 0; j < m; ++j) sample[j] += mu[j];

    const double count = sample_index + 1.0;
    std::vector<double> delta(m);
    std::vector<double> delta_after(m);
    for (int j = 0; j < m; ++j) {
      delta[j] = sample[j] - mean[j];
      mean[j] += delta[j] / count;
      delta_after[j] = sample[j] - mean[j];
    }
    for (int r = 0; r < m; ++r) {
      for (int c = 0; c < m; ++c) {
        covariance_sum(r, c) += delta[r] * delta_after[c];
      }
    }
  }

  NumericVector mean_out(m);
  NumericMatrix covariance(m, m);
  for (int r = 0; r < m; ++r) {
    mean_out[r] = mean[r];
    for (int c = 0; c < m; ++c) {
      covariance(r, c) = covariance_sum(r, c) / (n_mc - 1.0);
    }
  }
  return List::create(Named("mu") = mean_out, Named("Sigma") = covariance);
}

//' Exact HMC moments for reference-identified multinomial probit
//'
//' @param mu Mean vector of the relative latent utilities.
//' @param Precision Precision matrix of the relative latent utilities.
//' @param choice_index Choice index: 0 for the base choice, otherwise 1..m.
//' @param n_mc Number of retained Monte Carlo samples after 30 burn-in draws.
//' @return A list containing Monte Carlo estimates `mu` and `Sigma`.
// [[Rcpp::export]]
List hmcmnp(
    NumericVector mu,
    NumericMatrix Precision,
    int choice_index,
    int n_mc) {
  RNGScope scope;
  return hmcmnp_moments(mu, Precision, choice_index, n_mc);
}
