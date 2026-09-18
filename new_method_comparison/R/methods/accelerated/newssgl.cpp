#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using arma::mat;
using arma::vec;
using Rcpp::List;
using Rcpp::Named;

namespace {

double positive_floor(double x) {
  return std::max(x, 1e-12);
}

double logistic(double x) {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

double draw_gig(Rcpp::Function& rgig, double lambda,
                double chi, double psi) {
  Rcpp::NumericVector value =
    rgig(1, lambda, positive_floor(chi), positive_floor(psi));
  return positive_floor(value[0]);
}

vec draw_standard_normal(int n) {
  vec z(n);
  for (int i = 0; i < n; ++i) z(i) = R::rnorm(0.0, 1.0);
  return z;
}

vec precision_mean(const mat& chol_upper, const vec& rhs) {
  return arma::solve(
    arma::trimatu(chol_upper),
    arma::solve(arma::trimatl(chol_upper.t()), rhs)
  );
}

vec draw_from_precision(const mat& chol_upper, const vec& mean,
                        double sigma2) {
  return mean + std::sqrt(sigma2) *
    arma::solve(arma::trimatu(chol_upper),
                draw_standard_normal(mean.n_elem));
}

vec make_prior_diagonal(const vec& omega, const vec& tau,
                        int p, int h) {
  vec prior(p + p * h);
  for (int j = 0; j < p; ++j) {
    prior(j) = 1.0 / positive_floor(omega(j));
    const double tau_precision = 1.0 / positive_floor(tau(j));
    for (int k = 0; k < h; ++k) {
      prior(p + j * h + k) = tau_precision;
    }
  }
  return prior;
}

vec gamma_probabilities(const vec& alpha, double sigma2,
                        double pi_gamma, double lambda0,
                        double lambda1, int p, int h) {
  vec probability(p);
  for (int j = 0; j < p; ++j) {
    const vec alpha_j = alpha.subvec(j * h, (j + 1) * h - 1);
    const double norm_scaled =
      std::sqrt(arma::dot(alpha_j, alpha_j) / sigma2);
    const double log_slab =
      std::log(positive_floor(pi_gamma)) +
      h * std::log(lambda1) - lambda1 * norm_scaled;
    const double log_spike =
      std::log(positive_floor(1.0 - pi_gamma)) +
      h * std::log(lambda0) - lambda0 * norm_scaled;
    probability(j) = logistic(log_slab - log_spike);
  }
  return probability;
}

}  // namespace

// Deterministic conditional-parameter helper used by the validation suite.
// [[Rcpp::export]]
List newssgl_conditionals_cpp(
    const vec& y, const mat& W, int p, int h,
    const vec& xi, const vec& omega, const vec& tau,
    const vec& gamma, double sigma2, double lambda_theta2,
    double pi_gamma, double a_sigma, double b_sigma,
    double a_theta, double b_theta, double a_gamma,
    double b_gamma, double lambda0, double lambda1) {
  const int n = W.n_rows;
  const int q = W.n_cols;
  const vec prior = make_prior_diagonal(omega, tau, p, h);
  mat precision = W.t() * W;
  precision.diag() += prior;
  const vec rhs = W.t() * y;
  const mat chol_upper = arma::chol(precision);
  const vec mean = precision_mean(chol_upper, rhs);

  const vec theta = xi.subvec(0, p - 1);
  const vec alpha = xi.subvec(p, q - 1);
  const vec residual = y - W * xi;
  double penalty = 0.0;
  for (int j = 0; j < p; ++j) {
    penalty += theta(j) * theta(j) / positive_floor(omega(j));
    const vec alpha_j = alpha.subvec(j * h, (j + 1) * h - 1);
    penalty += arma::dot(alpha_j, alpha_j) / positive_floor(tau(j));
  }

  const vec gamma_prob = gamma_probabilities(
    alpha, sigma2, pi_gamma, lambda0, lambda1, p, h
  );

  return List::create(
    Named("prior_diagonal") = prior,
    Named("precision") = precision,
    Named("xi_mean") = mean,
    Named("xi_covariance") = sigma2 * arma::inv_sympd(precision),
    Named("sigma_shape") = a_sigma + 0.5 * (n + q),
    Named("sigma_rate") =
      b_sigma + 0.5 * (arma::dot(residual, residual) + penalty),
    Named("omega_lambda") = vec(p, arma::fill::value(0.5)),
    Named("omega_chi") = arma::square(theta) / sigma2,
    Named("omega_psi") = vec(p, arma::fill::value(lambda_theta2)),
    Named("lambda_theta_shape") = a_theta + p,
    Named("lambda_theta_rate") = b_theta + 0.5 * arma::sum(omega),
    Named("pi_shape1") = a_gamma + arma::sum(gamma),
    Named("pi_shape2") = b_gamma + p - arma::sum(gamma),
    Named("gamma_probability") = gamma_prob,
    Named("tau_lambda") = vec(p, arma::fill::value(0.5)),
    Named("tau_chi") = [&]() {
      vec out(p);
      for (int j = 0; j < p; ++j) {
        const vec alpha_j = alpha.subvec(j * h, (j + 1) * h - 1);
        out(j) = arma::dot(alpha_j, alpha_j) / sigma2;
      }
      return out;
    }(),
    Named("tau_psi_if_spike") =
      vec(p, arma::fill::value(lambda0 * lambda0)),
    Named("tau_psi_if_slab") =
      vec(p, arma::fill::value(lambda1 * lambda1))
  );
}

// [[Rcpp::export]]
List newssgl_gibbs_cpp(
    const vec& y, const mat& W, int p, int h,
    int n_iter, int burn_in,
    double lambda0, double lambda1,
    double a_sigma, double b_sigma,
    double a_theta, double b_theta,
    double a_gamma, double b_gamma,
    vec xi, vec omega, vec gamma, vec tau,
    double sigma2, double lambda_theta2, double pi_gamma) {
  Rcpp::RNGScope scope;
  const int n = W.n_rows;
  const int q = W.n_cols;
  const int keep = n_iter - burn_in;
  if (keep <= 0) Rcpp::stop("n_iter must exceed burn_in.");
  if (W.n_cols != static_cast<unsigned>(p + p * h)) {
    Rcpp::stop("W does not match p + p*h.");
  }

  const mat WtW = W.t() * W;
  const vec Wty = W.t() * y;
  mat precision = WtW;
  mat chol_upper(q, q);
  vec theta(p);
  vec alpha(p * h);

  mat theta_draws(p, keep);
  mat alpha_draws(p * h, keep);
  mat gamma_draws(p, keep);
  mat gamma_prob_draws(p, keep);
  mat omega_draws(p, keep);
  mat tau_draws(p, keep);
  vec sigma2_draws(keep);
  vec lambda_theta2_draws(keep);
  vec pi_gamma_draws(keep);

  Rcpp::Environment gigrvg =
    Rcpp::Environment::namespace_env("GIGrvg");
  Rcpp::Function rgig = gigrvg["rgig"];

  for (int iter = 0; iter < n_iter; ++iter) {
    precision = WtW;
    precision.diag() += make_prior_diagonal(omega, tau, p, h);
    if (!arma::chol(chol_upper, precision)) {
      Rcpp::stop("Precision Cholesky factorization failed.");
    }
    const vec xi_mean = precision_mean(chol_upper, Wty);
    xi = draw_from_precision(chol_upper, xi_mean, sigma2);
    theta = xi.subvec(0, p - 1);
    alpha = xi.subvec(p, q - 1);

    const vec residual = y - W * xi;
    double penalty = 0.0;
    for (int j = 0; j < p; ++j) {
      penalty += theta(j) * theta(j) / positive_floor(omega(j));
      const vec alpha_j = alpha.subvec(j * h, (j + 1) * h - 1);
      penalty += arma::dot(alpha_j, alpha_j) / positive_floor(tau(j));
    }
    const double sigma_shape = a_sigma + 0.5 * (n + q);
    const double sigma_rate =
      b_sigma + 0.5 * (arma::dot(residual, residual) + penalty);
    sigma2 = 1.0 / R::rgamma(sigma_shape, 1.0 / sigma_rate);

    for (int j = 0; j < p; ++j) {
      omega(j) = draw_gig(
        rgig, 0.5, theta(j) * theta(j) / sigma2, lambda_theta2
      );
    }
    lambda_theta2 = R::rgamma(
      a_theta + p, 1.0 / (b_theta + 0.5 * arma::sum(omega))
    );
    pi_gamma = R::rbeta(
      a_gamma + arma::sum(gamma),
      b_gamma + p - arma::sum(gamma)
    );

    const vec gamma_prob = gamma_probabilities(
      alpha, sigma2, pi_gamma, lambda0, lambda1, p, h
    );
    for (int j = 0; j < p; ++j) {
      gamma(j) = R::rbinom(1.0, gamma_prob(j));
      const double lambda = gamma(j) > 0.5 ? lambda1 : lambda0;
      const vec alpha_j = alpha.subvec(j * h, (j + 1) * h - 1);
      tau(j) = draw_gig(
        rgig, 0.5, arma::dot(alpha_j, alpha_j) / sigma2,
        lambda * lambda
      );
    }

    if (iter >= burn_in) {
      const int k = iter - burn_in;
      theta_draws.col(k) = theta;
      alpha_draws.col(k) = alpha;
      gamma_draws.col(k) = gamma;
      gamma_prob_draws.col(k) = gamma_prob;
      omega_draws.col(k) = omega;
      tau_draws.col(k) = tau;
      sigma2_draws(k) = sigma2;
      lambda_theta2_draws(k) = lambda_theta2;
      pi_gamma_draws(k) = pi_gamma;
    }
  }

  return List::create(
    Named("theta") = theta_draws,
    Named("alpha") = alpha_draws,
    Named("gamma") = gamma_draws,
    Named("gamma_prob") = gamma_prob_draws,
    Named("omega") = omega_draws,
    Named("tau") = tau_draws,
    Named("sigma2") = sigma2_draws,
    Named("lambda_theta2") = lambda_theta2_draws,
    Named("pi_gamma") = pi_gamma_draws
  );
}
