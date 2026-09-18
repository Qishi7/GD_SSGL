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

vec make_prior_diagonal_intercept(const vec& omega, const vec& tau,
                                  int p, int h) {
  vec prior(1 + p + p * h);
  prior(0) = 0.0;
  for (int j = 0; j < p; ++j) {
    prior(1 + j) = 1.0 / positive_floor(omega(j));
    const double tau_precision = 1.0 / positive_floor(tau(j));
    for (int k = 0; k < h; ++k) {
      prior(1 + p + j * h + k) = tau_precision;
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

// [[Rcpp::export]]
List newssgl_intercept_gibbs_cpp(
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
  if (W.n_cols != static_cast<unsigned>(1 + p + p * h)) {
    Rcpp::stop("W does not match 1 + p + p*h.");
  }
  if (xi.n_elem != static_cast<unsigned>(q)) {
    Rcpp::stop("xi does not match W columns.");
  }

  const mat WtW = W.t() * W;
  const vec Wty = W.t() * y;
  mat precision = WtW;
  mat chol_upper(q, q);
  vec theta(p);
  vec alpha(p * h);

  vec beta0_draws(keep);
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
    precision.diag() += make_prior_diagonal_intercept(omega, tau, p, h);
    if (!arma::chol(chol_upper, precision)) {
      Rcpp::stop("Precision Cholesky factorization failed.");
    }
    const vec xi_mean = precision_mean(chol_upper, Wty);
    xi = draw_from_precision(chol_upper, xi_mean, sigma2);
    const double beta0 = xi(0);
    theta = xi.subvec(1, p);
    alpha = xi.subvec(1 + p, q - 1);

    const vec residual = y - W * xi;
    double penalty = 0.0;
    for (int j = 0; j < p; ++j) {
      penalty += theta(j) * theta(j) / positive_floor(omega(j));
      const vec alpha_j = alpha.subvec(j * h, (j + 1) * h - 1);
      penalty += arma::dot(alpha_j, alpha_j) / positive_floor(tau(j));
    }
    const double sigma_shape = a_sigma + 0.5 * (n + q - 1);
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
      beta0_draws(k) = beta0;
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
    Named("beta0") = beta0_draws,
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

// [[Rcpp::export]]
List newssgl_intercept_collapsed_gibbs_cpp(
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
  if (xi.n_elem != static_cast<unsigned>(q)) {
    Rcpp::stop("xi does not match W columns.");
  }

  const vec y_centered = y - arma::mean(y);
  mat W_centered = W;
  for (int j = 0; j < q; ++j) {
    W_centered.col(j) -= arma::mean(W.col(j));
  }
  const mat WtW = W_centered.t() * W_centered;
  const vec Wty = W_centered.t() * y_centered;
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
  vec beta0_cond_mean_draws(keep);

  Rcpp::Environment gigrvg =
    Rcpp::Environment::namespace_env("GIGrvg");
  Rcpp::Function rgig = gigrvg["rgig"];

  for (int iter = 0; iter < n_iter; ++iter) {
    precision = WtW;
    precision.diag() += make_prior_diagonal_intercept(omega, tau, p, h).subvec(1, q);
    if (!arma::chol(chol_upper, precision)) {
      Rcpp::stop("Precision Cholesky factorization failed.");
    }
    const vec xi_mean = precision_mean(chol_upper, Wty);
    xi = draw_from_precision(chol_upper, xi_mean, sigma2);
    theta = xi.subvec(0, p - 1);
    alpha = xi.subvec(p, q - 1);

    const vec residual_centered = y_centered - W_centered * xi;
    double penalty = 0.0;
    for (int j = 0; j < p; ++j) {
      penalty += theta(j) * theta(j) / positive_floor(omega(j));
      const vec alpha_j = alpha.subvec(j * h, (j + 1) * h - 1);
      penalty += arma::dot(alpha_j, alpha_j) / positive_floor(tau(j));
    }
    const double sigma_shape = a_sigma + 0.5 * ((n - 1) + q);
    const double sigma_rate =
      b_sigma + 0.5 * (arma::dot(residual_centered, residual_centered) + penalty);
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
      beta0_cond_mean_draws(k) = arma::mean(y - W * xi);
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
    Named("pi_gamma") = pi_gamma_draws,
    Named("beta0_cond_mean") = beta0_cond_mean_draws
  );
}
