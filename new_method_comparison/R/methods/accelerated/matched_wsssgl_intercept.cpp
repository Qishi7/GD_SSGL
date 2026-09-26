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

double draw_gig(Rcpp::Function& rgig, double lambda, double chi, double psi) {
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

vec draw_from_precision(const mat& chol_upper, const vec& mean, double sigma2) {
  return mean + std::sqrt(sigma2) *
    arma::solve(arma::trimatu(chol_upper), draw_standard_normal(mean.n_elem));
}

vec make_prior_diagonal(double beta0_precision, const vec& tau, int p, int H) {
  vec prior(1 + p * H);
  prior(0) = beta0_precision;
  for (int j = 0; j < p; ++j) {
    const double tau_precision = 1.0 / positive_floor(tau(j));
    for (int k = 0; k < H; ++k) {
      prior(1 + j * H + k) = tau_precision;
    }
  }
  return prior;
}

vec gamma_probabilities(const vec& eta, double sigma2, double pi_gamma,
                        double lambda0, double lambda1, int p, int H) {
  vec probability(p);
  for (int j = 0; j < p; ++j) {
    const vec eta_j = eta.subvec(j * H, (j + 1) * H - 1);
    const double norm_scaled =
      std::sqrt(arma::dot(eta_j, eta_j) / positive_floor(sigma2));
    const double log_slab =
      std::log(positive_floor(pi_gamma)) +
      H * std::log(lambda1) - lambda1 * norm_scaled;
    const double log_spike =
      std::log(positive_floor(1.0 - pi_gamma)) +
      H * std::log(lambda0) - lambda0 * norm_scaled;
    probability(j) = logistic(log_slab - log_spike);
  }
  return probability;
}

}  // namespace

// [[Rcpp::export]]
List matched_wsssgl_intercept_gibbs_cpp(
    const vec& y, const mat& W, int p, int H,
    int n_iter, int burn_in,
    double lambda0, double lambda1,
    double a_sigma, double b_sigma,
    double a_gamma, double b_gamma,
    vec xi, vec gamma, vec tau,
    double sigma2, double pi_gamma) {
  Rcpp::RNGScope scope;
  const int n = W.n_rows;
  const int q = W.n_cols;
  const int keep = n_iter - burn_in;
  if (keep <= 0) Rcpp::stop("n_iter must exceed burn_in.");
  if (W.n_cols != static_cast<unsigned>(1 + p * H)) {
    Rcpp::stop("W does not match 1 + p*H.");
  }
  if (xi.n_elem != static_cast<unsigned>(q)) {
    Rcpp::stop("xi does not match W columns.");
  }
  if (gamma.n_elem != static_cast<unsigned>(p) ||
      tau.n_elem != static_cast<unsigned>(p)) {
    Rcpp::stop("gamma and tau must have length p.");
  }

  const mat WtW = W.t() * W;
  const vec Wty = W.t() * y;
  mat precision = WtW;
  mat chol_upper(q, q);
  vec eta(p * H);

  vec beta0_draws(keep);
  mat eta_draws(p * H, keep);
  mat gamma_draws(p, keep);
  mat gamma_prob_draws(p, keep);
  mat tau_draws(p, keep);
  vec sigma2_draws(keep);
  vec pi_gamma_draws(keep);

  Rcpp::Environment gigrvg =
    Rcpp::Environment::namespace_env("GIGrvg");
  Rcpp::Function rgig = gigrvg["rgig"];

  for (int iter = 0; iter < n_iter; ++iter) {
    precision = WtW;
    precision.diag() += make_prior_diagonal(0.0, tau, p, H);
    if (!arma::chol(chol_upper, precision)) {
      Rcpp::stop("Precision Cholesky factorization failed.");
    }
    const vec xi_mean = precision_mean(chol_upper, Wty);
    xi = draw_from_precision(chol_upper, xi_mean, sigma2);
    const double beta0 = xi(0);
    eta = xi.subvec(1, q - 1);

    const vec residual = y - W * xi;
    double penalty = 0.0;
    for (int j = 0; j < p; ++j) {
      const vec eta_j = eta.subvec(j * H, (j + 1) * H - 1);
      penalty += arma::dot(eta_j, eta_j) / positive_floor(tau(j));
    }
    const double sigma_shape = a_sigma + 0.5 * (n + p * H);
    const double sigma_rate =
      b_sigma + 0.5 * (arma::dot(residual, residual) + penalty);
    sigma2 = 1.0 / R::rgamma(sigma_shape, 1.0 / sigma_rate);

    pi_gamma = R::rbeta(
      a_gamma + arma::sum(gamma),
      b_gamma + p - arma::sum(gamma)
    );

    const vec gamma_prob = gamma_probabilities(
      eta, sigma2, pi_gamma, lambda0, lambda1, p, H
    );
    for (int j = 0; j < p; ++j) {
      gamma(j) = R::rbinom(1.0, gamma_prob(j));
      const double lambda = gamma(j) > 0.5 ? lambda1 : lambda0;
      const vec eta_j = eta.subvec(j * H, (j + 1) * H - 1);
      tau(j) = draw_gig(
        rgig, 0.5, arma::dot(eta_j, eta_j) / positive_floor(sigma2),
        lambda * lambda
      );
    }

    if (iter >= burn_in) {
      const int k = iter - burn_in;
      beta0_draws(k) = beta0;
      eta_draws.col(k) = eta;
      gamma_draws.col(k) = gamma;
      gamma_prob_draws.col(k) = gamma_prob;
      tau_draws.col(k) = tau;
      sigma2_draws(k) = sigma2;
      pi_gamma_draws(k) = pi_gamma;
    }
  }

  return List::create(
    Named("beta0") = beta0_draws,
    Named("eta") = eta_draws,
    Named("gamma") = gamma_draws,
    Named("gamma_prob") = gamma_prob_draws,
    Named("tau") = tau_draws,
    Named("sigma2") = sigma2_draws,
    Named("pi_gamma") = pi_gamma_draws
  );
}
