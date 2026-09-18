#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double dmvnorm0_log_var_cpp(const arma::vec& x, double variance) {
  int n = x.n_elem;
  double log_det = n * std::log(variance);
  double quad_form = sum(square(x)) / variance;
  return -0.5 * (n * std::log(2 * M_PI) + log_det + quad_form);
}

// [[Rcpp::export]]
List ssgl_intercept_cpp(arma::vec y, arma::mat X, arma::mat B,
                        int n_iter, int burn_in, double a, double b,
                        double a_sigma, double b_sigma,
                        double lambda0, double lambda1,
                        double beta0_init = NA_REAL,
                        double zeta0_init = 0.1,
                        double zeta1_init = 50) {
  int n = X.n_rows;
  int p = X.n_cols;
  int n_basis = B.n_cols;
  int eta_dim = p * n_basis;
  int q = 1 + eta_dim;

  arma::mat Z = zeros<mat>(n, eta_dim);
  for (int j = 0; j < p; j++) {
    Z.cols(j * n_basis, (j + 1) * n_basis - 1) = B.each_col() % X.col(j);
  }
  arma::mat W(n, q);
  W.col(0).ones();
  W.cols(1, q - 1) = Z;

  double sigma2 = var(y);
  if (!arma::is_finite(sigma2) || sigma2 <= 0) sigma2 = 1.0;
  double theta = 0.5;
  double zeta0 = zeta0_init;
  double zeta1 = zeta1_init;
  double beta0 = R_IsNA(beta0_init) ? mean(y) : beta0_init;

  arma::vec gamma = zeros<vec>(p);
  for (int j = 0; j < p; j++) gamma(j) = R::rbinom(1, 0.5);
  arma::vec xi = zeros<vec>(q);
  xi(0) = beta0;
  arma::vec eta = zeros<vec>(eta_dim);

  int n_samples = n_iter - burn_in;
  arma::vec beta0_samples = zeros<vec>(n_samples);
  arma::vec zeta0_samples = zeros<vec>(n_samples);
  arma::vec zeta1_samples = zeros<vec>(n_samples);
  arma::vec sigma2_samples = zeros<vec>(n_samples);
  arma::vec theta_samples = zeros<vec>(n_samples);
  arma::mat gamma_samples = zeros<mat>(p, n_samples);
  arma::mat gamma_prob_samples = zeros<mat>(p, n_samples);
  arma::mat eta_samples = zeros<mat>(eta_dim, n_samples);

  Environment gigrvg("package:GIGrvg");
  Function rgig = gigrvg["rgig"];

  arma::mat WtW = W.t() * W;
  arma::vec Wty = W.t() * y;

  for (int iter = 0; iter < n_iter; iter++) {
    arma::mat M_inv = zeros<mat>(q, q);
    for (int j = 0; j < p; j++) {
      double zeta_j = (gamma(j) == 1) ? zeta1 : zeta0;
      M_inv.submat(1 + j * n_basis, 1 + j * n_basis,
                   1 + (j + 1) * n_basis - 1,
                   1 + (j + 1) * n_basis - 1) =
        eye<mat>(n_basis, n_basis) / zeta_j;
    }

    arma::mat K = inv_sympd(M_inv + WtW / sigma2);
    arma::vec mean_xi = K * (Wty / sigma2);
    xi = mvnrnd(mean_xi, K);
    beta0 = xi(0);
    eta = xi.subvec(1, q - 1);

    arma::vec residuals = y - W * xi;
    double shape_sigma = n / 2.0 + a_sigma;
    double rate_sigma = dot(residuals, residuals) / 2.0 + b_sigma;
    sigma2 = 1.0 / R::rgamma(shape_sigma, 1.0 / rate_sigma);

    theta = R::rbeta(a + sum(gamma), b + p - sum(gamma));

    arma::vec eta_squares = zeros<vec>(p);
    for (int j = 0; j < p; j++) {
      arma::vec eta_j = eta.subvec(j * n_basis, (j + 1) * n_basis - 1);
      eta_squares(j) = sum(square(eta_j));
    }

    if (any(gamma == 1)) {
      double sum_eta_squared = 0;
      for (int j = 0; j < p; j++) {
        if (gamma(j) == 1) sum_eta_squared += eta_squares(j);
      }
      int n_selected = sum(gamma);
      double lambda = (n_basis + 1.0 - n_selected * n_basis) / 2.0;
      NumericVector gig = rgig(1, lambda, std::max(sum_eta_squared, 1e-12),
                               lambda1 * lambda1);
      zeta1 = gig[0];
    }

    if (any(gamma == 0)) {
      double sum_eta_squared = 0;
      for (int j = 0; j < p; j++) {
        if (gamma(j) == 0) sum_eta_squared += eta_squares(j);
      }
      int n_unselected = p - sum(gamma);
      double lambda = (n_basis + 1.0 - n_unselected * n_basis) / 2.0;
      NumericVector gig = rgig(1, lambda, std::max(sum_eta_squared, 1e-12),
                               lambda0 * lambda0);
      zeta0 = gig[0];
    }

    for (int j = 0; j < p; j++) {
      arma::vec eta_j = eta.subvec(j * n_basis, (j + 1) * n_basis - 1);
      double log_prob_slab = std::log(std::max(theta, 1e-12)) +
        dmvnorm0_log_var_cpp(eta_j, zeta1);
      double log_prob_spike = std::log(std::max(1.0 - theta, 1e-12)) +
        dmvnorm0_log_var_cpp(eta_j, zeta0);
      double log_odds = log_prob_slab - log_prob_spike;
      double prob_slab = 1.0 / (1.0 + std::exp(-log_odds));
      if (iter >= burn_in) {
        int idx = iter - burn_in;
        gamma_prob_samples(j, idx) = prob_slab;
      }
      gamma(j) = R::rbinom(1, prob_slab);
    }

    if (iter >= burn_in) {
      int idx = iter - burn_in;
      beta0_samples(idx) = beta0;
      zeta0_samples(idx) = zeta0;
      zeta1_samples(idx) = zeta1;
      sigma2_samples(idx) = sigma2;
      theta_samples(idx) = theta;
      gamma_samples.col(idx) = gamma;
      eta_samples.col(idx) = eta;
    }
  }

  return List::create(
    Named("beta0") = beta0_samples,
    Named("zeta0") = zeta0_samples,
    Named("zeta1") = zeta1_samples,
    Named("sigma2") = sigma2_samples,
    Named("theta") = theta_samples,
    Named("gamma") = gamma_samples,
    Named("gamma_prob") = gamma_prob_samples,
    Named("eta") = eta_samples
  );
}
