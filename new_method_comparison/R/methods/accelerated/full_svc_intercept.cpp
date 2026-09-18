#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List full_svc_intercept_cpp(arma::vec y, arma::mat W, int p, int h,
                            int n_iter, int burn_in,
                            double a_sigma, double b_sigma,
                            double a_theta, double b_theta,
                            double kappa2_alpha,
                            arma::vec xi, arma::vec omega,
                            double sigma2, double lambda2) {
  int n = W.n_rows;
  int q = W.n_cols;
  int keep = n_iter - burn_in;
  int alpha_dim = p * h;
  int q_penalized = q - 1;

  arma::mat WtW = W.t() * W;
  arma::vec Wty = W.t() * y;

  arma::vec beta0_draws = zeros<vec>(keep);
  arma::mat theta_draws = zeros<mat>(p, keep);
  arma::mat alpha_draws = zeros<mat>(alpha_dim, keep);
  arma::vec sigma2_draws = zeros<vec>(keep);
  arma::vec lambda2_draws = zeros<vec>(keep);
  arma::mat omega_draws = zeros<mat>(p, keep);

  Environment gigrvg("package:GIGrvg");
  Function rgig = gigrvg["rgig"];

  for (int iter = 0; iter < n_iter; iter++) {
    arma::vec prior(q, fill::zeros);
    for (int j = 0; j < p; j++) {
      prior(1 + j) = 1.0 / std::max(omega(j), 1e-10);
    }
    for (int k = 0; k < alpha_dim; k++) {
      prior(1 + p + k) = 1.0 / kappa2_alpha;
    }

    arma::mat precision = WtW + diagmat(prior);
    arma::mat K = inv_sympd(precision);
    arma::vec mean_xi = K * Wty;
    xi = mvnrnd(mean_xi, sigma2 * K);

    arma::vec theta = xi.subvec(1, p);
    arma::vec alpha = xi.subvec(1 + p, q - 1);
    arma::vec residual = y - W * xi;

    double penalty = 0.0;
    for (int j = 0; j < p; j++) {
      penalty += theta(j) * theta(j) / std::max(omega(j), 1e-10);
    }
    penalty += dot(alpha, alpha) / kappa2_alpha;

    sigma2 = 1.0 / R::rgamma(
      a_sigma + 0.5 * (n + q_penalized),
      1.0 / (b_sigma + 0.5 * (dot(residual, residual) + penalty))
    );

    for (int j = 0; j < p; j++) {
      NumericVector gig = rgig(1, 0.5,
                               std::max(theta(j) * theta(j) / sigma2, 1e-12),
                               lambda2);
      omega(j) = gig[0];
    }
    lambda2 = R::rgamma(a_theta + p,
                        1.0 / (b_theta + 0.5 * sum(omega)));

    if (iter >= burn_in) {
      int k = iter - burn_in;
      beta0_draws(k) = xi(0);
      theta_draws.col(k) = theta;
      alpha_draws.col(k) = alpha;
      sigma2_draws(k) = sigma2;
      lambda2_draws(k) = lambda2;
      omega_draws.col(k) = omega;
    }
  }

  return List::create(
    Named("beta0") = beta0_draws,
    Named("theta") = theta_draws,
    Named("alpha") = alpha_draws,
    Named("sigma2") = sigma2_draws,
    Named("lambda2") = lambda2_draws,
    Named("omega") = omega_draws
  );
}
