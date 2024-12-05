// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

struct alpha_draw
{
    vec alpha;
    vec mean_post;
    mat cov_post;
};

struct Sigma_draw
{
    mat Sigma;
    mat Sigma_post;
    int nu_post;
};

double loglik_normal_iwish(mat Sigma, vec alpha, vec y, mat X, int m, int T)
{
    mat Sigma_inv = inv_sympd(Sigma);
    vec alpha_hat = inv(kron(Sigma_inv, X.t() * X)) * kron(Sigma_inv, X).t() * y;
    mat eye_m(m, m, fill::eye);
    mat eye_T(T, T, fill::eye);
    mat X_cal = kron(eye_m, X);
    mat mat1 = (alpha - alpha_hat).t() * X_cal.t();
    mat mat2 = kron(Sigma_inv, eye_T) * X_cal * (alpha - alpha_hat);
    double normal = - T * det(Sigma) / 2 - as_scalar(mat1 * mat2) / 2;
    double iwish = - trace((y - X_cal * alpha_hat) * (y - X_cal * alpha_hat).t() * kron(Sigma_inv, eye_T)) / 2;
    double loglik = normal + iwish;
    return loglik;
}

alpha_draw updt_alpha(vec alpha, vec alpha_prior, mat sigma_prior_inv, mat Sigma, mat y, mat X, int m, int T)
{
    mat Sigma_inv = inv_sympd(Sigma);
    mat eye_m(m, m, fill::eye);
    mat eye_T(T, T, fill::eye);
    mat X_cal = kron(eye_m, X);
    mat Sigma_inv_cal = kron(Sigma_inv, eye_T);
    mat Sigma_post = inv(sigma_prior_inv + X_cal.t() * Sigma_inv_cal * X_cal);
    vec alpha_post = Sigma_post * (sigma_prior_inv * alpha_prior + X_cal.t() * Sigma_inv_cal * y);
    vec alpha_updt = mvnrnd(alpha_post, Sigma_post);
    alpha_draw out;
    out.alpha = alpha_updt;
    out.mean_post = alpha_post;
    out.cov_post = Sigma_post;
    return out;
}

Sigma_draw updt_Sigma(mat Sigma, mat Sigma_prior, int nu_prior, vec alpha, mat Y, mat X, int m, int T)
{
    mat beta = reshape(alpha, alpha.n_elem / m, m);
    mat residual = Y - X * beta;
    mat Sigma_post = Sigma_prior + residual.t() * residual;
    int nu_post = nu_prior + T;
    mat Sigma_updt = iwishrnd(Sigma_post, nu_post);
    Sigma_draw out;
    out.Sigma = Sigma_updt;
    out.Sigma_post = Sigma_post;
    out.nu_post = nu_post;
    return out;
}

// [[Rcpp::export]]
List gibbs_sampler(arma::mat Y, arma::mat X, List priors, int burn_in = 100, int draw = 1000, int thin = 1, int post_save = 500, int post_method = 0)
{
    // priors
    arma::vec alpha_bar = as<arma::vec>(priors[0]);
    arma::mat Sigma_alpha = as<arma::mat>(priors[1]);
    arma::mat Sigma_alpha_inv = arma::inv_sympd(Sigma_alpha);
    arma::mat Sigma_sigma = as<arma::mat>(priors[2]);
    int nu = as<int>(priors[3]);
    // data
    arma::vec y = vectorise(Y);
    // dimensions
    int m = Y.n_cols;
    int T_est = Y.n_rows;
    int K = alpha_bar.n_elem;
    // save draws
    arma::mat alpha_draws(K, draw, fill::zeros);
    arma::cube Sigma_draws(m, m, draw, fill::zeros);
    arma::vec loglik_draws(draw, fill::zeros);
    // save posterior draws
    arma::mat alpha_post(K, post_save, fill::zeros);
    arma::cube Sigma_alpha_post(K, K, post_save, fill::zeros);
    arma::cube Sigma_sigma_post(m, m, post_save, fill::zeros);
    arma::vec nu_post(post_save, fill::zeros);
    // save posterior out
    arma::vec alpha_mean_post(K, fill::zeros);
    arma::mat alpha_cov_post(K, K, fill::zeros);
    arma::mat Sigma_mean_post(m, m, fill::zeros);
    int nu_post1 = 0;
    // save updt function output
    alpha_draw alpha_out;
    Sigma_draw Sigma_out;

    // first draw and supervisor
    arma::vec alpha = arma::mvnrnd(alpha_bar, Sigma_alpha);
    arma::mat Sigma = arma::iwishrnd(Sigma_sigma, nu);
    int report = draw / 10;
    // burn in
    for (int i = 0; i < burn_in; i++)
    {
        alpha_out = updt_alpha(alpha, alpha_bar, Sigma_alpha_inv, Sigma, y, X, m, T_est);
        alpha = alpha_out.alpha;
        Sigma_out = updt_Sigma(Sigma, Sigma_sigma, nu, alpha, Y, X, m, T_est);
        Sigma = Sigma_out.Sigma;
    }
    for (int i = 0; i < draw; i++)
    {
        for (int j = 0; j < thin; j++)
        {
            alpha_out = updt_alpha(alpha, alpha_bar, Sigma_alpha_inv, Sigma, y, X, m, T_est);
            alpha = alpha_out.alpha;
            Sigma_out = updt_Sigma(Sigma, Sigma_sigma, nu, alpha, Y, X, m, T_est);
            Sigma = Sigma_out.Sigma;
        }
        alpha_draws.col(i) = alpha;
        Sigma_draws.slice(i) = Sigma;
        loglik_draws(i) = loglik_normal_iwish(Sigma, alpha, y, X, m, T_est);
        if (i % report == 0) {
            Rcout << i << " samples drawn from Gibbs sampler.\n";
        }
        if (i >= (draw - post_save))
        {
            int k = i - post_save;
            alpha_post.col(k) = alpha_out.mean_post;
            Sigma_alpha_post.slice(k) = alpha_out.cov_post;
            Sigma_sigma_post.slice(k) = Sigma_out.Sigma_post;
            nu_post(k) = Sigma_out.nu_post;
        }
    }
    if (post_method == 0)
    {
        alpha_mean_post = alpha_post.col(post_save - 1);
        alpha_cov_post = Sigma_alpha_post.slice(post_save - 1);
        Sigma_mean_post = Sigma_sigma_post.slice(post_save - 1);
        nu_post1 = nu_post(post_save - 1);
    }
    else
    {
        alpha_mean_post = mean(alpha_post, 1);
        alpha_cov_post = mean(Sigma_alpha_post, 2);
        Sigma_mean_post = mean(Sigma_sigma_post, 2);
        nu_post1 = mean(nu_post);
    }
    List out = List::create(
        _["alpha"] = alpha_draws,
        _["Sigma"] = Sigma_draws,
        _["loglik"] = loglik_draws,
        _["alpha_mean_post"] = alpha_mean_post,
        _["alpha_cov_post"] = alpha_cov_post,
        _["Sigma_mean_post"] = Sigma_mean_post,
        _["nu_post"] = nu_post1

    );
    return out;
}
