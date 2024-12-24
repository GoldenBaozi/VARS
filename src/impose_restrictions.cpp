// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <chrono>
#include "misc.h"
#include "VARtools.h"
#include "check_restrictions.h"

using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace chrono;

// [[Rcpp::export]]
List impose_SR_and_EBR(List alpha_post, List Sigma_post, List restrictions, arma::mat Y, arma::mat X, int draw, int plag, int hor)
{
    int nvar = Y.n_cols;
    // int T_est = Y.n_rows;
    // save draws
    arma::cube B_draw(nvar, nvar, draw);
    arma::cube beta_draw(nvar * plag + 1, nvar, draw);
    arma::cube Sigma_draw(nvar, nvar, draw);
    // unpack priors
    arma::vec alpha_post_mean = as<arma::vec>(alpha_post[0]);
    arma::mat alpha_post_cov = as<arma::mat>(alpha_post[1]);
    arma::mat Sigma_post_mean = as<arma::mat>(Sigma_post[0]);
    int nu_post = as<int>(Sigma_post[1]);
    // flag and supervisor
    int flag = 0;
    int report = draw / 10;
    int try_per_beta_sigma = 100;
    auto start_time = std::chrono::high_resolution_clock::now();
    // loop on draw
    while (flag < draw)
    {
        arma::mat Sigma_1draw = arma::iwishrnd(Sigma_post_mean, nu_post);
        arma::mat N0 = arma::inv(X.t() * X);
        arma::mat NT = arma::kron(Sigma_1draw, N0);
        arma::mat cov = arma::chol(NT, "lower");
        arma::vec alpha_1draw = alpha_post_mean + cov * randn(alpha_post_mean.n_elem, 1);
        arma::mat beta_1draw = arma::reshape(alpha_1draw, nvar * plag + 1, nvar);
        arma::mat u_1draw = Y - X * beta_1draw;
        arma::mat P = arma::chol(Sigma_1draw, "lower");
        double maxroot = max_root(beta_1draw);
        if (maxroot <= 1)
        {
            for (int i = 0; i < try_per_beta_sigma; i++)
            {
                arma::mat X = arma::randn(nvar, nvar);
                arma::mat Q, R;
                arma::qr(Q, R, X);
                Q = Q * arma::diagmat(sign(R.diag()));
                arma::mat B_1draw = P * Q;
                arma::mat B_inv_1draw = arma::inv(arma::trans(B_1draw));
                // mat eps_1draw = u_1draw * B_inv_1draw;
                Psi_IRF tmp = IRF_compute_internal(beta_1draw, B_1draw, hor, nvar, plag);
                arma::cube IRF_1draw = tmp.IRF;
                // check sign restrictions, to decide save or not
                bool save_me = check_SR_and_EBR(restrictions, IRF_1draw);
                if (save_me)
                {
                    B_draw.slice(flag) = B_1draw;
                    beta_draw.slice(flag) = beta_1draw;
                    Sigma_draw.slice(flag) = Sigma_1draw;
                    // Rcout << IRF_1draw.n_slices << "\n";
                    // IRF_draw.slice(flag) = flatten_cube(IRF_1draw);
                    // weights(flag) = compute_importance_weight(restrictions, IRF_1draw, M, T_est, nvar);
                    // Rcout << weights(flag) <<"\n";
                    flag++;
                    if (flag % report == 0)
                    {
                        auto current_time = std::chrono::high_resolution_clock::now();
                        auto elapsed_time = std::chrono::duration_cast<seconds>(current_time - start_time).count();
                        Rcout << "-> " << flag << " draws saved, total time usage: " << elapsed_time << " seconds.\n";
                        // p.increment();
                    }
                    break;
                }
            }
        }
        else
        {
            continue;
        }
    }
    Rcout << "* All draws are done\n";

    return List::create(
        _["B_saved"] = B_draw,
        _["beta_saved"] = beta_draw,
        _["Sigma_saved"] = Sigma_draw
        // _["IRF_saved"] = IRF_saved
    );
}
// [[Rcpp::export]]
List impose_NSR(List restrictions, arma::mat Y, arma::mat X, arma::cube B_draw, arma::cube beta_draw, arma::cube Sigma_draw, int plag, int hor, int M)
{
    int nvar = Y.n_cols;
    int T_est = Y.n_rows;
    int init_draw = B_draw.n_slices;
    // save draws
    arma::cube B_NSR(nvar, nvar, init_draw);
    arma::cube Sigma_NSR(nvar, nvar, init_draw);
    arma::cube beta_NSR(nvar * plag + 1, nvar, init_draw);
    arma::vec weights(init_draw);
    // flag and supervisor
    int flag = 0;
    // loop on draw
    for (int i = 0; i < init_draw; i++)
    {
        arma::mat beta_1draw = beta_draw.slice(i);
        arma::mat B_1draw = B_draw.slice(i);
        arma::mat Sigma_1draw = Sigma_draw.slice(i);
        arma::mat B_inv_1draw = inv(trans(B_1draw));
        arma::mat eps_1draw = (Y - X * beta_1draw) * B_inv_1draw;
        Psi_IRF tmp = IRF_compute_internal(beta_1draw, B_1draw, hor, nvar, plag);
        arma::cube IRF_1draw = tmp.IRF;
        bool save_me = check_NSR(restrictions, IRF_1draw, eps_1draw);
        if (save_me)
        {
            B_NSR.slice(flag) = B_1draw;
            beta_NSR.slice(flag) = beta_1draw;
            Sigma_NSR.slice(flag) = Sigma_1draw;
            weights(flag) = compute_importance_weight(restrictions, IRF_1draw, M, T_est, nvar);
            flag++;
        }
    }
    Rcout << "* " << flag << " draws in " << init_draw << " draws satisfy NSR, resampling using importance weights...\n";
    // re-weight using importance sampling
    if (flag >= 1)
    {
        arma::vec x = linspace(0, flag - 1, flag);
        NumericVector xx = wrap(x);
        NumericVector my_weights = wrap(weights.head(flag));
        NumericVector new_id = sample(xx, flag, true, my_weights);
        arma::uvec save_id = as<uvec>(new_id);
        arma::uvec unique_draw = arma::unique(save_id);
        int unique_draw_num = unique_draw.n_elem;
        arma::cube B_saved = B_NSR.slices(save_id);
        arma::cube beta_saved = beta_NSR.slices(save_id);
        arma::cube Sigma_saved = Sigma_NSR.slices(save_id);

        Rcout << "* " << unique_draw_num << " unique draws in resample" << "\n";
        return List::create(
            _["B_NSR"] = B_saved,
            _["beta_NSR"] = beta_saved,
            _["Sigma_NSR"] = Sigma_saved
        );
    }
    else
    {
        Rcout << "* " << "NSR fails to get solutions, return NA T_T \n";
        return List::create(
            _["B_NSR"] = B_draw,
            _["beta_NSR"] = beta_draw,
            _["Sigma_NSR"] = Sigma_draw);
    }
}

// // [[Rcpp::export]]
// List impose_all_restrictions(List alpha_post, List Sigma_post, List restrictions, mat Y, mat X, int draw, int plag, int hor, int M)
// {
//     int nvar = Y.n_cols;
//     int T_est = Y.n_rows;
//     // save draws
//     cube B_draw(nvar, nvar, draw);
//     cube beta_draw(nvar * plag + 1, nvar, draw);
//     cube Sigma_draw(nvar, nvar, draw);
//     vec weights(draw);
//     // unpack priors
//     vec alpha_post_mean = as<vec>(alpha_post[0]);
//     mat alpha_post_cov = as<mat>(alpha_post[1]);
//     mat Sigma_post_mean = as<mat>(Sigma_post[0]);
//     int nu_post = as<int>(Sigma_post[1]);
//     // flag and supervisor
//     int flag = 0;
//     int report = draw / 10;
//     int try_per_beta_sigma = 100;
//     auto start_time = high_resolution_clock::now();
//     // loop on draw
//     while (flag < draw)
//     {
//         mat Sigma_1draw = iwishrnd(Sigma_post_mean, nu_post);
//         mat N0 = inv(X.t() * X);
//         mat NT = kron(Sigma_1draw, N0);
//         mat cov = chol(NT, "lower");
//         vec alpha_1draw = alpha_post_mean + cov * randn(alpha_post_mean.n_elem, 1);
//         mat beta_1draw = reshape(alpha_1draw, nvar * plag + 1, nvar);
//         mat u_1draw = Y - X * beta_1draw;
//         mat P = chol(Sigma_1draw, "lower");
//         double maxroot = max_root(beta_1draw);
//         if (maxroot <= 1)
//         {
//             for (int i = 0; i < try_per_beta_sigma; i++)
//             {
//                 mat X = randn(nvar, nvar);
//                 mat Q, R;
//                 qr(Q, R, X);
//                 Q = Q * diagmat(sign(R.diag()));
//                 mat B_1draw = P * Q;
//                 mat B_inv_1draw = inv(trans(B_1draw));
//                 mat eps_1draw = u_1draw * B_inv_1draw;
//                 Psi_IRF tmp = IRF_compute_internal(beta_1draw, B_1draw, hor, nvar, plag);
//                 cube IRF_1draw = tmp.IRF;
//                 // check sign restrictions, to decide save or not
//                 bool save_me = check_all_restrictions(restrictions, IRF_1draw, eps_1draw);
//                 if (save_me)
//                 {
//                     B_draw.slice(flag) = B_1draw;
//                     beta_draw.slice(flag) = beta_1draw;
//                     Sigma_draw.slice(flag) = Sigma_1draw;
//                     weights(flag) = compute_importance_weight(restrictions, IRF_1draw, M, T_est, nvar);
//                     flag++;
//                     if (flag % report == 0)
//                     {
//                         auto current_time = high_resolution_clock::now();
//                         auto elapsed_time = duration_cast<seconds>(current_time - start_time).count();
//                         Rcout << "-> " << flag << " draws saved, total time usage: " << elapsed_time << " seconds.\n";
//                     }
//                 }
//             }
//         }
//     }
//     Rcout << "* All draws are done\n";
//     // re-weight using importance sampling
//     vec x = linspace(0, draw - 1, draw);
//     NumericVector xx = wrap(x);
//     NumericVector my_weights = wrap(weights);
//     NumericVector new_id = sample(xx, draw, true, my_weights);
//     uvec save_id = as<uvec>(new_id);
//     cube B_saved = B_draw.slices(save_id);
//     cube beta_saved = beta_draw.slices(save_id);
//     cube Sigma_saved = Sigma_draw.slices(save_id);
//
//     return List::create(
//         _["B_saved"] = B_draw,
//         _["beta_saved"] = beta_draw,
//         _["Sigma_saved"] = Sigma_draw);
// }
