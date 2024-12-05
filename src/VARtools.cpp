// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <chrono>
#include "misc.h"
#include "VARtools.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// struct Psi_IRF
// {
//     cube Psi;
//     cube IRF;
// };

// [[Rcpp::export]]
List IRF_compute(arma::mat beta, arma::mat B, int hor, int nvar, int plag)
{
    bool cons = (beta.n_rows % nvar == 1);
    arma::mat beta_1 = beta.t();
    arma::mat beta_lag;
    if (cons)
    {
        beta_lag = beta_1.submat(0, 1, beta_1.n_rows - 1, beta_1.n_cols - 1);
    }
    else
    {
        beta_lag = beta_1;
    }
    arma::mat large_eye = arma::eye(nvar * (plag - 1), nvar * (plag - 1));
    arma::mat large_zero = arma::zeros((plag - 1) * nvar, nvar);
    arma::mat beta_compact = arma::join_cols(beta_lag, arma::join_rows(large_eye, large_zero));
    arma::mat irf_trans = arma::join_rows(arma::eye<arma::mat>(nvar, nvar), arma::zeros<arma::mat>(nvar, nvar * (plag - 1)));
    arma::cube Psi(nvar, nvar, hor + 1);
    arma::cube IRF(nvar, nvar, hor + 1);
    for (int h = 0; h < hor + 1; h++)
    {
        Psi.slice(h) = irf_trans * arma::powmat(beta_compact, h) * irf_trans.t();
        IRF.slice(h) = Psi.slice(h) * B;
    }
    List out = List::create(
        _["Psi"] = Psi,
        _["IRF"] = IRF);
    return out;
}

Psi_IRF IRF_compute_internal(mat beta, mat B, int hor, int nvar, int plag)
{
    bool cons = (beta.n_rows % nvar == 1);
    mat beta_1 = beta.t();
    mat beta_lag;
    if (cons)
    {
        beta_lag = beta_1.submat(0, 1, beta_1.n_rows - 1, beta_1.n_cols - 1);
    }
    else
    {
        beta_lag = beta_1;
    }
    mat large_eye = eye(nvar * (plag - 1), nvar * (plag - 1));
    mat large_zero = zeros((plag - 1) * nvar, nvar);
    mat beta_compact = join_cols(beta_lag, join_rows(large_eye, large_zero));
    mat irf_trans = join_rows(eye<mat>(nvar, nvar), zeros<mat>(nvar, nvar * (plag - 1)));
    cube Psi(nvar, nvar, hor + 1);
    cube IRF(nvar, nvar, hor + 1);
    for (int h = 0; h < hor + 1; h++)
    {
        Psi.slice(h) = irf_trans * powmat(beta_compact, h) * irf_trans.t();
        IRF.slice(h) = Psi.slice(h) * B;
    }

    Psi_IRF out;
    out.Psi = Psi;
    out.IRF = IRF;

    return out;
}

// [[Rcpp::export]]
arma::cube FEVD_compute(arma::mat Sigma, arma::mat B, arma::cube Psi, int nvar, int hor)
{
    int nstep = hor + 1;
    arma::cube MSE(nvar, nvar, nstep);
    arma::cube MSE_shock(nvar, nvar, nstep);
    arma::cube FEVD(nvar, nvar, nstep);
    for (int i = 0; i < nvar; i++)
    {
        MSE.slice(0) = Sigma;
        MSE_shock.slice(0) = B.col(i) * trans(B.col(i));
        for (int j = 1; j < nstep; j++)
        {
            MSE.slice(j) = MSE.slice(j - 1) + Psi.slice(j) * Sigma * trans(Psi.slice(j));
            MSE_shock.slice(j) = MSE_shock.slice(j - 1) + Psi.slice(j) * MSE_shock.slice(0) * trans(Psi.slice(j));
            FEVD.slice(j).col(i) = arma::diagvec(MSE_shock.slice(j)) / arma::diagvec(MSE.slice(j));
        }
    }
    return FEVD;
}

// [[Rcpp::export]]
double HDC(int start, int end, int shock, int res_var, arma::cube IRF, arma::mat eps)
{
    int t = start - 1;
    int h = end - 1;
    int i = res_var - 1;
    int j = shock - 1;
    double HDC = 0.0;
    for (int k = t; k < h + 1; k++)
    {
        HDC += eps(k, j) * IRF(i, j, h - k);
    }
    return HDC;
}

// [[Rcpp::export]]
arma::vec HDC_ts(int start, int end, int shock, int res_var, arma::cube IRF, arma::mat eps)
{
    vec HDC_ts(end - start + 1);
    for (int t = start; t < end + 1; t++)
    {
        HDC_ts(t - start) = HDC(start, t, shock, res_var, IRF, eps);
    }
    return HDC_ts;
}

//' @details compute HD of all shocks to which_var from start to end.
vec get_HDs(int start, int end, int which_var, cube IRF, mat eps)
{
    int hmax = end - start + 1;
    int end_id = end - 1;
    int nvar = eps.n_cols;
    int var_id = which_var - 1;
    mat HD(nvar, hmax);

    for (int h = 0; h < hmax - 1; h++)
    {
        for (int j = 0; j < nvar - 1; j++)
        {
            HD(j, h) = IRF(var_id, j, h) * eps(end_id - h, j);
        }
    }

    vec HD_out = sum(HD, 1);
    return HD_out;
}

//' @details compute HD of all shocks to which_var from start to start, start + 1, ..., end.
//' @return a matrix, of which row j is HD of shock j to which_var during the given period (a time series)
// [[Rcpp::export]]
arma::mat get_HDs_ts(int start, int end, int which_var, arma::cube IRF, arma::mat eps)
{
    int nvar = eps.n_cols;
    int hmax = end - start + 1;
    arma::mat HDs(nvar, hmax);

    for (int h = 0; h < hmax - 1; h++)
    {
        HDs.col(h) = get_HDs(start, start + h, which_var, IRF, eps);
    }
    return HDs;
}

// [[Rcpp::export]]
List bootstrap_c(arma::mat Y_start, arma::mat beta, arma::mat u, arma::mat IV, std::string identify, std::string method, int save, int T)
{
    int T_est = u.n_rows;
    int nvar = u.n_cols;
    int n_para = beta.n_rows;
    int plag = Y_start.n_rows;
    arma::cube beta_save(n_para, nvar, save);
    arma::cube u_save(T_est, nvar, save);
    arma::cube Sigma_save(nvar, nvar, save);
    arma::cube B_save(nvar, nvar, save);
    arma::cube eps_save(T_est, nvar, save);
    bool cons = (n_para % nvar == 1);
    for (int i = 0; i < save; i++)
    {
        arma::mat u_tmp(T_est, nvar);
        arma::mat IV_tmp(T_est, nvar);
        arma::vec idx;
        if (method == "wild")
        {
            idx = wild_boot(u);
        }
        else
        {
            int size = plag;
            idx = block_boot(u, size);
        }
        u_tmp = shuffle_c(u, idx, method, plag);
        arma::mat data = construct_data(Y_start, beta, u_tmp, T);
        arma::mat Y = data.rows(plag, T - 1);
        arma::mat X(0, 0);
        for (int j = 1; j <= plag; j++)
        {
            X = arma::join_rows(X, data.rows(plag - j, T - 1 - j));
        }
        if (cons)
        {
            X = arma::join_rows(arma::ones(T_est), X);
        }
        arma::mat beta_tmp = arma::inv(X.t() * X) * X.t() * Y;
        // mat u_tmp = Y - X * beta_tmp;
        arma::mat Sigma_tmp = u_tmp.t() * u_tmp;
        beta_save.slice(i) = beta_tmp;
        u_save.slice(i) = u_tmp;
        Sigma_save.slice(i) = Sigma_tmp;
        if (identify == "IV")
        {
            IV_tmp = shuffle_c(IV, idx, method, plag);
            arma::mat B_tmp = u_tmp.t() * IV_tmp / T_est;
            arma::vec B_diag = arma::diagvec(B_tmp).t();
            for (int i = 0; i < nvar; i++)
            {
                if (B_diag(i) == 0)
                {
                    B_diag(i) = 1;
                }
            }
            B_save.slice(i) = B_tmp / repmat(B_diag, nvar, 1);
        }
        else
        {
            B_save.slice(i) = chol(Sigma_tmp, "lower");
        }
        eps_save.slice(i) = u_tmp * inv(B_save.slice(i).t());
    }

    return List::create(
        _["beta.set"] = beta_save,
        _["U.set"] = u_save,
        _["Sigma.set"] = Sigma_save,
        _["B.set"] = B_save,
        _["eps.set"] = eps_save
    );
}

// [[Rcpp::export]]
List IRF_compute_batch(arma::cube beta_draw, arma::cube B_draw, int hor)
{
    int nvar = B_draw.n_cols;
    int draw_num = B_draw.n_slices;
    int plag = (beta_draw.n_rows - 1) / nvar;
    arma::cube Psi_draws(nvar * nvar, hor + 1, draw_num);
    arma::cube IRF_draws(nvar * nvar, hor + 1, draw_num);
    for (int i = 0; i < draw_num; i++)
    {
        arma::mat beta = beta_draw.slice(i);
        arma::mat B = B_draw.slice(i);
        Psi_IRF IRF_tmp = IRF_compute_internal(beta, B, hor, nvar, plag);
        arma::cube Psi_1draw = IRF_tmp.Psi;
        arma::cube IRF_1draw = IRF_tmp.IRF;
        Psi_draws.slice(i) = flatten_cube(Psi_1draw);
        IRF_draws.slice(i) = flatten_cube(IRF_1draw);
    }
    return List::create(
        _["Psi"] = Psi_draws,
        _["IRF"] = IRF_draws
    );
}

// [[Rcpp::export]]
arma::cube FEVD_compute_batch(arma::cube beta_draw, arma::cube B_draw, arma::cube Sigma_draw, arma::cube Psi_draw, int hor)
{
    int nvar = B_draw.n_cols;
    int draw_num = B_draw.n_slices;
    // int plag = (beta_draw.n_rows - 1) / nvar;
    arma::cube FEVD_draws(nvar * nvar, hor + 1, draw_num);

    for (int i = 0; i < draw_num; i++)
    {
        arma::mat beta = beta_draw.slice(i);
        arma::mat B = B_draw.slice(i);
        arma::mat Sigma = Sigma_draw.slice(i);
        arma::cube Psi = mat_to_cube(Psi_draw.slice(i));
        arma::cube FEVD_1draw = FEVD_compute(Sigma, B, Psi, nvar, hor);
        FEVD_draws.slice(i) = flatten_cube(FEVD_1draw);
    }
    return FEVD_draws;
}

// [[Rcpp::export]]
arma::cube HDs_compute_batch(int start, int end, int which_var, arma::cube IRF_draw, arma::cube eps_draw)
{
    int nvar = eps_draw.n_cols;
    int hmax = end - start + 1;
    int draw = eps_draw.n_slices;
    arma::cube HDs_draw(nvar, hmax, draw);

    for (int i = 0; i < draw; i++)
    {
        arma::mat IRF_tmp = IRF_draw.slice(i);
        arma::mat eps = eps_draw.slice(i);
        arma::cube IRF = mat_to_cube(IRF_tmp);
        HDs_draw.slice(i) = get_HDs_ts(start, end, which_var, IRF, eps);
    }

    return HDs_draw;
}
