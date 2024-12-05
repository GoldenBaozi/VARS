// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "misc.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

double max_root(mat beta)
{
    int nvar = beta.n_cols;
    int plag = (beta.n_rows - 1) / nvar;
    mat beta_1 = beta.t();
    mat beta_lag = beta_1.submat(0, 1, beta_1.n_rows - 1, beta_1.n_cols - 1);
    mat large_eye = eye(nvar * (plag - 1), nvar * (plag - 1));
    mat large_zero = zeros((plag - 1) * nvar, nvar);
    mat beta_compact = join_cols(beta_lag, join_rows(large_eye, large_zero));
    cx_vec roots = eig_gen(beta_compact);
    double max_root = max(abs(roots));
    return max_root;
}

mat flatten_cube(cube obj)
{
    int hor = obj.n_slices;
    int nvar = obj.n_rows;
    mat out(nvar * nvar, hor);
    for (int i = 0; i < nvar * nvar; i++)
    {
        int row = i / nvar;
        int col = i % nvar;
        cube tmp = reshape(obj.subcube(row, col, 0, row, col, hor - 1), 1, hor, 1);
        out.row(i) = tmp.slice(0);
    }
    return out;
}

cube mat_to_cube(mat obj)
{
    int hor = obj.n_cols;
    int nvar = sqrt(obj.n_rows);
    cube out(nvar, nvar, hor);

    for (int i = 0; i < nvar * nvar - 1; i++)
    {
        int row = i / nvar;
        int col = i % nvar;
        cube tmp(1, hor, 1);
        tmp.slice(0) = obj.row(i);
        tmp = reshape(tmp, 1, 1, hor);
        out.tube(row, col) = tmp;
    }

    return out;
}

vec wild_boot(mat eps)
{
    int T_est = eps.n_rows;
    // int nvar = eps.n_cols;
    NumericVector boot_tmp = 2 * rbinom(T_est, 1, 0.5);
    vec boot_idx = as<vec>(boot_tmp);
    // mat boot_out = repmat(boot_core, 1, nvar) % eps;
    return boot_idx;
}

vec block_boot(mat eps, int size)
{
    int T_est = eps.n_rows;
    // int nvar = eps.n_cols;
    int block_num = T_est / size;
    int flag = 0;
    mat block_tmp(0, 0);
    if (T_est % size != 0)
    {
        block_num++;
        flag++; // block_num * size > T_est, be careful with the last block
    }
    vec pool = linspace(0, block_num - 1, block_num);
    vec boot_idx = sort(randi<vec>(block_num, distr_param(0, block_num - 1)));
    // for (int i = 0; i < draw.n_elem; i++)
    // {
    //     int start = draw(i) * size;
    //     int end = start + size - 1;
    //     if (draw(i) == block_num - 1)
    //     {
    //         end = T_est - 1;
    //     }

    //     block_tmp = join_cols(block_tmp, eps.rows(start, end));
    // }
    // mat block_out = block_tmp.rows(0, T_est - 1);
    return boot_idx;
}

mat shuffle_c(mat obj, vec idx, string method, int block_size)
{
    int T_est = obj.n_rows;
    int nvar = obj.n_cols;
    mat boot_out(T_est, nvar);
    int block_num = T_est / block_size;
    if (T_est % block_size != 0)
    {
        block_num++;
    }
    if (method == "wild")
    {
        boot_out = repmat(idx, 1, nvar) % obj;
    }
    else
    {
        mat block_tmp(0, 0);
        for (int i = 0; i < idx.n_elem; i++)
        {
            int start = idx(i) * block_size;
            int end = start + block_size - 1;
            if (idx(i) == block_num - 1)
            {
                end = T_est - 1;
            }

            block_tmp = join_cols(block_tmp, obj.rows(start, end));
        }
        boot_out = block_tmp.rows(0, T_est-1);
    }
    return boot_out;
}

mat construct_data(mat Y_start, mat beta, mat u, int T)
{
    int plag = Y_start.n_rows;
    mat Y = Y_start;
    // bool cons = (beta.n_rows % plag == 1);

    for (int i = 0; i < T - plag; i++)
    {
        mat X_1 = trans((Y.rows(i, plag - 1 + i).t()));
        mat Y_1 = X_1 * beta + u.row(i);
        Y = join_cols(Y, Y_1);
    }
    return Y;
}
