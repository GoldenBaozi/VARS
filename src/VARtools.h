#ifndef VARTOOLS_H
#define VARTOOLS_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

struct Psi_IRF
{
    cube Psi;
    cube IRF;
};

List IRF_compute(mat beta, mat B, int hor, int nvar, int plag);
Psi_IRF IRF_compute_internal(mat beta, mat B, int hor, int nvar, int plag);
cube FEVD_compute(mat Sigma, mat B, cube Psi, int nvar, int hor);
double HDC(int start, int end, int shock, int res_var, cube IRF, mat eps);
vec HDC_ts(int start, int end, int shock, int res_var, cube IRF, mat eps);
vec get_HDs(int start, int end, int which_var, cube IRF, mat eps);
mat get_HDs_ts(int start, int end, int which_var, cube IRF, mat eps);
List bootstrap_c(mat Y_start, mat beta, mat u, mat IV, string identify, string method, int save, int T);
List IRF_compute_batch(cube beta_draw, cube B_draw, int hor);
cube FEVD_compute_batch(cube beta_draw, cube B_draw, cube Sigma_draw, cube Psi_draw, int hor);
cube HDs_compute_batch(int start, int end, int which_var, cube IRF_draw, cube eps_draw);

#endif
