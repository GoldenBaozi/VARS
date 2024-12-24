library(Rcpp)
library(RcppArmadillo)

sourceCpp(code="
  // [[Rcpp::depends(RcppArmadillo)]]
  #include <RcppArmadillo.h>
  using namespace arma;
  using namespace Rcpp;

  // [[Rcpp::export()]]
  vec get_HDs(int start, int end, int which_var, arma::cube IRF, arma::mat eps)
  {
      int hmax = end - start + 1;
      int end_id = end - 1;
      int nvar = eps.n_cols;
      int var_id = which_var - 1;
      arma::mat HD(nvar, hmax);

      for (int h = 0; h < hmax; h++)
      {
          for (int j = 0; j < nvar; j++)
          {
              HD(j, h) = IRF(var_id, j, h) * eps(end_id - h, j);
          }
      }

      arma::vec HD_out = sum(HD, 1);
      return HD_out;
  }
  ")

IRFs <- array(1, dim = c(3,3,2))
eps <- matrix(c(1,2,3), 1, 3)
start <- end <- 1
which.var <- 3

get_HDs(start, end, which.var, IRFs, eps)
