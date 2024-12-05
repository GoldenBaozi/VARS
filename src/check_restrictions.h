#ifndef CHECK_RESTRICTIONS_H
#define CHECK_RESTRICTIONS_H

#include <RcppArmadillo.h>
#include "misc.h"
#include "VARtools.h"

using namespace Rcpp;
using namespace arma;

bool check_SR_and_EBR(List restrictions, cube IRF);
bool check_NSR(List restrictions, cube IRF, mat eps);
double compute_importance_weight(List restrictions, cube IRF, int M, int row, int col);
double check_all_restrictions(List restrictions, cube IRF, mat eps);

#endif
