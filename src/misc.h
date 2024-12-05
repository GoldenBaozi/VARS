#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h>

using namespace arma;
using namespace std;

double max_root(mat beta);
mat flatten_cube(cube obj);
cube mat_to_cube(mat obj);
vec wild_boot(mat eps);
vec block_boot(mat eps, int size);
mat shuffle_c(mat obj, vec idx, string method, int block_size);
mat construct_data(mat Y_start, mat beta, mat u, int T);

#endif
