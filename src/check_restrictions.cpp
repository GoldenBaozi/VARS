#include <RcppArmadillo.h>
#include "misc.h"
#include "VARtools.h"
#include "check_restrictions.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

bool check_SR_and_EBR(List restrictions, cube IRF)
{
    int n = restrictions.length();
    int flag = 0;
    int num = 0;
    for (int i = 0; i < n; i++)
    {
        List res_1 = restrictions[i];
        string type = as<string>(res_1[0]);
        if (type == "SR")
        {
            num++;
            int shock = as<int>(res_1[1]) - 1;
            int var = as<int>(res_1[2]) - 1;
            int hor = as<int>(res_1[3]);
            double sign = as<double>(res_1[4]);
            mat IRF_test = IRF.slice(hor);
            double SR_test = IRF_test(var, shock) * sign;
            if (SR_test > 0)
            {
                flag++;
            }
        }
        if (type == "EBR")
        {
            num++;
            int shock = as<int>(res_1[1]) - 1;
            int var_1 = as<int>(res_1[2]) - 1;
            int var_2 = as<int>(res_1[3]) - 1;
            int hor = as<int>(res_1[4]);
            double mb = as<double>(res_1[5]);
            double lb = as<double>(res_1[6]);
            mat IRF_test = IRF.slice(hor);
            double EBR_test = IRF_test(var_1, shock) / IRF_test(var_2, shock);
            // Rcout << mb << " " << lb << " " << EBR_test << "\n";
            if (EBR_test <= mb && EBR_test >= lb)
            {
                flag++;
            }
        }
    }
    // Rcout << flag << "\n";
    bool out = (flag == num);
    return out;
}

bool check_NSR(List restrictions, cube IRF, mat eps)
{
    int n = restrictions.length();
    int flag = 0;
    int NSR_num = 0;
    for (int i = 0; i < n; i++)
    {
        List res_1 = restrictions[i];
        string type = as<string>(res_1[0]);
        if (type == "NSR")
        {
            NSR_num++;
            string NSR_type = as<string>(res_1[1]);
            if (NSR_type == "sign")
            {
                int shock = as<int>(res_1[2]) - 1;
                uvec period = as<uvec>(res_1[3]) - 1;
                double sign = as<double>(res_1[4]);
                mat eps_test = eps.submat(period(0), shock, period(period.n_elem - 1), shock) * sign;
                bool test = all(vectorise(eps_test) > 0);
                if (test)
                {
                    flag++;
                }
            }
            else
            {
                int shock = as<int>(res_1[2]);
                int var = as<int>(res_1[3]);
                int nvar = eps.n_cols;
                uvec period = as<uvec>(res_1[4]);
                int start = period(0);
                int end = period(period.n_elem - 1);
                string inten = as<string>(res_1[6]);
                vec HDs = get_HDs(start, end, var, IRF, eps);
                double target = HDs(shock - 1);
                vec HDs_test1(nvar);
                HDs_test1.fill(target);
                double HDs_test2 = sum(HDs) - target;
                bool cond_most = all(HDs_test1 - HDs >= 0);
                bool cond_least = all(HDs_test1 - HDs <= 0);
                bool cond_overwhelm = (target >= HDs_test2);
                bool cond = ((inten == "most" && cond_most) || (inten == "least" && cond_least) || (inten == "overwhelm" && cond_overwhelm) || (inten == "negligible" && !cond_overwhelm));
                if (cond)
                {
                    flag++;
                }
            }
        }
    }
    bool out = (flag == NSR_num);
    return out;
}

double check_all_restrictions(List restrictions, cube IRF, mat eps)
{
    int n = restrictions.length();
    int flag = 0;
    for (int i = 0; i < n; i++)
    {
        List res_1 = restrictions[i];
        string type = as<string>(res_1[0]);
        if (type == "SR")
        {
            int shock = as<int>(res_1[1]) - 1;
            int var = as<int>(res_1[2]) - 1;
            int hor = as<int>(res_1[3]);
            double sign = as<double>(res_1[4]);
            mat IRF_test = IRF.slice(hor);
            double SR_test = IRF_test(var, shock) * sign;
            if (SR_test > 0)
            {
                flag++;
            }
        }
        if (type == "EBR")
        {
            int shock = as<int>(res_1[1]) - 1;
            int var_1 = as<int>(res_1[2]) - 1;
            int var_2 = as<int>(res_1[3]) - 1;
            int hor = as<int>(res_1[4]);
            double mb = as<double>(res_1[5]);
            double lb = as<double>(res_1[6]);
            mat IRF_test = IRF.slice(hor);
            double EBR_test = IRF_test(var_1, shock) / IRF_test(var_2, shock);
            // Rcout << mb << " " << lb << " " << EBR_test << "\n";
            if (EBR_test <= mb && EBR_test >= lb)
            {
                flag++;
            }

            if (type == "NSR")
            {
                string NSR_type = as<string>(res_1[1]);
                if (NSR_type == "sign")
                {
                    int shock = as<int>(res_1[2]) - 1;
                    uvec period = as<uvec>(res_1[3]) - 1;
                    double sign = as<double>(res_1[4]);
                    mat eps_test = eps.submat(period(0), shock, period(period.n_elem - 1), shock) * sign;
                    bool test = all(vectorise(eps_test) > 0);
                    if (test)
                    {
                        flag++;
                    }
                }
                else
                {
                    int shock = as<int>(res_1[2]);
                    int var = as<int>(res_1[3]);
                    int nvar = eps.n_cols;
                    uvec period = as<uvec>(res_1[4]);
                    int start = period(0);
                    int end = period(period.n_elem - 1);
                    string inten = as<string>(res_1[6]);
                    vec HDs = get_HDs(start, end, var, IRF, eps);
                    double target = HDs(shock - 1);
                    vec HDs_test1(nvar);
                    HDs_test1.fill(target);
                    double HDs_test2 = sum(HDs) - target;
                    bool cond_most = all(HDs_test1 - HDs >= 0);
                    bool cond_least = all(HDs_test1 - HDs <= 0);
                    bool cond_overwhelm = (target >= HDs_test2);
                    bool cond = ((inten == "most" && cond_most) || (inten == "least" && cond_least) || (inten == "overwhelm" && cond_overwhelm) || (inten == "negligible" && !cond_overwhelm));
                    if (cond)
                    {
                        flag++;
                    }
                }
            }
        }
    }
    bool out = (flag == n);
    return out;
}

double compute_importance_weight(List restrictions, cube IRF, int M, int row, int col)
{
    int satisfy = 1; // to avoid zero weight
    for (int i = 0; i < M; i++)
    {
        mat eps_v = randn(row, col);
        bool cond = check_NSR(restrictions, IRF, eps_v);
        if (cond)
        {
            satisfy++;
        }
    }
    double weight = 1.0 / (satisfy * 1.0 / M);
    return weight;
}

