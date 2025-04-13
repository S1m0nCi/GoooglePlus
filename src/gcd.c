#include <math.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
// #include "gplus_init.c"

double norm(double *x, int p);
double S(double, double);
double F(double, double, double);
double Fs(double, double, double);
double gLoss(double *r, int n);
double dotproduct(double *u, double *v, int);

SEXP groupCoordinateDescent(
    SEXP groupMatOrth_, 
    SEXP lambda_,
    SEXP penaltyFactors_,
    SEXP initialB_,
    SEXP epsilon_,
    SEXP groupStartIndices_,
    SEXP pseudoY_,
    SEXP gamma_,
    SEXP penalty_,
    SEXP maxIterations_) {

    // get lengths
    // number of observations
    int n = length(pseudoY_);
    // number of lambda values to try
    int L = length(lambda_);
    // number of groups
    int J = length(groupStartIndices_) - 1;
    // number of coefficients needed (variables including intercepts)
    int p = length(groupMatOrth_)/n;

    // create pointers for translation from R to C
    double *X = REAL(groupMatOrth_);
    double *lambda = REAL(lambda_);
    double *penalty_factors = REAL(penaltyFactors_);
    double *start = REAL(initialB_);
    double eps = REAL(epsilon_)[0];
    int *markers = INTEGER(groupStartIndices_);
    double *Y = REAL(pseudoY_);
    double gamma = REAL(gamma_)[0];
    const char *penalty = CHAR(STRING_ELT(penalty_, 0));
    int max_iterations = INTEGER(maxIterations_)[0];
    
    // initialise the return result of all coefficient solutions
    SEXP coefficients;
    PROTECT(coefficients = allocVector(REALSXP, L*p));

    for (int j=0; j<(L*p); j++) REAL(coefficients)[j] = 0;

    // group coordinate descent
    for (int j = 0; j <= J; j++) {
        markers[j]--;
    }
    
    // get initial residuals
    double res[n];
    for (int j = 0; j < n; j++) {
        // get the jth row of X
        double row[p];
        for (int k = 0; k < p; k++) {
            row[k] = X[j + k * n];
        }
        res[j] = Y[j] - dotproduct(row, start, p);
    }

    double rss = gLoss(res, n);
    double sdy = sqrt(rss/n);

    // loop over the lambdas
    for(int i = 0; i < L; i++) {
        // initialise penalty factors
        double pf[J];
        for (int j = 0; j < J; j++) {
            pf[j] = lambda[i] * penalty_factors[j];
        }
        
        // initialise coefficients for this lambda
        double b[p];
        for (int j = 0; j < p; j++) {
            b[j] = start[j];
        }
        
        // initialise updating residuals for this lambda
        double r[n];
        for (int j = 0; j < n; j++) {
            r[j] = res[j];
        }
        
        // gcd steps
        for (int j = 0; j < max_iterations; j++) {
            double max_change = 0;
            double l0 = pf[0];
            for (int k = markers[0]; k < markers[1]; k++) {
                // get the intercept coefficient to be updated
                double b1 = b[k];
                // get its corresponding column (kth) in design matrix
                double col[n];
                for (int q = 0; q < n; q++) {
                    col[q] = X[k * n + q];
                }
                // find z
                double z = dotproduct(col, r, n)/((double)n) + b1;
                // apply thresholding operator to z
                double b2  = 0;
                if ((strcmp(penalty, "grLasso")==0) | (strcmp(penalty, "grALasso")==0)) b2 = S(z, l0);
                if (strcmp(penalty, "grMCP")==0) b2 = F(z, l0, gamma);
                if (strcmp(penalty, "grSCAD")==0) b2 = Fs(z, l0, gamma);
                // update r
                double bdiff = b2 - b1;
                for (int q = 0; q < n; q++) {
                    r[q] -= col[q] * bdiff;
                }
                if (fabs(bdiff) > max_change) {
                    max_change = fabs(bdiff);
                }
                // update coefficients for this lambda
                b[k] = b2;
            }
            for (int k = 1; k < J; k++) {
                // get the penalty factor for the kth group
                double l1 = pf[k];
                int group_size = markers[k+1] - markers[k];
                double b1[group_size];
                double z[group_size];
                // get the group coefficients and calculate z
                for (int q = markers[k]; q < markers[k+1]; q++) {
                    b1[q - markers[k]] = b[q];
                    // get the qth column
                    double col[n];
                    for (int t = 0; t < n; t++) {
                        col[t] = X[q * n + t];
                    }
                    z[q - markers[k]] = dotproduct(col, r, n)/((double)n) + b[q];
                }
                // get the norm of z
                double z_norm = norm(z, group_size);
                // find the result of its soft-thresholding operator
                double thresh = 0;
                if ((strcmp(penalty, "grLasso")==0) | (strcmp(penalty, "grALasso")==0)) thresh = S(z_norm, l1);
                if (strcmp(penalty, "grMCP")==0) thresh = F(z_norm, l1, gamma);
                if (strcmp(penalty, "grSCAD")==0) thresh = Fs(z_norm, l1, gamma);
                // find the new group coefficients
                double b2[group_size];
                double bdiff[group_size];
                for (int q = 0; q < group_size; q++) {
                    // what if z_norm is zero?
                    if (z_norm > 0) {
                        b2[q] = thresh*z[q]/z_norm;
                        bdiff[q] = b2[q] - b1[q];
                    } else {
                        b2[q] = 0;
                        bdiff[q] = 0;
                    }
                }

                // update r
                for (int q = 0; q < n; q++) {
                    double rowsect[group_size];
                    for (int t = markers[k]; t < markers[k+1]; t++) {
                        rowsect[t - markers[k]] = X[n*t + q];
                    }
                    r[q] -= dotproduct(rowsect, bdiff, group_size);
                }
                // update max_change
                for (int q = 0; q < group_size; q++) {
                    if (fabs(bdiff[q]) > max_change) {
                        max_change = fabs(bdiff[q]);
                    }
                }
                // update coefficients
                for (int q = markers[k]; q < markers[k+1]; q++) {
                    b[q] = b2[q - markers[k]];
                }
            }
            if (max_change < eps*sdy) {
                break;
            }
        }
        // append the coefficients b for each lambda to the array of coefficients
        for (int j = (i * p); j < ((i + 1) * p); j++) {
            REAL(coefficients)[j] = b[j - (i * p)];
        }
    }
    UNPROTECT(1);
    return(coefficients);
}