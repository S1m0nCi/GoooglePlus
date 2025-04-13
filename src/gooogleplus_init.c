#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>

extern SEXP groupCoordinateDescent(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"groupCoordinateDescent", (DL_FUNC) &groupCoordinateDescent, 10},
    {NULL, NULL, 0}
};

void R_init_gooogleplus(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

// Euclidean norm
double norm(double *x, int p) {
    double x_norm = 0;
    for (int j=0; j<p; j++) x_norm = x_norm + pow(x[j],2);
    x_norm = sqrt(x_norm);
    return(x_norm);
}

// Soft-thresholding operator
double S(double z, double l) {
    if (z > l) return(z-l);
    if (z < -l) return(z+l);
    return(0);
}

// Firm-thresholding operator for MCP penalty
double F(double z, double l, double gamma) {
  if (fabs(z) <= l * gamma) {
      return(S(z, l) / (1 - 1 / gamma));
  }
  return(z);
}

// SCAD-modified firm-thresholding operator
double Fs(double z, double l, double gamma) {
    if (fabs(z) <= 2 * l) {
        return(S(z, l));
    } if (fabs(z) <= gamma * l) {
        return(S(z, l * gamma / (gamma - 1)) / (1 - 1 / (gamma - 1)));
    } 
    return(z);
}

// Gaussian loss
double gLoss(double *r, int n) {
    double l = 0;
    for (int i=0;i<n;i++) l = l + pow(r[i],2);
    return(l);
}

// Dot product for matrix multiplication
double dotproduct(double *u, double *v, int len) {
    double w = 0;
    for (int i = 0; i < len; i++) {
        w += u[i]*v[i];
    }
    return(w);
}