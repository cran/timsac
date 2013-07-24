#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(mulcorf) (double*, int*, int*, int*, double*, double*, double*);

SEXP mulcor(SEXP y, SEXP n, SEXP d, SEXP lag1)
{
    double *d1, *d2, *d3, *d4;
    int    *i1,*i2,*i3;

    SEXP ans =  R_NilValue,  mean = R_NilValue, cov = R_NilValue, cor = R_NilValue;
    double *xmean, *xcov, *xcor = NULL;
    int   i, nd, ldd;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(d);
    i3 = INTEGER_POINTER(lag1);

    nd = *i2;
    ldd = (*i3) * nd * nd;
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 1, cov = allocVector(REALSXP, ldd));
    SET_VECTOR_ELT(ans, 2, cor = allocVector(REALSXP, ldd)); 

    d2 = NUMERIC_POINTER(mean);
    d3 = NUMERIC_POINTER(cov);
    d4 = NUMERIC_POINTER(cor);

    F77_CALL(mulcorf) (d1,i1,i2,i3,d2,d3,d4);

    xmean = REAL(mean);
    xcov = REAL(cov);
    xcor = REAL(cor);

    for(i=0; i<nd; i++) xmean[i] = d2[i];
    for(i=0; i<ldd; i++) xcov[i] = d3[i];
    for(i=0; i<ldd; i++) xcor[i] = d4[i];

    UNPROTECT(1);

    return ans;
}

