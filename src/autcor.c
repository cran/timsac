#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(autcorf)(double*, int*, double*, double*, int*, double*);

SEXP autcor(SEXP y, SEXP n, SEXP lag1)
{
    double *d1, *d2, *d3, *d4;
    int *i1, *i2;

    SEXP ans =  R_NilValue, acov = R_NilValue, acor = R_NilValue, mean = R_NilValue;
    double *xacov, *xacor, *xmean = NULL;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(lag1);

    int len = *i2;
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, acov = allocVector(REALSXP, len));
    SET_VECTOR_ELT(ans, 1, acor = allocVector(REALSXP, len));
    SET_VECTOR_ELT(ans, 2, mean = allocVector(REALSXP, 1)); 

    d2 = NUMERIC_POINTER(acov);
    d3 = NUMERIC_POINTER(acor);
    d4 = NUMERIC_POINTER(mean);

    F77_CALL(autcorf)(d1,i1,d2,d3,i2,d4);

    xacov = REAL(acov);
    xacor = REAL(acor);
    xmean = REAL(mean);

    for(int i=0; i<len; i++) xacov[i] = d2[i];
    for(int i=0; i<len; i++) xacor[i] = d3[i];
    *xmean = *d4;

    UNPROTECT(1);

    return ans;
}
