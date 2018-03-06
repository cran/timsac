#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(auspecf)(int*, int*, double*, double*, double*, double*);

SEXP AuspecC(SEXP n, SEXP lag1, SEXP acov)
{
    double *d1, *d2, *d3, *d4;
    int *i1, *i2;
    int i;

    SEXP ans =  R_NilValue, spec1 = R_NilValue, spec2 = R_NilValue, stat = R_NilValue;
    double *xspec1, *xspec2, *xstat = NULL;

    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(lag1);
    d1 = NUMERIC_POINTER(acov);

    int lag = *i2;
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, spec1 = allocVector(REALSXP, lag));
    SET_VECTOR_ELT(ans, 1, spec2 = allocVector(REALSXP, lag));
    SET_VECTOR_ELT(ans, 2, stat = allocVector(REALSXP, lag));

    d2 = NUMERIC_POINTER(spec1);
    d3 = NUMERIC_POINTER(spec2);
    d4 = NUMERIC_POINTER(stat);

    F77_CALL(auspecf)(i1,i2,d1,d2,d3,d4);

    xspec1 = REAL(spec1);
    xspec2 = REAL(spec2);
    xstat = REAL(stat);

    for(i=0; i<lag; i++) xspec1[i] = d2[i];
    for(i=0; i<lag; i++) xspec2[i] = d3[i];
    for(i=0; i<lag; i++) xstat[i] = d4[i];

    UNPROTECT(1);

    return ans;
}

