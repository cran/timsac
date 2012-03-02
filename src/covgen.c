#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(covgenf) (int*,int*,double*,double*,double*,double*);

SEXP covgen(SEXP lag, SEXP k, SEXP f, SEXP gain)
{
    double *d1,*d2,*d3,*d4;
    int *i1,*i2;
    int i, lag1;

    SEXP ans =  R_NilValue, acov = R_NilValue, acor = R_NilValue;
    double *xacov, *xacor = NULL;

    i1 = INTEGER_POINTER(lag);
    i2 = INTEGER_POINTER(k);
    d1 = NUMERIC_POINTER(f);
    d2 = NUMERIC_POINTER(gain);

    lag1 = *i1+1;
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, acov = allocVector(REALSXP, lag1));
    SET_VECTOR_ELT(ans, 1, acor = allocVector(REALSXP, lag1));

    d3 = NUMERIC_POINTER(acov);
    d4 = NUMERIC_POINTER(acor);

    F77_CALL(covgenf) (i1,i2,d1,d2,d3,d4);

    xacov = REAL(acov);
    xacor = REAL(acor);

    for(i=0; i<lag1; i++) xacov[i] = d3[i];
    for(i=0; i<lag1; i++) xacor[i] = d4[i];

    UNPROTECT(1);

    return ans;
}
