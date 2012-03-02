#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(thirmof) (int*,int*,double*,double*,double*,double*,double*);

SEXP thirmo(SEXP n, SEXP lag, SEXP y)
{
    double *d1,*d2,*d3,*d4,*d5;
    int    *i1,*i2;

    SEXP ans =  R_NilValue,  mean = R_NilValue, acov = R_NilValue, acor = R_NilValue, mnt = R_NilValue;
    double *xmean, *xcov, *xcor, *xmnt = NULL;
    int   i, lg1;

    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(lag);
    d1 = NUMERIC_POINTER(y);

    lg1 = *i2 + 1;
    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, acov = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 2, acor = allocVector(REALSXP, lg1)); 
    SET_VECTOR_ELT(ans, 3, mnt = allocVector(REALSXP, lg1*lg1)); 

    d2 = NUMERIC_POINTER(mean);
    d3 = NUMERIC_POINTER(acov);
    d4 = NUMERIC_POINTER(acor);
    d5 = NUMERIC_POINTER(mnt);

    F77_CALL(thirmof) (i1,i2,d1,d2,d3,d4,d5);

    xmean = REAL(mean);
    xcov = REAL(acov);
    xcor = REAL(acor);
    xmnt = REAL(mnt);

    *xmean = *d2;
    for(i=0; i<lg1; i++) xcov[i] = d3[i];
    for(i=0; i<lg1; i++) xcor[i] = d4[i];
    for(i=0; i<lg1*lg1; i++) xmnt[i] = d5[i];

    UNPROTECT(1);

    return ans;
}
