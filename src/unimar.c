#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(unimarf)  (double*,int*,int*,double*,double*,double*,double*,double*,int*,double*,double*,double*);

SEXP UnimarC( SEXP y, SEXP n, SEXP morder )
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int    *i1,*i2,*i3;

    SEXP ans =  R_NilValue, mean = R_NilValue, var = R_NilValue, v = R_NilValue, aic = R_NilValue, daic = R_NilValue, m =  R_NilValue, aicm = R_NilValue, vm = R_NilValue, a = R_NilValue;

    double *xmean, *xvar, *xv, *xaic, *xdaic, *xaicm, *xvm, *xa = NULL;
    int      *xm = NULL;
    int      i, mo, mo1;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(morder);

    mo = *i2;
    mo1 = mo + 1;
    PROTECT(ans = allocVector(VECSXP, 9));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, var = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, v = allocVector(REALSXP, mo1)); 
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, mo1));
    SET_VECTOR_ELT(ans, 4, daic = allocVector(REALSXP, mo1));
    SET_VECTOR_ELT(ans, 5, m = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 6, aicm = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 7, vm = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 8, a = allocVector(REALSXP, mo)); 

    d2 = NUMERIC_POINTER(mean);
    d3 = NUMERIC_POINTER(var);
    d4 = NUMERIC_POINTER(v);
    d5 = NUMERIC_POINTER(aic);
    d6 = NUMERIC_POINTER(daic);
    i3 = INTEGER_POINTER(m);
    d7 = NUMERIC_POINTER(aicm);
    d8 = NUMERIC_POINTER(vm);
    d9 = NUMERIC_POINTER(a);

    F77_CALL(unimarf) (d1,i1,i2,d2,d3,d4,d5,d6,i3,d7,d8,d9);

    xmean = REAL(mean);
    xvar = REAL(var);
    xv = REAL(v);
    xaic = REAL(aic);
    xdaic = REAL(daic);
    xm = INTEGER(m);
    xaicm = REAL(aicm);
    xvm = REAL(vm);
    xa = REAL(a);

    *xmean = *d2;
    *xvar = *d3;
    for(i=0; i<mo1; i++) xv[i] = d4[i];
    for(i=0; i<mo1; i++) xaic[i] = d5[i];
    for(i=0; i<mo1; i++) xdaic[i] = d6[i];
    *xm = *i3;
    *xaicm = *d7;
    *xvm = *d8;
    for(i=0; i<mo; i++) xa[i] = d9[i];

    UNPROTECT(1);
    return ans; 
}
