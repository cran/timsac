#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(mulspef) (int*, int*, int*, int*, double*, double*, double*, double*, double*, double*);

SEXP mulspe(SEXP n, SEXP d, SEXP lag1, SEXP lag3, SEXP cov)
{
    double *d1,*d2,*d3,*d4,*d5,*d6;
    int *i1,*i2,*i3,*i4;

    SEXP ans =  R_NilValue,  spec1 = R_NilValue, spec2 = R_NilValue, stat = R_NilValue,  coh1 = R_NilValue, coh2 = R_NilValue;
    double *xspec1, *xspec2, *xstat, *xcoh1, *xcoh2 = NULL;
    int   i, nd, nd2, lg1;

    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(d);
    i3 = INTEGER_POINTER(lag1);
    i4 = INTEGER_POINTER(lag3);
    d1 = NUMERIC_POINTER(cov);

    nd = *i2;
    nd2 = nd * nd;
    lg1 = *i3;
    PROTECT(ans = allocVector(VECSXP, 5));
    SET_VECTOR_ELT(ans, 0, spec1 = allocVector(REALSXP, lg1*nd2));
    SET_VECTOR_ELT(ans, 1, spec2 = allocVector(REALSXP, lg1*nd2));
    SET_VECTOR_ELT(ans, 2, stat = allocVector(REALSXP, lg1*nd)); 
    SET_VECTOR_ELT(ans, 3, coh1 = allocVector(REALSXP, lg1*nd2));
    SET_VECTOR_ELT(ans, 4, coh2 = allocVector(REALSXP, lg1*nd2)); 

    d2 = NUMERIC_POINTER(spec1);
    d3 = NUMERIC_POINTER(spec2);
    d4 = NUMERIC_POINTER(stat);
    d5 = NUMERIC_POINTER(coh1);
    d6 = NUMERIC_POINTER(coh2);

    F77_CALL(mulspef) (i1,i2,i3,i4,d1,d2,d3,d4,d5,d6);

    xspec1 = REAL(spec1);
    xspec2 = REAL(spec2);
    xstat = REAL(stat);
    xcoh1 = REAL(coh1);
    xcoh2 = REAL(coh2);

    for(i=0; i<lg1*nd2; i++) xspec1[i] = d2[i];
    for(i=0; i<lg1*nd2; i++) xspec2[i] = d3[i];
    for(i=0; i<lg1*nd; i++) xstat[i] = d4[i];
    for(i=0; i<lg1*nd2; i++) xcoh1[i] = d5[i];
    for(i=0; i<lg1*nd2; i++) xcoh2[i] = d6[i];

    UNPROTECT(1);

    return ans;
}

