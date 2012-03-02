#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(exsarf)(double*,int*,int*,double*,double*,double*,double*,double*,int*,double*,double*,double*,double*,double*,int*);

SEXP exsar( SEXP y, SEXP n, SEXP morder )
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11;
    int    *i1,*i2,*i3,*i4;

    SEXP ans =  R_NilValue, mean = R_NilValue, var = R_NilValue, v = R_NilValue, aic = R_NilValue, daic = R_NilValue, m =  R_NilValue, aicm = R_NilValue, sdm1 = R_NilValue, a1 = R_NilValue, sdm2 = R_NilValue, a2 = R_NilValue, ier = R_NilValue;

    double *xmean, *xvar, *xv, *xaic, *xdaic, *xaicm, *xsdm1, *xa1, *xsdm2, *xa2 = NULL;
    int      *xm, *xier = NULL;
    int      i, mm, mm1;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(morder);

    mm = *i2;
    mm1 = mm + 1;
    PROTECT(ans = allocVector(VECSXP, 12));
    SET_VECTOR_ELT(ans, 0,mean = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 1, var = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, v = allocVector(REALSXP, mm1));
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, mm1)); 
    SET_VECTOR_ELT(ans, 4, daic = allocVector(REALSXP, mm1));
    SET_VECTOR_ELT(ans, 5, m = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 6, aicm = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 7, sdm1 = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 8, a1 = allocVector(REALSXP, mm));
    SET_VECTOR_ELT(ans, 9, sdm2 = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 10, a2 = allocVector(REALSXP, mm));
    SET_VECTOR_ELT(ans, 11, ier = allocVector(INTSXP, 1));

    d2 = NUMERIC_POINTER(mean);
    d3 = NUMERIC_POINTER(var);
    d4 = NUMERIC_POINTER(v);
    d5 = NUMERIC_POINTER(aic);
    d6 = NUMERIC_POINTER(daic);
    i3 = INTEGER_POINTER(m);
    d7 = NUMERIC_POINTER(aicm);
    d8 = NUMERIC_POINTER(sdm1);
    d9 = NUMERIC_POINTER(a1);
    d10 = NUMERIC_POINTER(sdm2);
    d11 = NUMERIC_POINTER(a2);
    i4 = INTEGER_POINTER(ier);

    F77_CALL(exsarf) (d1,i1,i2,d2,d3,d4,d5,d6,i3,d7,d8,d9,d10,d11,i4);

    xmean = REAL(mean);
    xvar = REAL(var);
    xv = REAL(v);
    xaic = REAL(aic);
    xdaic = REAL(daic);
    xm = INTEGER(m);
    xaicm = REAL(aicm);
    xsdm1 = REAL(sdm1);
    xa1 = REAL(a1);
    xsdm2 = REAL(sdm2);
    xa2 = REAL(a2);
    xier = INTEGER(ier);

    *xmean = *d2;
    *xvar = *d3;
    for(i=0; i<mm1; i++) xv[i] = d4[i];
    for(i=0; i<mm1; i++) xaic[i] = d5[i];
    for(i=0; i<mm1; i++) xdaic[i] = d6[i];
    *xm = *i3;
    *xaicm = *d7;
    *xsdm1 = *d8;
    for(i=0; i<mm; i++) xa1[i] = d9[i];
    *xsdm2 = *d10;
    for(i=0; i<mm; i++) xa2[i] = d11[i];
    *xier = *i4;

    UNPROTECT(1);
    return ans; 
}


