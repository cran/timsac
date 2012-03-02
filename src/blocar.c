#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void  F77_NAME(blocarf)(double*,int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,int*,int*,double*);

SEXP blocar(SEXP y, SEXP n, SEXP sorder, SEXP span, SEXP ns)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int    *i1,*i2,*i3,*i4,*i5,*i6;
    int    i, odr, nsp, nsp2;

    SEXP ans =  R_NilValue, mean = R_NilValue, var = R_NilValue, aic = R_NilValue, bw = R_NilValue, b = R_NilValue, a = R_NilValue, v = R_NilValue, ks = R_NilValue, ke = R_NilValue, pxx = R_NilValue;
    double *xmean, *xvar, *xaic, *xbw, *xb, *xa, *xv, *xpxx = NULL;
    int   *xks, *xke = NULL;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(sorder);
    i3 = INTEGER_POINTER(span);
    i4 = INTEGER_POINTER(ns);

    odr = *i2;
    nsp = *i4;
    nsp2 = nsp * nsp;
    PROTECT(ans = allocVector(VECSXP, 10));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, var = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, aic = allocVector(REALSXP, nsp2)); 
    SET_VECTOR_ELT(ans, 3, bw = allocVector(REALSXP, nsp2));
    SET_VECTOR_ELT(ans, 4, b = allocVector(REALSXP, nsp*odr));
    SET_VECTOR_ELT(ans, 5, a = allocVector(REALSXP, nsp*odr)); 
    SET_VECTOR_ELT(ans, 6, v = allocVector(REALSXP, nsp)); 
    SET_VECTOR_ELT(ans, 7, ks = allocVector(INTSXP, nsp));
    SET_VECTOR_ELT(ans, 8, ke = allocVector(INTSXP, nsp));
    SET_VECTOR_ELT(ans, 9, pxx = allocVector(REALSXP, nsp*121)); 

    d2 = NUMERIC_POINTER(mean);
    d3 = NUMERIC_POINTER(var);
    d4 = NUMERIC_POINTER(aic);
    d5 = NUMERIC_POINTER(bw);
    d6 = NUMERIC_POINTER(b);
    d7 = NUMERIC_POINTER(a);
    d8 = NUMERIC_POINTER(v);
    i5 = INTEGER_POINTER(ks);
    i6 = INTEGER_POINTER(ke);
    d9 = NUMERIC_POINTER(pxx);

    F77_CALL(blocarf) (d1,i1,i2,i3,i4,d2,d3,d4,d5,d6,d7,d8,i5,i6,d9);

    xmean = REAL(mean);
    xvar = REAL(var);
    xaic = REAL(aic);
    xbw = REAL(bw);
    xb = REAL(b);
    xa = REAL(a);
    xv = REAL(v);
    xks = INTEGER(ks);
    xke = INTEGER(ke);
    xpxx = REAL(pxx);

    *xmean = *d2;
    *xvar = *d3;
    for(i=0; i<nsp2; i++) xaic[i] = d4[i];
    for(i=0; i<nsp2; i++) xbw[i] = d5[i];
    for(i=0; i<nsp*odr; i++) xb[i] = d6[i];
    for(i=0; i<nsp*odr; i++) xa[i] = d7[i];
    for(i=0; i<nsp; i++) xv[i] = d8[i];
    for(i=0; i<nsp; i++) xks[i] = i5[i];
    for(i=0; i<nsp; i++) xke[i] = i6[i];
    for(i=0; i<nsp*121; i++) xpxx[i] = d9[i];

    UNPROTECT(1);

    return ans;
}
