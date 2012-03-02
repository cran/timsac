#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void  F77_NAME(blomarf) (double*,int*,int*,double*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,int*,int*);

SEXP blomar(SEXP y, SEXP n, SEXP d, SEXP calb, SEXP morder, SEXP span, SEXP ns)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int    *i1,*i2,*i3,*i4,*i5,*i6,*i7;
    int    i, nd, nd2, odr, nsp, nsp2;

    SEXP ans =  R_NilValue, mean = R_NilValue, var = R_NilValue, bw = R_NilValue, raic = R_NilValue, a = R_NilValue, e = R_NilValue, aic = R_NilValue, ks = R_NilValue, ke = R_NilValue;
    double *xmean, *xvar, *xbw, *xraic, *xa, *xe, *xaic = NULL;
    int   *xks, *xke = NULL;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(d);
    d2 = NUMERIC_POINTER(calb);
    i3 = INTEGER_POINTER(morder);
    i4 = INTEGER_POINTER(span);
    i5 = INTEGER_POINTER(ns);

    nd = *i2;
    nd2 = nd*nd;
    odr = *i3;
    nsp = *i5;
    nsp2 = nsp * nsp;
    PROTECT(ans = allocVector(VECSXP, 9));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 1, var = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 2, bw = allocVector(REALSXP, nsp2));
    SET_VECTOR_ELT(ans, 3, raic = allocVector(REALSXP, nsp2)); 
    SET_VECTOR_ELT(ans, 4, a = allocVector(REALSXP, nd2*odr*nsp)); 
    SET_VECTOR_ELT(ans, 5, e = allocVector(REALSXP, nd2*nsp)); 
    SET_VECTOR_ELT(ans, 6, aic = allocVector(REALSXP, nsp)); 
    SET_VECTOR_ELT(ans, 7, ks = allocVector(INTSXP, nsp));
    SET_VECTOR_ELT(ans, 8, ke = allocVector(INTSXP, nsp));

    d3 = NUMERIC_POINTER(mean);
    d4 = NUMERIC_POINTER(var);
    d5 = NUMERIC_POINTER(bw);
    d6 = NUMERIC_POINTER(raic);
    d7 = NUMERIC_POINTER(a);
    d8 = NUMERIC_POINTER(e);
    d9 = NUMERIC_POINTER(aic);
    i6 = INTEGER_POINTER(ks);
    i7 = INTEGER_POINTER(ke);

    F77_CALL(blomarf) (d1,i1,i2,d2,i3,i4,i5,d3,d4,d5,d6,d7,d8,d9,i6,i7);

    xmean = REAL(mean);
    xvar = REAL(var);
    xbw = REAL(bw);
    xraic = REAL(raic);
    xa = REAL(a);
    xe = REAL(e);
    xaic = REAL(aic);
    xks = INTEGER(ks);
    xke = INTEGER(ke);

    for(i=0; i<nd; i++) xmean[i] = d3[i];
    for(i=0; i<nd; i++) xvar[i] = d4[i];
    for(i=0; i<nsp2; i++) xbw[i] = d5[i];
    for(i=0; i<nsp2; i++) xraic[i] = d6[i];
    for(i=0; i<nd2*odr*nsp; i++) xa[i] = d7[i];
    for(i=0; i<nd2*nsp; i++) xe[i] = d8[i];
    for(i=0; i<nsp; i++) xaic[i] = d9[i];
    for(i=0; i<nsp; i++) xks[i] = i6[i];
    for(i=0; i<nsp; i++) xke[i] = i7[i];

    UNPROTECT(1);

    return ans;
}

