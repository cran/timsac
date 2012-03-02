#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(mlomarf) (double*,int*,int*,double*,int*,int*,int*,int*,double*,double*,int*,int*,int*,double*,int*,double*,int*,double*,double*,double*,int*,int*);

SEXP mlomar(SEXP y, SEXP n, SEXP d, SEXP calb, SEXP morder, SEXP span, SEXP cnst, SEXP ns)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13;

    SEXP ans =  R_NilValue, mean = R_NilValue, var = R_NilValue, ld1 = R_NilValue, ld2 = R_NilValue;
    SEXP ms = R_NilValue, aicm = R_NilValue, mp = R_NilValue, aicc = R_NilValue, mf = R_NilValue, aic = R_NilValue,  a = R_NilValue,  e = R_NilValue,  ks = R_NilValue, ke = R_NilValue;

    double *xmean, *xvar, *xaicm, *xaicc, *xaic, *xa, *xe = NULL;
    int    *xld1, *xld2, *xms, *xmp , *xmf, *xks, *xke= NULL;
    int    i, nd, mo, nns, k;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(d);
    d2 = NUMERIC_POINTER(calb);
    i3 = INTEGER_POINTER(morder);
    i4 = INTEGER_POINTER(span);
    i5 = INTEGER_POINTER(cnst);
    i6 = INTEGER_POINTER(ns);

    nd = *i2;
    mo = *i3;
    nns = *i6;
    k = nd*nd*nns;
    PROTECT(ans = allocVector(VECSXP, 14));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 1, var = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 2, ld1 = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 3, ld2 = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 4, ms = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 5, aicm = allocVector(REALSXP, nns)); 
    SET_VECTOR_ELT(ans, 6, mp = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 7, aicc = allocVector(REALSXP, nns));
    SET_VECTOR_ELT(ans, 8, mf = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 9, aic = allocVector(REALSXP, nns)); 
    SET_VECTOR_ELT(ans, 10, a = allocVector(REALSXP, k*mo));
    SET_VECTOR_ELT(ans, 11, e = allocVector(REALSXP, k)); 
    SET_VECTOR_ELT(ans, 12, ks = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 13, ke = allocVector(INTSXP, nns));

    d3 = NUMERIC_POINTER(mean);
    d4 = NUMERIC_POINTER(var);
    i7 = INTEGER_POINTER(ld1);
    i8 = INTEGER_POINTER(ld2);
    i9 = INTEGER_POINTER(ms);
    d5 = NUMERIC_POINTER(aicm);
    i10 = INTEGER_POINTER(mp);
    d6 = NUMERIC_POINTER(aicc);
    i11 = INTEGER_POINTER(mf);
    d7 = NUMERIC_POINTER(aic);
    d8 = NUMERIC_POINTER(a);
    d9 = NUMERIC_POINTER(e);
    i12 = INTEGER_POINTER(ks);
    i13 = INTEGER_POINTER(ke);

    F77_CALL(mlomarf) (d1,i1,i2,d2,i3,i4,i5,i6,d3,d4,i7,i8,i9,d5,i10,d6,i11,d7,d8,d9,i12,i13);

    xmean = REAL(mean);
    xvar = REAL(var);
    xld1 = INTEGER(ld1);
    xld2 = INTEGER(ld2);
    xms = INTEGER(ms);
    xaicm = REAL(aicm);
    xmp = INTEGER(mp);
    xaicc = REAL(aicc);
    xmf = INTEGER(mf);
    xaic = REAL(aic);
    xa = REAL(a);
    xe = REAL(e);
    xks = INTEGER(ks);
    xke = INTEGER(ke);

    *xmean = *d3;
    *xvar = *d4;
    for(i=0; i<nd; i++) xld1[i] = i7[i];
    for(i=0; i<nd; i++) xld2[i] = i8[i];
    for(i=0; i<nns; i++) xms[i] = i9[i];
    for(i=0; i<nns; i++) xaicm[i] = d5[i];
    for(i=0; i<nns; i++) xmp[i] = i10[i];
    for(i=0; i<nns; i++) xaicc[i] = d6[i];
    for(i=0; i<nns; i++) xmf[i] = i11[i];
    for(i=0; i<nns; i++) xaic[i] = d7[i];
    for(i=0; i<k*mo; i++) xa[i] = d8[i];
    for(i=0; i<k; i++) xe[i] = d9[i];
    for(i=0; i<nns; i++) xks[i] = i12[i];
    for(i=0; i<nns; i++) xke[i] = i13[i];

    UNPROTECT(1);

    return ans;
}
