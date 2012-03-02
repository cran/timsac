#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(mlocarf) (double*,int*,int*,int*,int*,int*,double*,double*,double*,int*,double*,int*,int*,double*,int*,int*,int*,double*,double*,int*,double*,double*);

SEXP mlocar(SEXP y, SEXP n, SEXP morder, SEXP span, SEXP cnst, SEXP ns)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12;

    SEXP ans =  R_NilValue, mean = R_NilValue, var = R_NilValue,  a = R_NilValue,  mf = R_NilValue, sdf = R_NilValue, ks = R_NilValue, ke = R_NilValue, pxx = R_NilValue, ld1 = R_NilValue, ld2 = R_NilValue;
    SEXP ms = R_NilValue, sdms = R_NilValue, aics = R_NilValue, mp = R_NilValue, sdmp = R_NilValue, aicp = R_NilValue;

    double *xmean, *xvar, *xa, *xsdf, *xpxx, *xsdms, *xaics, *xsdmp, *xaicp = NULL;
    int    *xmf, *xks, *xke, *xld1, *xld2, *xms, *xmp = NULL;
    int    i, mo, cnt, nns, k;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(morder);
    i3 = INTEGER_POINTER(span);
    i4 = INTEGER_POINTER(cnst);
    i5 = INTEGER_POINTER(ns);

    mo = *i2;
    cnt = *i4;
    nns = *i5;
    k = (mo+cnt)*nns;
    PROTECT(ans = allocVector(VECSXP, 16));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, var = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, a = allocVector(REALSXP, k)); 
    SET_VECTOR_ELT(ans, 3, mf = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 4, sdf = allocVector(REALSXP, nns));
    SET_VECTOR_ELT(ans, 5, ks = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 6, ke = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 7, pxx = allocVector(REALSXP, 121*nns));
    SET_VECTOR_ELT(ans, 8, ld1 = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 9, ld2 = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 10, ms = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 11, sdms = allocVector(REALSXP, nns)); 
    SET_VECTOR_ELT(ans, 12, aics = allocVector(REALSXP, nns));
    SET_VECTOR_ELT(ans, 13, mp = allocVector(INTSXP, nns));
    SET_VECTOR_ELT(ans, 14, sdmp = allocVector(REALSXP, nns)); 
    SET_VECTOR_ELT(ans, 15, aicp = allocVector(REALSXP, nns));

    d2 = NUMERIC_POINTER(mean);
    d3 = NUMERIC_POINTER(var);
    d4 = NUMERIC_POINTER(a);
    i6 = INTEGER_POINTER(mf);
    d5 = NUMERIC_POINTER(sdf);
    i7 = INTEGER_POINTER(ks);
    i8 = INTEGER_POINTER(ke);
    d6 = NUMERIC_POINTER(pxx);
    i9 = INTEGER_POINTER(ld1);
    i10 = INTEGER_POINTER(ld2);
    i11 = INTEGER_POINTER(ms);
    d7 = NUMERIC_POINTER(sdms);
    d8 = NUMERIC_POINTER(aics);
    i12 = INTEGER_POINTER(mp);
    d9 = NUMERIC_POINTER(sdmp);
    d10 = NUMERIC_POINTER(aicp);

    F77_CALL(mlocarf) (d1,i1,i2,i3,i4,i5,d2,d3,d4,i6,d5,i7,i8,d6,i9,i10,i11,d7,d8,i12,d9,d10);

    xmean = REAL(mean);
    xvar = REAL(var);
    xa = REAL(a);
    xmf = INTEGER(mf);
    xsdf = REAL(sdf);
    xks = INTEGER(ks);
    xke = INTEGER(ke);
    xpxx = REAL(pxx);
    xld1 = INTEGER(ld1);
    xld2 = INTEGER(ld2);
    xms = INTEGER(ms);
    xsdms = REAL(sdms);
    xaics = REAL(aics);
    xmp = INTEGER(mp);
    xsdmp = REAL(sdmp);
    xaicp = REAL(aicp);

    *xmean = *d2;
    *xvar = *d3;
    for(i=0; i<k; i++) xa[i] = d4[i];
    for(i=0; i<nns; i++) xmf[i] = i6[i];
    for(i=0; i<nns; i++) xsdf[i] = d5[i];
    for(i=0; i<nns; i++) xks[i] = i7[i];
    for(i=0; i<nns; i++) xke[i] = i8[i];
    for(i=0; i<121*nns; i++) xpxx[i] = d6[i];
    for(i=0; i<nns; i++) xld1[i] = i9[i];
    for(i=0; i<nns; i++) xld2[i] = i10[i];
    for(i=0; i<nns; i++) xms[i] = i11[i];
    for(i=0; i<nns; i++) xsdms[i] = d7[i];
    for(i=0; i<nns; i++) xaics[i] = d8[i];
    for(i=0; i<nns; i++) xmp[i] = i12[i];
    for(i=0; i<nns; i++) xsdmp[i] = d9[i];
    for(i=0; i<nns; i++) xaicp[i] = d10[i];

    UNPROTECT(1);

    return ans;
}
