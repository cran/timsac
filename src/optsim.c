#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(optsimf) (int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

SEXP OptsimC(SEXP ns, SEXP order, SEXP ir, SEXP il, SEXP trans, SEXP gamma, SEXP gain, SEXP noise)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14;
    int *i1,*i2,*i3,*i4;

    SEXP ans = R_NilValue, x =  R_NilValue, y =  R_NilValue, xmean =  R_NilValue, ymean =  R_NilValue, x2sum =  R_NilValue, y2sum =  R_NilValue, x2mean =  R_NilValue, y2mean =  R_NilValue, xvar =  R_NilValue, yvar =  R_NilValue;
    double *xx, *xy, *xxm, *xym, *xx2s, *xy2s, *xx2m, *xy2m, *xxv, *xyv = NULL;
    int   i, nns, nr, nl;

    i1 = INTEGER_POINTER(ns);
    i2 = INTEGER_POINTER(order);
    i3 = INTEGER_POINTER(ir);
    i4 = INTEGER_POINTER(il);
    d1 = NUMERIC_POINTER(trans);
    d2 = NUMERIC_POINTER(gamma);
    d3 = NUMERIC_POINTER(gain);
    d4 = NUMERIC_POINTER(noise);

    nns = *i1;
    nr = *i3;
    nl = *i4;
    PROTECT(ans = allocVector(VECSXP, 10));
    SET_VECTOR_ELT(ans, 0, x = allocVector(REALSXP, nr*nns));
    SET_VECTOR_ELT(ans, 1, y = allocVector(REALSXP, nl*nns));
    SET_VECTOR_ELT(ans, 2, xmean = allocVector(REALSXP, nr));
    SET_VECTOR_ELT(ans, 3, ymean = allocVector(REALSXP, nl));
    SET_VECTOR_ELT(ans, 4, x2sum = allocVector(REALSXP, nr));
    SET_VECTOR_ELT(ans, 5, y2sum = allocVector(REALSXP, nl));
    SET_VECTOR_ELT(ans, 6, x2mean = allocVector(REALSXP, nr));
    SET_VECTOR_ELT(ans, 7, y2mean = allocVector(REALSXP, nl));
    SET_VECTOR_ELT(ans, 8, xvar = allocVector(REALSXP, nr));
    SET_VECTOR_ELT(ans, 9, yvar = allocVector(REALSXP, nl));

    d5 = NUMERIC_POINTER(x);
    d6 = NUMERIC_POINTER(y);
    d7 = NUMERIC_POINTER(xmean);
    d8 = NUMERIC_POINTER(ymean);
    d9 = NUMERIC_POINTER(x2sum);
    d10 = NUMERIC_POINTER(y2sum);
    d11 = NUMERIC_POINTER(x2mean);
    d12 = NUMERIC_POINTER(y2mean);
    d13 = NUMERIC_POINTER(xvar);
    d14 = NUMERIC_POINTER(yvar);

    F77_CALL(optsimf) (i1,i2,i3,i4,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14);

    xx = REAL(x);
    xy = REAL(y);
    xxm = REAL(xmean);
    xym = REAL(ymean);
    xx2s = REAL(x2sum);
    xy2s = REAL(y2sum);
    xx2m = REAL(x2mean);
    xy2m = REAL(y2mean);
    xxv = REAL(xvar);
    xyv = REAL(yvar);

    for(i=0; i<nr*nns; i++) xx[i] = d5[i];
    for(i=0; i<nl*nns; i++) xy[i] = d6[i];
    for(i=0; i<nr; i++) xxm[i] = d7[i];
    for(i=0; i<nl; i++) xym[i] = d8[i];
    for(i=0; i<nr; i++) xx2s[i] = d9[i];
    for(i=0; i<nl; i++) xy2s[i] = d10[i];
    for(i=0; i<nr; i++) xx2m[i] = d11[i];
    for(i=0; i<nl; i++) xy2m[i] = d12[i];
    for(i=0; i<nr; i++) xxv[i] = d13[i];
    for(i=0; i<nl; i++) xyv[i] = d14[i];

    UNPROTECT(1);

    return ans;
}



