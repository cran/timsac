#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(mulbarf) (double*,int*,int*,double*,int*,double*,double*,double*,double*,double*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

SEXP mulbar(SEXP y, SEXP n, SEXP d, SEXP calb, SEXP morder)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17;
    int *i1,*i2,*i3,*i4;

    SEXP ans =  R_NilValue, mean = R_NilValue, var = R_NilValue, v = R_NilValue, aic = R_NilValue, daic = R_NilValue;
    SEXP m = R_NilValue, aicm = R_NilValue, vm = R_NilValue, w1 = R_NilValue, w2 = R_NilValue,  a = R_NilValue, b = R_NilValue,  g = R_NilValue, h = R_NilValue,  e= R_NilValue, aicb = R_NilValue;

    double *xmean, *xvar, *xv, *xaic, *xdaic, *xaicm, *xvm, *xw1, *xw2, *xa, *xb, *xg, *xh, *xe, *xaicb = NULL;
    int    *xm= NULL;
    int    i, nd, nd2, mo, mo1;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(d);
    d2 = NUMERIC_POINTER(calb);
    i3 = INTEGER_POINTER(morder);

    nd = *i2;
    nd2 = nd * nd;
    mo = *i3;
    mo1 = mo + 1;
    PROTECT(ans = allocVector(VECSXP, 16));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 1, var = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 2, v = allocVector(REALSXP, mo1));
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, mo1));
    SET_VECTOR_ELT(ans, 4, daic = allocVector(REALSXP, mo1));
    SET_VECTOR_ELT(ans, 5, m = allocVector(INTSXP, 1)); 
    SET_VECTOR_ELT(ans, 6, aicm = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 7, vm = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 8, w1 = allocVector(REALSXP, mo1));
    SET_VECTOR_ELT(ans, 9, w2 = allocVector(REALSXP, mo)); 
    SET_VECTOR_ELT(ans, 10, a = allocVector(REALSXP, nd2*mo));
    SET_VECTOR_ELT(ans, 11, b = allocVector(REALSXP, nd2*mo)); 
    SET_VECTOR_ELT(ans, 12, g = allocVector(REALSXP, nd2*mo));
    SET_VECTOR_ELT(ans, 13, h = allocVector(REALSXP, nd2*mo));
    SET_VECTOR_ELT(ans, 14, e = allocVector(REALSXP, nd2));
    SET_VECTOR_ELT(ans, 15, aicb = allocVector(REALSXP, 1));

    d3 = NUMERIC_POINTER(mean);
    d4 = NUMERIC_POINTER(var);
    d5 = NUMERIC_POINTER(v);
    d6 = NUMERIC_POINTER(aic);
    d7 = NUMERIC_POINTER(daic);
    i4 = INTEGER_POINTER(m);
    d8 = NUMERIC_POINTER(aicm);
    d9 = NUMERIC_POINTER(vm);
    d10 = NUMERIC_POINTER(w1);
    d11 = NUMERIC_POINTER(w2);
    d12 = NUMERIC_POINTER(a);
    d13 = NUMERIC_POINTER(b);
    d14 = NUMERIC_POINTER(g);
    d15 = NUMERIC_POINTER(h);
    d16 = NUMERIC_POINTER(e);
    d17 = NUMERIC_POINTER(aicb);

    F77_CALL(mulbarf) (d1,i1,i2,d2,i3,d3,d4,d5,d6,d7,i4,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17);

    xmean = REAL(mean);
    xvar = REAL(var);
    xv = REAL(v);
    xaic = REAL(aic);
    xdaic = REAL(daic);
    xm = INTEGER(m);
    xaicm = REAL(aicm);
    xvm = REAL(vm);
    xw1 = REAL(w1);
    xw2 = REAL(w2);
    xa = REAL(a);
    xb = REAL(b);
    xg = REAL(g);
    xh = REAL(h);
    xe = REAL(e);
    xaicb = REAL(aicb);

    for(i=0; i<nd; i++) xmean[i] = d3[i];
    for(i=0; i<nd; i++) xvar[i] = d4[i];
    for(i=0; i<mo1; i++) xv[i] = d5[i];
    for(i=0; i<mo1; i++) xaic[i] = d6[i];
    for(i=0; i<mo1; i++) xdaic[i] = d7[i];
    *xm = *i4;
    *xaicm = *d8;
    *xvm = *d9;
    for(i=0; i<mo1; i++) xw1[i] = d10[i];
    for(i=0; i<mo; i++) xw2[i] = d11[i];
    for(i=0; i<nd2*mo; i++) xa[i] = d12[i];
    for(i=0; i<nd2*mo; i++) xb[i] = d13[i];
    for(i=0; i<nd2*mo; i++) xg[i] = d14[i];
    for(i=0; i<nd2*mo; i++) xh[i] = d15[i];
    for(i=0; i<nd2; i++) xe[i] = d16[i];
    *xaicb = *d17;

    UNPROTECT(1);

    return ans;
}
