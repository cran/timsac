#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(unibarf)  (double*,int*,int*,double*,double*,double*,double*,double*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

SEXP unibar(SEXP y, SEXP n, SEXP arorder)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17;
    int    *i1,*i2,*i3;

    SEXP ans =  R_NilValue, mean = R_NilValue, var = R_NilValue, v = R_NilValue, aic = R_NilValue, daic = R_NilValue, m = R_NilValue, aicm = R_NilValue, vm = R_NilValue, pa = R_NilValue, bw = R_NilValue,  sbw = R_NilValue, pab = R_NilValue,  aicb = R_NilValue, vb = R_NilValue,  pn = R_NilValue, a = R_NilValue, pxx = R_NilValue;

    double *xmean, *xvar, *xv, *xaic, *xdaic, *xaicm, *xvm, *xpa, *xbw, *xsbw, *xpab, *xaicb, *xvb, *xpn, *xa , *xpxx= NULL;
    int    *xm= NULL;
    int    i, na, na1;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(arorder);

    na = *i2;
    na1 = na + 1;
    PROTECT(ans = allocVector(VECSXP, 17));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, var = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, v = allocVector(REALSXP, na1));
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, na1));
    SET_VECTOR_ELT(ans, 4, daic = allocVector(REALSXP, na1));
    SET_VECTOR_ELT(ans, 5, m = allocVector(INTSXP, 1)); 
    SET_VECTOR_ELT(ans, 6, aicm = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 7, vm = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 8, pa = allocVector(REALSXP, na));
    SET_VECTOR_ELT(ans, 9, bw = allocVector(REALSXP, na1)); 
    SET_VECTOR_ELT(ans, 10, sbw = allocVector(REALSXP, na));
    SET_VECTOR_ELT(ans, 11, pab = allocVector(REALSXP, na)); 
    SET_VECTOR_ELT(ans, 12, aicb = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 13, vb = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 14, pn = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 15, a = allocVector(REALSXP, na));
    SET_VECTOR_ELT(ans, 16, pxx = allocVector(REALSXP, 121));


    d2 = NUMERIC_POINTER(mean);
    d3 = NUMERIC_POINTER(var);
    d4 = NUMERIC_POINTER(v);
    d5 = NUMERIC_POINTER(aic);
    d6 = NUMERIC_POINTER(daic);
    i3 = INTEGER_POINTER(m);
    d7 = NUMERIC_POINTER(aicm);
    d8 = NUMERIC_POINTER(vm);
    d9 = NUMERIC_POINTER(pa);
    d10 = NUMERIC_POINTER(bw);
    d11 = NUMERIC_POINTER(sbw);
    d12 = NUMERIC_POINTER(pab);
    d13 = NUMERIC_POINTER(aicb);
    d14 = NUMERIC_POINTER(vb);
    d15 = NUMERIC_POINTER(pn);
    d16 = NUMERIC_POINTER(a);
    d17 = NUMERIC_POINTER(pxx);

    F77_CALL(unibarf) (d1,i1,i2,d2,d3,d4,d5,d6,i3,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17);

    xmean = REAL(mean);
    xvar = REAL(var);
    xv = REAL(v);
    xaic = REAL(aic);
    xdaic = REAL(daic);
    xm = INTEGER(m);
    xaicm = REAL(aicm);
    xvm = REAL(vm);
    xpa = REAL(pa);
    xbw = REAL(bw);
    xsbw = REAL(sbw);
    xpab = REAL(pab);
    xaicb = REAL(aicb);
    xvb = REAL(vb);
    xpn = REAL(pn);
    xa = REAL(a);
    xpxx = REAL(pxx);

    *xmean = *d2;
    *xvar = *d3;
    for(i=0; i<na1; i++) xv[i] = d4[i];
    for(i=0; i<na1; i++) xaic[i] = d5[i];
    for(i=0; i<na1; i++) xdaic[i] = d6[i];
    *xm = *i3;
    *xaicm = *d7;
    *xvm = *d8;
    for(i=0; i<na; i++) xpa[i] = d9[i];
    for(i=0; i<na1; i++) xbw[i] = d10[i];
    for(i=0; i<na; i++) xsbw[i] = d11[i];
    for(i=0; i<na; i++) xpab[i] = d12[i];
    *xaicb = *d13;
    *xvb = *d14;
    *xpn = *d15;
    for(i=0; i<na; i++) xa[i] = d16[i];
    for(i=0; i<121; i++) xpxx[i] = d17[i];

    UNPROTECT(1);

    return ans;
}

