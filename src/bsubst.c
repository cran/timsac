#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void  F77_NAME(bsubstf) (double*,int*,int*,int*,int*,int*,int*,int*,double*,double*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

SEXP bsubst(SEXP y, SEXP n, SEXP mtype, SEXP lag, SEXP nreg, SEXP cstep, SEXP reg, SEXP tlag)
{
    double *d1,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17,*d18,*d19,*d20;
    double *d21,*d22,*d23,*d24,*d25,*d26,*d27,*d28,*d29;
    int    *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9;
    int    i, nn, k, k1, nc;

    SEXP ans =  R_NilValue, ymean = R_NilValue, yvar = R_NilValue, m = R_NilValue, aicm = R_NilValue;
    SEXP vm =  R_NilValue, a1 = R_NilValue, v = R_NilValue, aic = R_NilValue, daic = R_NilValue;
    SEXP aicb = R_NilValue, vb =  R_NilValue, ek = R_NilValue, a2 = R_NilValue, ind = R_NilValue;
    SEXP c = R_NilValue, c1 = R_NilValue, c2 = R_NilValue, b = R_NilValue, eicmin = R_NilValue;
    SEXP esum = R_NilValue, npm = R_NilValue, npmnreg = R_NilValue, e = R_NilValue, mean = R_NilValue;
    SEXP var = R_NilValue, skew = R_NilValue, peak = R_NilValue, cov = R_NilValue, pxx = R_NilValue;
    double *xymean, *xyvar, *xaicm, *xvm, *xa1, *xv, *xaic, *xdaic, *xaicb, *xvb, *xek, *xa2 = NULL;
    double *xc, *xc1, *xc2, *xb, *xeicmin, *xesum,  *xnpm, *xnpmnreg, *xe, *xmean, *xvar = NULL;
    double *xskew, *xpeak, *xcov, *xpxx = NULL;
    int   *xm, *xind = NULL;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(mtype);
    i3 = INTEGER_POINTER(lag);
    i4 = INTEGER_POINTER(nreg);
    i5 = INTEGER_POINTER(cstep);
    i6 = INTEGER_POINTER(reg);
    i7 = INTEGER_POINTER(tlag);

    nn = *i1;
    k = *i4;
    k1 = k+1;
    nc = *i5;
    PROTECT(ans = allocVector(VECSXP, 29));
    SET_VECTOR_ELT(ans, 0, ymean = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, yvar = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, m = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 3, aicm = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 4, vm = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 5, a1 = allocVector(REALSXP, k)); 
    SET_VECTOR_ELT(ans, 6, v = allocVector(REALSXP, k1)); 
    SET_VECTOR_ELT(ans, 7, aic = allocVector(REALSXP, k1));
    SET_VECTOR_ELT(ans, 8, daic = allocVector(REALSXP, k1));
    SET_VECTOR_ELT(ans, 9, aicb = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 10, vb = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 11, ek = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 12, a2 = allocVector(REALSXP, k));
    SET_VECTOR_ELT(ans, 13, ind = allocVector(INTSXP, k)); 
    SET_VECTOR_ELT(ans, 14, c = allocVector(REALSXP, k)); 
    SET_VECTOR_ELT(ans, 15, c1 = allocVector(REALSXP, k1)); 
    SET_VECTOR_ELT(ans, 16, c2 = allocVector(REALSXP, k)); 
    SET_VECTOR_ELT(ans, 17, b = allocVector(REALSXP, k));
    SET_VECTOR_ELT(ans, 18, eicmin = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 19, esum = allocVector(REALSXP, k1));
    SET_VECTOR_ELT(ans, 20, npm = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 21, npmnreg = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 22, e = allocVector(REALSXP, nn*nc));
    SET_VECTOR_ELT(ans, 23, mean = allocVector(REALSXP, nc)); 
    SET_VECTOR_ELT(ans, 24, var = allocVector(REALSXP, nc)); 
    SET_VECTOR_ELT(ans, 25, skew = allocVector(REALSXP, nc)); 
    SET_VECTOR_ELT(ans, 26, peak = allocVector(REALSXP, nc)); 
    SET_VECTOR_ELT(ans, 27, cov = allocVector(REALSXP, 101));
    SET_VECTOR_ELT(ans, 28, pxx = allocVector(REALSXP, 121));

    d3 = NUMERIC_POINTER(ymean);
    d4 = NUMERIC_POINTER(yvar);
    i8 = INTEGER_POINTER(m);
    d5 = NUMERIC_POINTER(aicm);
    d6 = NUMERIC_POINTER(vm);
    d7 = NUMERIC_POINTER(a1);
    d8 = NUMERIC_POINTER(v);
    d9 = NUMERIC_POINTER(aic);
    d10 = NUMERIC_POINTER(daic);
    d11 = NUMERIC_POINTER(aicb);
    d12 = NUMERIC_POINTER(vb);
    d13 = NUMERIC_POINTER(ek);
    d14 = NUMERIC_POINTER(a2);
    i9 = INTEGER_POINTER(ind);
    d15 = NUMERIC_POINTER(c);
    d16 = NUMERIC_POINTER(c1);
    d17 = NUMERIC_POINTER(c2);
    d18 = NUMERIC_POINTER(b);
    d19 = NUMERIC_POINTER(eicmin);
    d20 = NUMERIC_POINTER(esum);
    d21 = NUMERIC_POINTER(npm);
    d22 = NUMERIC_POINTER(npmnreg);
    d23 = NUMERIC_POINTER(e);
    d24 = NUMERIC_POINTER(mean);
    d25 = NUMERIC_POINTER(var);
    d26 = NUMERIC_POINTER(skew);
    d27 = NUMERIC_POINTER(peak);
    d28 = NUMERIC_POINTER(cov);
    d29 = NUMERIC_POINTER(pxx);

    F77_CALL(bsubstf) (d1,i1,i2,i3,i4,i5,i6,i7,d3,d4,i8,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,i9,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26,d27,d28,d29);

    xymean = REAL(ymean);
    xyvar = REAL(yvar);
    xm = INTEGER(m);
    xaicm = REAL(aicm);
    xvm = REAL(vm);
    xa1 = REAL(a1);
    xv = REAL(v);
    xaic = REAL(aic);
    xdaic = REAL(daic);
    xaicb = REAL(aicb);
    xvb = REAL(vb);
    xek = REAL(ek);
    xa2 = REAL(a2);
    xind = INTEGER(ind);
    xc = REAL(c);
    xc1 = REAL(c1);
    xc2 = REAL(c2);
    xb = REAL(b);
    xeicmin = REAL(eicmin);
    xesum = REAL(esum);
    xnpm = REAL(npm);
    xnpmnreg = REAL(npmnreg);
    xe = REAL(e);
    xmean = REAL(mean);
    xvar = REAL(var);
    xskew = REAL(skew);
    xpeak = REAL(peak);
    xcov = REAL(cov);
    xpxx = REAL(pxx);

    *xymean = *d3;
    *xyvar = *d4;
    *xm = *i8;
    *xaicm = *d5;
    *xvm = *d6;
    for(i=0; i<k; i++) xa1[i] = d7[i];
    for(i=0; i<k1; i++) xv[i] = d8[i];
    for(i=0; i<k1; i++) xaic[i] = d9[i];
    for(i=0; i<k1; i++) xdaic[i] = d10[i];
    *xaicb = *d11;
    *xvb = *d12;
    *xek = *d13;
    for(i=0; i<k; i++) xa2[i] = d14[i];
    for(i=0; i<k; i++) xind[i] = i9[i];
    for(i=0; i<k; i++) xc[i] = d15[i];
    for(i=0; i<k1; i++) xc1[i] = d16[i];
    for(i=0; i<k; i++) xc2[i] = d17[i];
    for(i=0; i<k; i++) xb[i] = d18[i];
    *xeicmin = *d19;
    for(i=0; i<k1; i++) xesum[i] = d20[i];
    *xnpm = *d21;
    *xnpmnreg = *d22;
    for(i=0; i<nn*nc; i++) xe[i] = d23[i];
    for(i=0; i<nc; i++) xmean[i] = d24[i];
    for(i=0; i<nc; i++) xvar[i] = d25[i];
    for(i=0; i<nc; i++) xskew[i] = d26[i];
    for(i=0; i<nc; i++) xpeak[i] = d27[i];
    for(i=0; i<101; i++) xcov[i] = d28[i];
    for(i=0; i<121; i++) xpxx[i] = d29[i];

    UNPROTECT(1);

    return ans;
}


