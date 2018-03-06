#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(mulmarf)  (double*,int*,int*,double*,int*,double*,double*,double*,double*,double*,int*,double*,double*,int*,int*,double*,double*,double*,double*,double*,double*,double*,int*,double*);

SEXP MulmarC( SEXP y, SEXP n, SEXP d, SEXP calb, SEXP lag )
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17;
    int    *i1,*i2,*i3,*i4,*i5,*i6,*i7;

    SEXP ans =  R_NilValue, mean = R_NilValue, var = R_NilValue, v = R_NilValue, aic = R_NilValue, daic = R_NilValue, m =  R_NilValue, aicm = R_NilValue, vm = R_NilValue, npr = R_NilValue, jnd = R_NilValue, a = R_NilValue, rv = R_NilValue, aicf = R_NilValue, ei = R_NilValue, bi = R_NilValue, matv = R_NilValue, arcoef = R_NilValue, morder = R_NilValue, aics = R_NilValue;

    double *xmean, *xvar, *xv, *xaic, *xdaic, *xaicm, *xvm, *xa, *xrv, *xaicf, *xei, *xbi, *xmatv, *xarcoef,  *xaics = NULL;
    int      *xm, *xnpr, *xjnd, *xmorder = NULL;
    int      i, nd, nd2, lg, lg1;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(d);
    d2 = NUMERIC_POINTER(calb);
    i3 = INTEGER_POINTER(lag);

    nd = *i2;
    nd2 = nd * nd;
    lg = *i3;
    lg1 = lg + 1;
    PROTECT(ans = allocVector(VECSXP, 19));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 1, var = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 2, v = allocVector(REALSXP, lg1*nd)); 
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, lg1*nd));
    SET_VECTOR_ELT(ans, 4, daic = allocVector(REALSXP, lg1*nd));
    SET_VECTOR_ELT(ans, 5, m = allocVector(INTSXP, nd));
    SET_VECTOR_ELT(ans, 6, aicm = allocVector(REALSXP, nd)); 
    SET_VECTOR_ELT(ans, 7, vm = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 8, npr = allocVector(INTSXP, nd));
    SET_VECTOR_ELT(ans, 9, jnd = allocVector(INTSXP, lg1*nd2));
    SET_VECTOR_ELT(ans, 10, a = allocVector(REALSXP, lg1*nd2)); 
    SET_VECTOR_ELT(ans, 11, rv = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 12, aicf = allocVector(REALSXP, nd)); 
    SET_VECTOR_ELT(ans, 13, ei = allocVector(REALSXP, nd2)); 
    SET_VECTOR_ELT(ans, 14, bi = allocVector(REALSXP, nd2*lg));
    SET_VECTOR_ELT(ans, 15, matv = allocVector(REALSXP, nd2)); 
    SET_VECTOR_ELT(ans, 16, arcoef = allocVector(REALSXP, nd2*lg)); 
    SET_VECTOR_ELT(ans, 17, morder = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 18, aics = allocVector(REALSXP, 1));

    d3 = NUMERIC_POINTER(mean);
    d4 = NUMERIC_POINTER(var);
    d5 = NUMERIC_POINTER(v);
    d6 = NUMERIC_POINTER(aic);
    d7 = NUMERIC_POINTER(daic);
    i4 = INTEGER_POINTER(m);
    d8 = NUMERIC_POINTER(aicm);
    d9 = NUMERIC_POINTER(vm);
    i5 = INTEGER_POINTER(npr);
    i6 = INTEGER_POINTER(jnd);
    d10 = NUMERIC_POINTER(a);
    d11= NUMERIC_POINTER(rv);
    d12= NUMERIC_POINTER(aicf);
    d13 = NUMERIC_POINTER(ei);
    d14 = NUMERIC_POINTER(bi);
    d15 = NUMERIC_POINTER(matv);
    d16 = NUMERIC_POINTER(arcoef);
    i7 = INTEGER_POINTER(morder);
    d17 = NUMERIC_POINTER(aics);

    F77_CALL(mulmarf)  (d1,i1,i2,d2,i3,d3,d4,d5,d6,d7,i4,d8,d9,i5,i6,d10,d11,d12,d13,d14,d15,d16,i7,d17);

    xmean = REAL(mean);
    xvar = REAL(var);
    xv = REAL(v);
    xaic = REAL(aic);
    xdaic = REAL(daic);
    xm = INTEGER(m);
    xaicm = REAL(aicm);
    xvm = REAL(vm);
    xnpr = INTEGER(npr);
    xjnd = INTEGER(jnd);
    xa = REAL(a);
    xrv = REAL(rv);
    xaicf = REAL(aicf);
    xei = REAL(ei);
    xbi = REAL(bi);
    xmatv = REAL(matv);
    xarcoef = REAL(arcoef);
    xmorder = INTEGER(morder);
    xaics = REAL(aics);

    for(i=0; i<nd; i++) xmean[i] = d3[i];
    for(i=0; i<nd; i++) xvar[i] = d4[i];
    for(i=0; i<lg1*nd; i++) xv[i] = d5[i];
    for(i=0; i<lg1*nd; i++) xaic[i] = d6[i];
    for(i=0; i<lg1*nd; i++) xdaic[i] = d7[i];
    for(i=0; i<nd; i++) xm[i] = i4[i];
    for(i=0; i<nd; i++) xaicm[i] = d8[i];
    for(i=0; i<nd; i++) xvm[i] = d9[i];
    for(i=0; i<nd; i++) xnpr[i] = i5[i];
    for(i=0; i<lg1*nd2; i++) xjnd[i] = i6[i];
    for(i=0; i<lg1*nd2; i++) xa[i] = d10[i];
    for(i=0; i<nd; i++) xrv[i] = d11[i];
    for(i=0; i<nd; i++) xaicf[i] = d12[i];
    for(i=0; i<nd2; i++) xei[i] = d13[i];
    for(i=0; i<nd2*lg; i++) xbi[i] = d14[i];
    for(i=0; i<nd2; i++) xmatv[i] = d15[i];
    for(i=0; i<nd2*lg; i++) xarcoef[i] = d16[i];
    *xmorder = *i7;
    *xaics = *d17;

    UNPROTECT(1);
    return ans; 
}




