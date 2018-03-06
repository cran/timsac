
#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(perarsf)  (double*,int*,int*,int*,int*,double*,double*,int*,int*,double*,double*,double*,double*,double*,double*,int*);

SEXP PerarsC(SEXP y, SEXP n, SEXP ni, SEXP lag, SEXP ksw)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7;

    SEXP ans =  R_NilValue, mean = R_NilValue, var = R_NilValue, np = R_NilValue, jnd = R_NilValue, a = R_NilValue, aic = R_NilValue, b = R_NilValue, v = R_NilValue, c = R_NilValue, osd = R_NilValue, mo = R_NilValue;
    double *xmean, *xvar, *xa, *xaic, *xb, *xv, *xc, *xosd = NULL;
    int    *xnp, *xjnd, *xmo = NULL;
    int   i, nn, lg, ks, k;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(ni);
    i3 = INTEGER_POINTER(lag);
    i4 = INTEGER_POINTER(ksw);

    nn = *i2;
    lg = *i3;
    ks = *i4;
    k = ((lg+1)*nn+ks)*nn;
    PROTECT(ans = allocVector(VECSXP, 11));
    SET_VECTOR_ELT(ans, 0, mean = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, var = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, np = allocVector(INTSXP, nn)); 
    SET_VECTOR_ELT(ans, 3, jnd = allocVector(INTSXP, k));
    SET_VECTOR_ELT(ans, 4, a = allocVector(REALSXP, k)); 
    SET_VECTOR_ELT(ans, 5, aic = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 6, b = allocVector(REALSXP, nn*nn*lg));
    SET_VECTOR_ELT(ans, 7, v = allocVector(REALSXP, nn*nn)); 
    SET_VECTOR_ELT(ans, 8, c = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 9, osd = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 10, mo = allocVector(INTSXP, 1)); 

    d2 = NUMERIC_POINTER(mean);
    d3 = NUMERIC_POINTER(var);
    i5 = INTEGER_POINTER(np);
    i6 = INTEGER_POINTER(jnd);
    d4 = NUMERIC_POINTER(a);
    d5 = NUMERIC_POINTER(aic);
    d6 = NUMERIC_POINTER(b);
    d7 = NUMERIC_POINTER(v);
    d8 = NUMERIC_POINTER(c);
    d9 = NUMERIC_POINTER(osd);
    i7 = INTEGER_POINTER(mo);

    F77_CALL(perarsf) (d1,i1,i2,i3,i4,d2,d3,i5,i6,d4,d5,d6,d7,d8,d9,i7);

    xmean = REAL(mean);
    xvar = REAL(var);
    xnp = INTEGER(np);
    xjnd = INTEGER(jnd);
    xa = REAL(a);
    xaic = REAL(aic);
    xb = REAL(b);
    xv = REAL(v);
    xc = REAL(c);
    xosd = REAL(osd);
    xmo = INTEGER(mo);

    *xmean = *d2;
    *xvar = *d3;
    for(i=0; i<nn; i++) xnp[i] = i5[i];
    for(i=0; i<k; i++) xjnd[i] = i6[i];
    for(i=0; i<k; i++) xa[i] = d4[i];
    for(i=0; i<nn; i++) xaic[i] = d5[i];
    for(i=0; i<nn*nn*lg; i++) xb[i] = d6[i];
    for(i=0; i<nn*nn; i++) xv[i] = d7[i];
    for(i=0; i<nn; i++) xc[i] = d8[i];
    for(i=0; i<nn; i++) xosd[i] = d9[i];
    *xmo = *i7;

    UNPROTECT(1);

    return ans;
}


