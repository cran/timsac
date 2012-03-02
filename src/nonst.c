#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(nonstf) (int*,int*,double*,int*,int*,int*,double*,double*,double*,double*,double*,int*,int*,double*);

SEXP nonst( SEXP n, SEXP span, SEXP y, SEXP ns, SEXP morder )
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7;

    SEXP ans =  R_NilValue, p = R_NilValue, coef = R_NilValue, v = R_NilValue, aic = R_NilValue, daic21 = R_NilValue, daic = R_NilValue, ks = R_NilValue, ke = R_NilValue, pspec = R_NilValue;

    double *xcoef, *xv, *xaic, *xdaic21, *xdaic, *xpspec = NULL;
    int      *xp, *xks, *xke = NULL;
    int      i, nn, mo;

    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(span);
    d1 = NUMERIC_POINTER(y);
    i3 = INTEGER_POINTER(ns);
    i4 = INTEGER_POINTER(morder);

    nn = *i3;
    mo = *i4;
    PROTECT(ans = allocVector(VECSXP, 9));
    SET_VECTOR_ELT(ans, 0, p = allocVector(INTSXP, nn));
    SET_VECTOR_ELT(ans, 1, coef = allocVector(REALSXP, nn*mo));
    SET_VECTOR_ELT(ans, 2, v = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, nn)); 
    SET_VECTOR_ELT(ans, 4, daic21 = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 5, daic = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 6, ks = allocVector(INTSXP, nn)); 
    SET_VECTOR_ELT(ans, 7, ke = allocVector(INTSXP, nn));
    SET_VECTOR_ELT(ans, 8, pspec = allocVector(REALSXP, 121*nn)); 

    i5 = INTEGER_POINTER(p);
    d2 = NUMERIC_POINTER(coef);
    d3 = NUMERIC_POINTER(v);
    d4 = NUMERIC_POINTER(aic);
    d5 = NUMERIC_POINTER(daic21);
    d6 = NUMERIC_POINTER(daic);
    i6 = INTEGER_POINTER(ks);
    i7 = INTEGER_POINTER(ke);
    d7 = NUMERIC_POINTER(pspec);

    F77_CALL(nonstf) (i1,i2,d1,i3,i4,i5,d2,d3,d4,d5,d6,i6,i7,d7);

    xp = INTEGER(p);
    xcoef = REAL(coef);
    xv = REAL(v);
    xaic = REAL(aic);
    xdaic21 = REAL(daic21);
    xdaic = REAL(daic);
    xks = INTEGER(ks);
    xke = INTEGER(ke);
    xpspec = REAL(pspec);

    for(i=0; i<nn; i++) xp[i] = i5[i];
    for(i=0; i<nn*mo; i++) xcoef[i] = d2[i];
    for(i=0; i<nn; i++) xv[i] = d3[i];
    for(i=0; i<nn; i++) xaic[i] = d4[i];
    for(i=0; i<nn; i++) xdaic21[i] = d5[i];
    for(i=0; i<nn; i++) xdaic[i] = d6[i];
    for(i=0; i<nn; i++) xks[i] = i6[i];
    for(i=0; i<nn; i++) xke[i] = i7[i];
    for(i=0; i<121*nn; i++) xpspec[i] = d7[i];

    UNPROTECT(1);
    return ans; 
}


