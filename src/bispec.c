#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void  F77_NAME(bispecf)(int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

SEXP BispecC(SEXP n, SEXP lag, SEXP cv, SEXP tmnt)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int *i1, *i2;
    int i, lag1, lag12;

    SEXP ans =  R_NilValue, pspec1 = R_NilValue, pspec2 = R_NilValue, sig = R_NilValue, ch = R_NilValue, br = R_NilValue, bi = R_NilValue, rat = R_NilValue;
    double *xpspc1, *xpspc2, *xsig, *xch, *xbr, *xbi, *xrat = NULL;

    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(lag);
    d1 = NUMERIC_POINTER(cv);
    d2 = NUMERIC_POINTER(tmnt);

    lag1 = *i2 + 1;
    lag12 = lag1 * lag1;
    PROTECT(ans = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(ans, 0, pspec1 = allocVector(REALSXP, lag1));
    SET_VECTOR_ELT(ans, 1, pspec2 = allocVector(REALSXP, lag1));
    SET_VECTOR_ELT(ans, 2, sig = allocVector(REALSXP, lag1)); 
    SET_VECTOR_ELT(ans, 3, ch = allocVector(REALSXP, lag12));
    SET_VECTOR_ELT(ans, 4, br = allocVector(REALSXP, lag12));
    SET_VECTOR_ELT(ans, 5, bi = allocVector(REALSXP, lag12)); 
    SET_VECTOR_ELT(ans, 6, rat = allocVector(REALSXP, 1)); 

    d3 = NUMERIC_POINTER(pspec1);
    d4 = NUMERIC_POINTER(pspec2);
    d5 = NUMERIC_POINTER(sig);
    d6 = NUMERIC_POINTER(ch);
    d7 = NUMERIC_POINTER(br);
    d8 = NUMERIC_POINTER(bi);
    d9 = NUMERIC_POINTER(rat);

    F77_CALL(bispecf) (i1,i2,d1,d2,d3,d4,d5,d6,d7,d8,d9);

    xpspc1 = REAL(pspec1);
    xpspc2 = REAL(pspec2);
    xsig = REAL(sig);
    xch = REAL(ch);
    xbr = REAL(br);
    xbi = REAL(bi);
    xrat = REAL(rat);

    for(i=0; i<lag1; i++) xpspc1[i] = d3[i];
    for(i=0; i<lag1; i++) xpspc2[i] = d4[i];
    for(i=0; i<lag1; i++) xsig[i] = d5[i];
    for(i=0; i<lag12; i++) xch[i] = d6[i];
    for(i=0; i<lag12; i++) xbr[i] = d7[i];
    for(i=0; i<lag12; i++) xbi[i] = d8[i];
    *xrat = *d9;

    UNPROTECT(1);

    return ans;
}

