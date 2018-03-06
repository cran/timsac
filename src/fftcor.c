#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(fftcorf) (int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

SEXP FftcorC(SEXP ld, SEXP lag1, SEXP n, SEXP n2p, SEXP isw, SEXP x1, SEXP y1)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int *i1,*i2,*i3,*i4,*i5;

    SEXP ans =  R_NilValue, acov = R_NilValue, cov21 = R_NilValue,  cov12 = R_NilValue,  acor = R_NilValue, cor21 = R_NilValue, cor12 = R_NilValue, mean = R_NilValue;

    double *xacov, *xcov21, *xcov12, *xacor, *xcor21, *xcor12, *xmean = NULL;
    int    i, nn, lg1;

    i1 = INTEGER_POINTER(ld);
    i2 = INTEGER_POINTER(lag1);
    i3 = INTEGER_POINTER(n);
    i4 = INTEGER_POINTER(n2p);
    i5 = INTEGER_POINTER(isw);
    d1 = NUMERIC_POINTER(x1);
    d2 = NUMERIC_POINTER(y1);

    lg1 = *i2;
    nn = *i3;
    PROTECT(ans = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(ans, 0, acov = allocVector(REALSXP, nn*2));
    SET_VECTOR_ELT(ans, 1, cov21 = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 2, cov12 = allocVector(REALSXP, nn)); 
    SET_VECTOR_ELT(ans, 3, acor = allocVector(REALSXP, lg1*2));
    SET_VECTOR_ELT(ans, 4, cor21 = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 5, cor12 = allocVector(REALSXP, lg1)); 
    SET_VECTOR_ELT(ans, 6, mean = allocVector(REALSXP, 2));

    d3 = NUMERIC_POINTER(acov);
    d4 = NUMERIC_POINTER(cov21);
    d5 = NUMERIC_POINTER(cov12);
    d6 = NUMERIC_POINTER(acor);
    d7 = NUMERIC_POINTER(cor21);
    d8 = NUMERIC_POINTER(cor12);
    d9 = NUMERIC_POINTER(mean);

    F77_CALL(fftcorf) (i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9);

    xacov = REAL(acov);
    xcov21 = REAL(cov21);
    xcov12 = REAL(cov12);
    xacor = REAL(acor);
    xcor21 = REAL(cor21);
    xcor12 = REAL(cor12);
    xmean = REAL(mean);

    for(i=0; i<nn*2; i++) xacov[i] = d3[i];
    for(i=0; i<nn; i++) xcov21[i] = d4[i];
    for(i=0; i<nn; i++) xcov12[i] = d5[i];
    for(i=0; i<lg1*2; i++) xacor[i] = d6[i];
    for(i=0; i<lg1; i++) xcor21[i] = d7[i];
    for(i=0; i<lg1; i++) xcor12[i] = d8[i];
    for(i=0; i<2; i++) xmean[i] = d9[i];

    UNPROTECT(1);

    return ans;
}

