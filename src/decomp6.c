

#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(decompf) (double*, int*, int*, double*, double*, double*, double*, double*, double*, int*, double*);

SEXP decomp(SEXP y, SEXP n, SEXP ipar, SEXP miss, SEXP omax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8;
    int *i1,*i2,*i3;

    SEXP ans = R_NilValue, trend = R_NilValue, seasonal = R_NilValue, ar = R_NilValue, trad = R_NilValue, noise = R_NilValue, para = R_NilValue;
    double *xtrend, *xseasonal, *xar, *xtrad, *xnoise, *xpara = NULL;
    int i, nn;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(ipar);
    i3 = INTEGER_POINTER(miss);
    d8 = NUMERIC_POINTER(omax);

    nn = *i1;
    PROTECT(ans = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, trend = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 1, seasonal = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 2, ar = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 3, trad = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 4, noise = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 5, para = allocVector(REALSXP, 26));

    d2 = NUMERIC_POINTER(trend);
    d3 = NUMERIC_POINTER(seasonal);
    d4 = NUMERIC_POINTER(ar);
    d5 = NUMERIC_POINTER(trad);
    d6 = NUMERIC_POINTER(noise);
    d7 = NUMERIC_POINTER(para);

    F77_CALL(decompf) (d1,i1,i2,d2,d3,d4,d5,d6,d7,i3,d8);

    xtrend = REAL(trend);
    xseasonal = REAL(seasonal);
    xar = REAL(ar);
    xtrad = REAL(trad);
    xnoise = REAL(noise);
    xpara = REAL(para);

    for(i=0; i<nn; i++) xtrend[i] = d2[i];
    for(i=0; i<nn; i++) xseasonal[i] = d3[i];
    for(i=0; i<nn; i++) xar[i] = d4[i];
    for(i=0; i<nn; i++) xtrad[i] = d5[i];
    for(i=0; i<nn; i++) xnoise[i] = d6[i];
    for(i=0; i<26; i++) xpara[i] = d7[i];

    UNPROTECT(1);

    return ans;
}



