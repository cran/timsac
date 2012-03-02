#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(spgrhf) (double*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*);

SEXP spgrh(SEXP y, SEXP n, SEXP lag1, SEXP ifpl1, SEXP mode, SEXP period)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8;
    int *i1,*i2,*i3,*i4,*i5,*i6;

    SEXP ans = R_NilValue, cxx = R_NilValue, cn = R_NilValue, xmean = R_NilValue, sd = R_NilValue, aic = R_NilValue, parcor = R_NilValue, pxx = R_NilValue, ier = R_NilValue;
    double *xcxx, *xcn, *xxmean, *xsd, *xaic, *xparcor, *xpxx = NULL;
    int   *xier = NULL;
    int   i, lg1, ip1;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(lag1);
    i3 = INTEGER_POINTER(ifpl1);
    i4 = INTEGER_POINTER(mode);
    i5 = INTEGER_POINTER(period);

    lg1 = *i2;
    ip1 = *i3;
    PROTECT(ans = allocVector(VECSXP, 8));
    SET_VECTOR_ELT(ans, 0, cxx = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 1, cn = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 2, xmean = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 3, sd = allocVector(REALSXP, ip1));
    SET_VECTOR_ELT(ans, 4, aic = allocVector(REALSXP, ip1));
    SET_VECTOR_ELT(ans, 5, parcor = allocVector(REALSXP, ip1-1));
    SET_VECTOR_ELT(ans, 6, pxx = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 7, ier = allocVector(INTSXP, 1));

    d2 = NUMERIC_POINTER(cxx);
    d3 = NUMERIC_POINTER(cn);
    d4 = NUMERIC_POINTER(xmean);
    d5 = NUMERIC_POINTER(sd);
    d6 = NUMERIC_POINTER(aic);
    d7 = NUMERIC_POINTER(parcor);
    d8 = NUMERIC_POINTER(pxx);
    i6 = INTEGER_POINTER(ier);

    F77_CALL(spgrhf) (d1,i1,i2,i3,i4,i5,d2,d3,d4,d5,d6,d7,d8,i6);

    xcxx = REAL(cxx);
    xcn = REAL(cn);
    xxmean = REAL(xmean);
    xsd = REAL(sd);
    xaic = REAL(aic);
    xparcor = REAL(parcor);
    xpxx = REAL(pxx);
    xier = INTEGER(ier);

    for(i=0; i<lg1; i++) xcxx[i] = d2[i];
    for(i=0; i<lg1; i++) xcn[i] = d3[i];
    *xxmean = *d4;
    for(i=0; i<ip1; i++) xsd[i] = d5[i];
    for(i=0; i<ip1; i++) xaic[i] = d6[i];
    for(i=0; i<(ip1-1); i++) xparcor[i] = d7[i];
    for(i=0; i<lg1; i++) xpxx[i] = d8[i];
    *xier = *i6;

    UNPROTECT(1);

    return ans;
}
