#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(fpec7f) (int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, int*, double*, double*, double*, double*, double*);

SEXP Fpec7C(SEXP n, SEXP morder, SEXP ncon, SEXP ip, SEXP d, SEXP inw, SEXP cov)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7;

    SEXP ans =  R_NilValue, ccv = R_NilValue, fpec = R_NilValue,  rfpec = R_NilValue,  aic = R_NilValue, mo = R_NilValue, fpecm = R_NilValue, rfpecm = R_NilValue, aicm = R_NilValue, perr = R_NilValue, arcoef = R_NilValue;

    double *xccv, *xfpec, *xrfpec, *xaic, *xfpecm, *xrfpecm, *xaicm, *xperr, *xarcoef = NULL;
    int    *xmo = NULL;
    int    i, m1, nc, jp;

    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(morder);
    i3 = INTEGER_POINTER(ncon);
    i4 = INTEGER_POINTER(ip);
    i5 = INTEGER_POINTER(d);
    i6 = INTEGER_POINTER(inw);
    d1 = NUMERIC_POINTER(cov);

    m1 = *i1 + 1;
    nc = *i3;
    jp = *i4;
    PROTECT(ans = allocVector(VECSXP, 10));
    SET_VECTOR_ELT(ans, 0, ccv = allocVector(REALSXP, m1*jp*jp));
    SET_VECTOR_ELT(ans, 1, fpec = allocVector(REALSXP, m1));
    SET_VECTOR_ELT(ans, 2, rfpec = allocVector(REALSXP, m1)); 
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, m1));
    SET_VECTOR_ELT(ans, 4, mo = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 5, fpecm = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 6, rfpecm = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 7, aicm = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 8, perr = allocVector(REALSXP, nc*nc)); 
    SET_VECTOR_ELT(ans, 9, arcoef = allocVector(REALSXP, m1*nc*jp));

    d2 = NUMERIC_POINTER(ccv);
    d3 = NUMERIC_POINTER(fpec);
    d4 = NUMERIC_POINTER(rfpec);
    d5 = NUMERIC_POINTER(aic);
    i7 = INTEGER_POINTER(mo);
    d6 = NUMERIC_POINTER(fpecm);
    d7 = NUMERIC_POINTER(rfpecm);
    d8 = NUMERIC_POINTER(aicm);
    d9 = NUMERIC_POINTER(perr);
    d10 = NUMERIC_POINTER(arcoef);

    F77_CALL(fpec7f) (i1,i2,i3,i4,i5,i6,d1,d2,d3,d4,d5,i7,d6,d7,d8,d9,d10);

    xccv = REAL(ccv);
    xfpec = REAL(fpec);
    xrfpec = REAL(rfpec);
    xaic = REAL(aic);
    xmo = INTEGER(mo);
    xfpecm = REAL(fpecm);
    xrfpecm = REAL(rfpecm);
    xaicm = REAL(aicm);
    xperr = REAL(perr);
    xarcoef = REAL(arcoef);

    for(i=0; i<m1*jp*jp; i++) xccv[i] = d2[i];
    for(i=0; i<m1; i++) xfpec[i] = d3[i];
    for(i=0; i<m1; i++) xrfpec[i] = d4[i];
    for(i=0; i<m1; i++) xaic[i] = d5[i];
    *xmo = *i7;
    *xfpecm = *d6;
    *xrfpecm = *d7;
    *xaicm = *d8;
    for(i=0; i<nc*nc; i++) xperr[i] = d9[i];
    for(i=0; i<m1*nc*jp; i++) xarcoef[i] = d10[i];

    UNPROTECT(1);

    return ans;
}

