
#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(fpeautf) (int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*, double*);

SEXP FpeautC(SEXP morder, SEXP n, SEXP sd, SEXP cxx)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13;
    int *i1,*i2,*i3;

    SEXP ans =  R_NilValue, sig2 = R_NilValue, fpe = R_NilValue,  rfpe = R_NilValue,  parcor = R_NilValue, chi2 = R_NilValue, ofpe = R_NilValue, fpem = R_NilValue, rfpem = R_NilValue, mo = R_NilValue, sig2m = R_NilValue, a = R_NilValue, ao = R_NilValue;

    double *xsig2, *xfpe, *xrfpe, *xparcor, *xchi2, *xofpe, *xfpem, *xrfpem, *xsig2m, *xa, *xao = NULL;
    int    *xmo = NULL;
    int    i, m;

    i1 = INTEGER_POINTER(morder);
    i2 = INTEGER_POINTER(n);
    d1 = NUMERIC_POINTER(sd);
    d2 = NUMERIC_POINTER(cxx);

    m = *i1;
    PROTECT(ans = allocVector(VECSXP, 12));
    SET_VECTOR_ELT(ans, 0, sig2 = allocVector(REALSXP, m));
    SET_VECTOR_ELT(ans, 1, fpe = allocVector(REALSXP, m));
    SET_VECTOR_ELT(ans, 2, rfpe = allocVector(REALSXP, m)); 
    SET_VECTOR_ELT(ans, 3, parcor = allocVector(REALSXP, m));
    SET_VECTOR_ELT(ans, 4, chi2 = allocVector(REALSXP, m));
    SET_VECTOR_ELT(ans, 5, ofpe = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 6, fpem = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 7, rfpem = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 8, mo = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 9, sig2m = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 10, a = allocVector(REALSXP, m*m)); 
    SET_VECTOR_ELT(ans, 11, ao = allocVector(REALSXP, m));

    d3 = NUMERIC_POINTER(sig2);
    d4 = NUMERIC_POINTER(fpe);
    d5 = NUMERIC_POINTER(rfpe);
    d6 = NUMERIC_POINTER(parcor);
    d7 = NUMERIC_POINTER(chi2);
    d8 = NUMERIC_POINTER(ofpe);
    d9 = NUMERIC_POINTER(fpem);
    d10 = NUMERIC_POINTER(rfpem);
    i3 = INTEGER_POINTER(mo);
    d11= NUMERIC_POINTER(sig2m);
    d12 = NUMERIC_POINTER(a);
    d13 = NUMERIC_POINTER(ao);

    F77_CALL(fpeautf) (i1,i2,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,i3,d11,d12,d13);

    xsig2 = REAL(sig2);
    xfpe = REAL(fpe);
    xrfpe = REAL(rfpe);
    xparcor = REAL(parcor);
    xchi2 = REAL(chi2);
    xofpe = REAL(ofpe);
    xfpem = REAL(fpem);
    xrfpem = REAL(rfpem);
    xmo = INTEGER(mo);
    xsig2m = REAL(sig2m);
    xa = REAL(a);
    xao = REAL(ao);

    for(i=0; i<m; i++) xsig2[i] = d3[i];
    for(i=0; i<m; i++) xfpe[i] = d4[i];
    for(i=0; i<m; i++) xrfpe[i] = d5[i];
    for(i=0; i<m; i++) xparcor[i] = d6[i];
    for(i=0; i<m; i++) xchi2[i] = d7[i];
    *xofpe = *d8;
    *xfpem = *d9;
    *xrfpem = *d10;
    *xmo = *i3;
    *xsig2m = *d11;
    for(i=0; i<2; i++) xa[i] = d12[i];
    for(i=0; i<2; i++) xao[i] = d13[i];

    UNPROTECT(1);

    return ans;
}
