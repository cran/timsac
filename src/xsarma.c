#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(xsarmaf) (double*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

SEXP XsarmaC(SEXP y, SEXP n, SEXP p, SEXP q, SEXP p01)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
    int    *i1,*i2,*i3;

    SEXP ans =  R_NilValue, g1 = R_NilValue, tl1 = R_NilValue, p02 = R_NilValue, g2 = R_NilValue, alpar = R_NilValue, alpma = R_NilValue, tl2 = R_NilValue, sigma2 = R_NilValue;

    double *xg1, *xtl1, *xp02, *xg2, *xalpar, *xalpma, *xtl2, *xsigma2 = NULL;
    int   i, np, nq, npq;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(p);
    i3 = INTEGER_POINTER(q);
    d2 = NUMERIC_POINTER(p01);

    np = *i2;
    nq = *i3;
    npq = np + nq;
    PROTECT(ans = allocVector(VECSXP, 8));
    SET_VECTOR_ELT(ans, 0, g1 = allocVector(REALSXP, npq));
    SET_VECTOR_ELT(ans, 1, tl1 = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, p02 = allocVector(REALSXP, npq));
    SET_VECTOR_ELT(ans, 3, g2 = allocVector(REALSXP, npq));
    SET_VECTOR_ELT(ans, 4, alpar = allocVector(REALSXP, np));
    SET_VECTOR_ELT(ans, 5, alpma = allocVector(REALSXP, nq));
    SET_VECTOR_ELT(ans, 6, tl2 = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 7, sigma2 = allocVector(REALSXP, 1));

    d3 = NUMERIC_POINTER(g1);
    d4 = NUMERIC_POINTER(tl1);
    d5 = NUMERIC_POINTER(p02);
    d6 = NUMERIC_POINTER(g2);
    d7 = NUMERIC_POINTER(alpar);
    d8 = NUMERIC_POINTER(alpma);
    d9 = NUMERIC_POINTER(tl2);
    d10 = NUMERIC_POINTER(sigma2);

    F77_CALL(xsarmaf) (d1,i1,i2,i3,d2,d3,d4,d5,d6,d7,d8,d9,d10);

    xg1 = REAL(g1);
    xtl1 = REAL(tl1);
    xp02 = REAL(p02);
    xg2 = REAL(g2);
    xalpar = REAL(alpar);
    xalpma = REAL(alpma);
    xtl2 = REAL(tl2);
    xsigma2 = REAL(sigma2);

    for(i=0; i<npq; i++) xg1[i] = d3[i];
    *xtl1 = *d4;
    for(i=0; i<npq; i++) xp02[i] = d5[i];
    for(i=0; i<npq; i++) xg2[i] = d6[i];
    for(i=0; i<np; i++) xalpar[i] = d7[i];
    for(i=0; i<nq; i++) xalpma[i] = d8[i];
    *xtl2 = *d9;
    *xsigma2 = *d10;

    UNPROTECT(1);

    return ans;
}

