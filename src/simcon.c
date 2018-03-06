#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(simconf)  (int*,int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

SEXP SimconC(SEXP d, SEXP k, SEXP span, SEXP len, SEXP r, SEXP arcoef, SEXP impuls, SEXP v, SEXP weight)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
    int *i1,*i2,*i3,*i4,*i5;

    SEXP ans = R_NilValue, bc = R_NilValue, bd = R_NilValue, g = R_NilValue, av = R_NilValue,  si = R_NilValue, s2 = R_NilValue;
    double *xbc, *xbd, *xg, *xav, *xsi, *xs2 = NULL;
    int   i, nd, nk, nr, nkd;

    i1 = INTEGER_POINTER(d);
    i2 = INTEGER_POINTER(k);
    i3 = INTEGER_POINTER(span);
    i4 = INTEGER_POINTER(len);
    i5 = INTEGER_POINTER(r);
    d1 = NUMERIC_POINTER(arcoef);
    d2 = NUMERIC_POINTER(impuls);
    d3 = NUMERIC_POINTER(v);
    d4 = NUMERIC_POINTER(weight);

    nd = *i1;
    nk = *i2;
    nr = *i5;
    nkd = nk * nd;
    PROTECT(ans = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, bc = allocVector(REALSXP, nr*nkd));
    SET_VECTOR_ELT(ans, 1, bd = allocVector(REALSXP, nkd*(nd-nr)));
    SET_VECTOR_ELT(ans, 2, g = allocVector(REALSXP, nr*nkd));
    SET_VECTOR_ELT(ans, 3, av = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 4, si = allocVector(REALSXP, nd));
    SET_VECTOR_ELT(ans, 5, s2 = allocVector(REALSXP, nd));

    d5 = NUMERIC_POINTER(bc);
    d6 = NUMERIC_POINTER(bd);
    d7 = NUMERIC_POINTER(g);
    d8 = NUMERIC_POINTER(av);
    d9 = NUMERIC_POINTER(si);
    d10 = NUMERIC_POINTER(s2);

    F77_CALL(simconf) (i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10);

    xbc = REAL(bc);
    xbd = REAL(bd);
    xg = REAL(g);
    xav = REAL(av);
    xsi = REAL(si);
    xs2 = REAL(s2);

    for(i=0; i<nr*nkd; i++) xbc[i] = d5[i];
    for(i=0; i<nkd*(nd-nr); i++) xbd[i] = d6[i];
    for(i=0; i<nr*nkd; i++) xg[i] = d7[i];
    for(i=0; i<nd; i++) xav[i] = d8[i];
    for(i=0; i<nd; i++) xsi[i] = d9[i];
    for(i=0; i<nd; i++) xs2[i] = d10[i];

    UNPROTECT(1);

    return ans;
}
