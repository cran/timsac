
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Complex.h>
#include "timsac.h"

   extern void F77_NAME(mulrspf) (int*, int*, int*, int*, double*, double*, double*, Rcomplex*, double*);

SEXP MulrspC(SEXP h, SEXP l, SEXP d, SEXP k, SEXP cov, SEXP arcoef, SEXP macoef)
{
    double *d1, *d2, *d3, *d4;
    int    *i1, *i2, *i3, *i4;
    Rcomplex *c1;

    SEXP ans =  R_NilValue, rspec = R_NilValue, scoh = R_NilValue;
    Rcomplex *xrspec = NULL;
    double *xscoh = NULL;
    int   i, h1, nd2;

    i1 = INTEGER_POINTER(h);
    i2 = INTEGER_POINTER(l);
    i3 = INTEGER_POINTER(d);
    i4 = INTEGER_POINTER(k);
    d1 = NUMERIC_POINTER(cov);
    d2 = NUMERIC_POINTER(arcoef);
    d3 = NUMERIC_POINTER(macoef);

    h1 = *i1 + 1;
    nd2 = (*i3) * (*i3);
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, rspec = allocVector(CPLXSXP, nd2*h1)); 
    SET_VECTOR_ELT(ans, 1, scoh = allocVector(REALSXP, nd2*h1));

    xrspec = COMPLEX(rspec);
    xscoh = REAL(scoh);

    c1 = COMPLEX_POINTER(rspec);
    d4 = NUMERIC_POINTER(scoh);

    F77_CALL(mulrspf) (i1,i2,i3,i4,d1,d2,d3,c1,d4);

    for(i=0; i<nd2*h1; i++) xrspec[i] = c1[i];
    for(i=0; i<nd2*h1; i++) xscoh[i] = d4[i];

    UNPROTECT(1);

    return ans;
}
