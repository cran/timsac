#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(raspecf) (int*, int*, int*, double*, double*, double*, double*);

SEXP RaspecC(SEXP h, SEXP l, SEXP k, SEXP var, SEXP arcoef, SEXP macoef)
{
    double *d1,*d2,*d3,*d4;
    int    *i1,*i2,*i3;

    SEXP ans =  R_NilValue,  rspec = R_NilValue;
    double *xrspec = NULL;
    int   i, h1;

    i1 = INTEGER_POINTER(h);
    i2 = INTEGER_POINTER(l);
    i3 = INTEGER_POINTER(k);
    d1 = NUMERIC_POINTER(var);
    d2 = NUMERIC_POINTER(arcoef);
    d3 = NUMERIC_POINTER(macoef);

    h1 = *i1 + 1;
    PROTECT(ans = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(ans, 0, rspec = allocVector(REALSXP, h1));

    d4 = NUMERIC_POINTER(rspec);

    F77_CALL(raspecf) (i1,i2,i3,d1,d2,d3,d4);

    xrspec = REAL(rspec);

    for(i=0; i<h1; i++) xrspec[i] = d4[i];

    UNPROTECT(1);

    return ans;
}
