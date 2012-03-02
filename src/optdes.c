#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(optdesf) (int*, int*, int*, int*, double*, double*, double*, double*, double*, double*);

SEXP optdes(SEXP ncon, SEXP nman, SEXP ns, SEXP order, SEXP q, SEXP r, SEXP osd, SEXP a, SEXP b)
{
    double *d1,*d2,*d3,*d4,*d5,*d6;
    int *i1,*i2,*i3,*i4;

    SEXP ans = R_NilValue, gain =  R_NilValue;
    double *xgain = NULL;
    int   i, nn;

    i1 = INTEGER_POINTER(ncon);
    i2 = INTEGER_POINTER(nman);
    i3 = INTEGER_POINTER(ns);
    i4 = INTEGER_POINTER(order);
    d1 = NUMERIC_POINTER(q);
    d2 = NUMERIC_POINTER(r);
    d3 = NUMERIC_POINTER(osd);
    d4 = NUMERIC_POINTER(a);
    d5 = NUMERIC_POINTER(b);

    nn = (*i1) * (*i2) * (*i4);
    PROTECT(ans = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(ans, 0, gain = allocVector(REALSXP, nn));

    d6 = NUMERIC_POINTER(gain);

    F77_CALL(optdesf) (i1,i2,i3,i4,d1,d2,d3,d4,d5,d6);

    xgain = REAL(gain);
    for(i=0; i<nn; i++) xgain[i] = d6[i];

    UNPROTECT(1);

    return ans;
}


