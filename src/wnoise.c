#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(wnoisef) (int*, int*, double*, double*);

SEXP wnoise(SEXP len, SEXP ir, SEXP perr)
{
    double *d1,*d2;
    int *i1,*i2;

    SEXP ans =  R_NilValue,  wn = R_NilValue;
    double *xwn = NULL;
    int   i, k;

    i1 = INTEGER_POINTER(len);
    i2 = INTEGER_POINTER(ir);
    d1 = NUMERIC_POINTER(perr);

    k = (*i1) * (*i2);
    PROTECT(ans = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(ans, 0, wn = allocVector(REALSXP, k));

    d2 = NUMERIC_POINTER(wn);

    F77_CALL(wnoisef) (i1,i2,d1,d2);

    xwn = REAL(wn);

    for(i=0; i<k; i++) xwn[i] = d2[i];

    UNPROTECT(1);

    return ans;
}

