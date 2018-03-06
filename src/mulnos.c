
#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(mulnosf) (int*, int*, int*, double*, double*, double*, double*, double*);

SEXP MulnosC(SEXP h, SEXP l, SEXP ip, SEXP sd, SEXP arcoef)
{
    double *d1,*d2,*d3,*d4,*d5;
    int    *i1,*i2,*i3;

    SEXP ans =  R_NilValue,  nsd = R_NilValue, drpc = R_NilValue, irpc = R_NilValue;
    double *xnsd, *xdrpc, *xirpc = NULL;
    int   i, ip2, h1;

    i1 = INTEGER_POINTER(h);
    i2 = INTEGER_POINTER(l);
    i3 = INTEGER_POINTER(ip);
    d1 = NUMERIC_POINTER(sd);
    d2 = NUMERIC_POINTER(arcoef);

    h1 = *i1 + 1;
    ip2 = (*i3) * (*i3);
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, nsd = allocVector(REALSXP, ip2));
    SET_VECTOR_ELT(ans, 1, drpc = allocVector(REALSXP, ip2*h1));
    SET_VECTOR_ELT(ans, 2, irpc = allocVector(REALSXP, ip2*h1)); 

    d3 = NUMERIC_POINTER(nsd);
    d4 = NUMERIC_POINTER(drpc);
    d5 = NUMERIC_POINTER(irpc);

    F77_CALL(mulnosf) (i1,i2,i3,d1,d2,d3,d4,d5);

    xnsd = REAL(nsd);
    xdrpc = REAL(drpc);
    xirpc = REAL(irpc);

    for(i=0; i<ip2; i++) xnsd[i] = d3[i];
    for(i=0; i<ip2*h1; i++) xdrpc[i] = d4[i];
    for(i=0; i<ip2*h1; i++) xirpc[i] = d5[i];

    UNPROTECT(1);

    return ans;
}

