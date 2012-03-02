#include <R.h>
#include <Rdefines.h>
#include <R_ext/Complex.h>
#include "timsac.h"

   extern void F77_NAME(mulfrff) (int*, int*, int*, int*, int*, double*, Rcomplex*, double*, double*, double*, double*, double*, double*, double*);

SEXP mulfrf(SEXP nv, SEXP iovar, SEXP n, SEXP lag1, SEXP d, SEXP spec)
{
    double *d1, *d2, *d3, *d4, *d5, *d6, *d7, *d8;
    int *i1, *i2, *i3, *i4, *i5;
    Rcomplex *c1;
    int   i;


 SEXP ans =  R_NilValue, cospec = R_NilValue, freqr = R_NilValue, freqi = R_NilValue;
    SEXP gain = R_NilValue, phase = R_NilValue, pcoh = R_NilValue, errstat = R_NilValue, mcoh = R_NilValue;
    Rcomplex *xcospec = NULL;
    double *xfreqr, *xfreqi, *xgain, *xphase, *xpcoh, *xerrstat, *xmcoh = NULL;

    i1 = INTEGER_POINTER(nv);
    i2 = INTEGER_POINTER(iovar);
    i3 = INTEGER_POINTER(n);
    i4 = INTEGER_POINTER(lag1);
    i5 = INTEGER_POINTER(d);
    d1 = NUMERIC_POINTER(spec);

    int inv = *i1;
    int ilag1 = *i4;
    int id = *i5;

    PROTECT(ans = allocVector(VECSXP, 8));
    SET_VECTOR_ELT(ans, 0, cospec = allocVector(CPLXSXP, id*id*ilag1)); 
    SET_VECTOR_ELT(ans, 1, freqr = allocVector(REALSXP, inv*ilag1));
    SET_VECTOR_ELT(ans, 2, freqi = allocVector(REALSXP, inv*ilag1));
    SET_VECTOR_ELT(ans, 3, gain = allocVector(REALSXP, inv*ilag1));
    SET_VECTOR_ELT(ans, 4, phase = allocVector(REALSXP, inv*ilag1));
    SET_VECTOR_ELT(ans, 5, pcoh = allocVector(REALSXP, inv*ilag1));
    SET_VECTOR_ELT(ans, 6, errstat = allocVector(REALSXP, inv*ilag1));
    SET_VECTOR_ELT(ans, 7, mcoh = allocVector(REALSXP, ilag1));

    xcospec = COMPLEX(cospec);
    xfreqr = REAL(freqr);
    xfreqi = REAL(freqi);
    xgain = REAL(gain);
    xphase = REAL(phase);
    xpcoh = REAL(pcoh);
    xerrstat = REAL(errstat);
    xmcoh = REAL(mcoh);

    c1 = COMPLEX_POINTER(cospec);
    d2 = NUMERIC_POINTER(freqr);
    d3 = NUMERIC_POINTER(freqi);
    d4 = NUMERIC_POINTER(gain);
    d5 = NUMERIC_POINTER(phase);
    d6 = NUMERIC_POINTER(pcoh);
    d7 = NUMERIC_POINTER(errstat);
    d8 = NUMERIC_POINTER(mcoh);

    F77_CALL(mulfrff)(i1,i2,i3,i4,i5,d1,c1,d2,d3,d4,d5,d6,d7,d8);

    for( i=0; i<ilag1*id*id; i++)  xcospec[i] = c1[i];
    for( i=0; i<ilag1; i++) xfreqr[i] = d2[i];
    for( i=0; i<ilag1; i++) xfreqi[i] = d3[i];
    for( i=0; i<ilag1; i++) xgain[i] = d4[i];
    for( i=0; i<ilag1; i++) xphase[i] = d5[i];
    for( i=0; i<ilag1; i++) xpcoh[i] = d6[i];
    for( i=0; i<ilag1; i++) xerrstat[i] = d7[i];
    for( i=0; i<ilag1; i++) xmcoh[i] = d8[i];

    UNPROTECT(1);

    return ans;
}

