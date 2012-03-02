#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(sglfref) (int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

SEXP sglfre(SEXP invar, SEXP outvar, SEXP n, SEXP lag1, SEXP d, SEXP spec)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11;
    int    *i1,*i2,*i3,*i4,*i5;

    SEXP ans = R_NilValue, spec1 = R_NilValue, spec2 = R_NilValue, cspec = R_NilValue, qspec = R_NilValue,  gain = R_NilValue, coh = R_NilValue, freqr = R_NilValue, freqi = R_NilValue, err = R_NilValue, phase =  R_NilValue;

    double *xspec1, *xspec2, *xcspec, *xqspec, *xgain, *xcoh, *xfreqr, *xfreqi, *xerr, *xphase = NULL;
    int   i, lg1;

    i1 = INTEGER_POINTER(invar);
    i2 = INTEGER_POINTER(outvar);
    i3 = INTEGER_POINTER(n);
    i4 = INTEGER_POINTER(lag1);
    i5 = INTEGER_POINTER(d);
    d1 = NUMERIC_POINTER(spec);

    lg1 = *i4;
    PROTECT(ans = allocVector(VECSXP, 10));
    SET_VECTOR_ELT(ans, 0, spec1 = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 1, spec2 = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 2, cspec = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 3, qspec = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 4, gain = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 5, coh = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 6, freqr = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 7, freqi = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 8, err = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 9, phase = allocVector(REALSXP, lg1));

    d2 = NUMERIC_POINTER(spec1);
    d3 = NUMERIC_POINTER(spec2);
    d4 = NUMERIC_POINTER(cspec);
    d5 = NUMERIC_POINTER(qspec);
    d6 = NUMERIC_POINTER(gain);
    d7 = NUMERIC_POINTER(coh);
    d8 = NUMERIC_POINTER(freqr);
    d9 = NUMERIC_POINTER(freqi);
    d10 = NUMERIC_POINTER(err);
    d11 = NUMERIC_POINTER(phase);

    F77_CALL(sglfref) (i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11);

    xspec1 = REAL(spec1);
    xspec2 = REAL(spec2);
    xcspec = REAL(cspec);
    xqspec = REAL(qspec);
    xgain = REAL(gain);
    xcoh = REAL(coh);
    xfreqr = REAL(freqr);
    xfreqi = REAL(freqi);
    xerr = REAL(err);
    xphase = REAL(phase);

    for(i=0; i<lg1; i++) xspec1[i] = d2[i];
    for(i=0; i<lg1; i++) xspec2[i] = d3[i];
    for(i=0; i<lg1; i++) xcspec[i] = d4[i];
    for(i=0; i<lg1; i++) xqspec[i] = d5[i];
    for(i=0; i<lg1; i++) xgain[i] = d6[i];
    for(i=0; i<lg1; i++) xcoh[i] = d7[i];
    for(i=0; i<lg1; i++) xfreqr[i] = d8[i];
    for(i=0; i<lg1; i++) xfreqi[i] = d9[i];
    for(i=0; i<lg1; i++) xerr[i] = d10[i];
    for(i=0; i<lg1; i++) xphase[i] = d11[i];

    UNPROTECT(1);

    return ans;
}
