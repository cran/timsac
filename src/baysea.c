#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

void F77_NAME(bayseaf)(double*, int*, int*, double*, double*, double*, double*, double*, double*, double*,
 double*, double*, double*, double*, int*, double*, double*, double*, double*, int*, int*, int*);

SEXP baysea(SEXP y, SEXP ndata, SEXP forecast, SEXP ipara, SEXP para, SEXP arft, SEXP arfs, SEXP arfn, SEXP iart, SEXP iars, SEXP iarn)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16;
    int *i1,*i2,*i3,*i4,*i5,*i6;
    int i;

    SEXP ans =  R_NilValue, outlier = R_NilValue, dmoi = R_NilValue, trend = R_NilValue;
    SEXP season = R_NilValue, tdcmp = R_NilValue, irreg = R_NilValue, adjust = R_NilValue;
    SEXP est = R_NilValue, psds = R_NilValue, psdt = R_NilValue, avabic = R_NilValue;
    double *xoutlier, *xdmoi, *xtrend, *xseason, *xtdcmp, *xirreg, *xadjust, *xest, *xpsds, *xpsdt, *xavabic = NULL;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(ndata);
    i2 = INTEGER_POINTER(forecast);
    i3 = INTEGER_POINTER(ipara);
    d13 = NUMERIC_POINTER(para);
    d14 = NUMERIC_POINTER(arft);
    d15 = NUMERIC_POINTER(arfs);
    d16 = NUMERIC_POINTER(arfn);
    i4 = INTEGER_POINTER(iart);
    i5 = INTEGER_POINTER(iars);
    i6 = INTEGER_POINTER(iarn);	

    int mdata = *i1;
    int mfocast = *i2;
    int npf = mdata + mfocast;
    PROTECT(ans = allocVector(VECSXP, 11));
    SET_VECTOR_ELT(ans, 0, outlier = allocVector(REALSXP, mdata));
    SET_VECTOR_ELT(ans, 1, dmoi = allocVector(REALSXP, mdata));
    SET_VECTOR_ELT(ans, 2, trend = allocVector(REALSXP, npf));
    SET_VECTOR_ELT(ans, 3, season = allocVector(REALSXP, npf));
    SET_VECTOR_ELT(ans, 4, tdcmp = allocVector(REALSXP, npf));
    SET_VECTOR_ELT(ans, 5, irreg = allocVector(REALSXP, mdata));
    SET_VECTOR_ELT(ans, 6, adjust = allocVector(REALSXP, mdata));
    SET_VECTOR_ELT(ans, 7, est = allocVector(REALSXP, npf));
    SET_VECTOR_ELT(ans, 8, psds = allocVector(REALSXP, npf));
    SET_VECTOR_ELT(ans, 9, psdt = allocVector(REALSXP, npf));
    SET_VECTOR_ELT(ans, 10, avabic = allocVector(REALSXP, 1));

    d2 = NUMERIC_POINTER(outlier);
    d3 = NUMERIC_POINTER(dmoi);
    d4 = NUMERIC_POINTER(trend);
    d5 = NUMERIC_POINTER(season);
    d6 = NUMERIC_POINTER(tdcmp);
    d7 = NUMERIC_POINTER(irreg);
    d8 = NUMERIC_POINTER(adjust);
    d9 = NUMERIC_POINTER(est);
    d10 = NUMERIC_POINTER(psds);
    d11 = NUMERIC_POINTER(psdt);
    d12 = NUMERIC_POINTER(avabic);

    F77_CALL(bayseaf) (d1,i1,i2,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,i3,d13,d14,d15,d16,i4,i5,i6);

    xoutlier = REAL(outlier);
    xdmoi = REAL(dmoi);
    xtrend = REAL(trend);
    xseason = REAL(season);
    xtdcmp = REAL(tdcmp);
    xirreg = REAL(irreg);
    xadjust = REAL(adjust);
    xest = REAL(est);
    xpsds = REAL(psds);
    xpsdt = REAL(psdt);
    xavabic = REAL(avabic);

    int iout = i3[6];
    if( iout == 0 ) xoutlier = NULL;
    if( iout != 0 ) for(i=0; i<mdata; i++) xoutlier[i] = d2[i];
    for(i=0; i<mdata; i++) xdmoi[i] = d3[i];
    for(i=0; i<npf; i++) xtrend[i] = d4[i];
    for(i=0; i<npf; i++) xseason[i] = d5[i];
    int iyear = i3[10];
    if( iyear == 0 ) xtdcmp = NULL;
    if( iyear != 0 ) for(i=0; i<npf; i++) xtdcmp[i] = d6[i];
    for(i=0; i<mdata; i++) xirreg[i] = d7[i];
    for(i=0; i<mdata; i++) xadjust[i] = d8[i];
    for(i=0; i<npf; i++) xest[i] = d9[i];
    for(i=0; i<npf; i++) xpsds[i] = d10[i];
    for(i=0; i<npf; i++) xpsdt[i] = d11[i];
    *xavabic = *d12;



    UNPROTECT(1);

    return ans;
}
