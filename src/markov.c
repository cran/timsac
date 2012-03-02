#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(markovf) (int*,int*,int*,double*,int*,int*,int*,double*,double*,int*,int*,int*,int*,int*,int*,double*,double*,double*,double*,double*,int*,double*,double*,double*,double*,double*,int*,int*,int*,int*);

SEXP markov( SEXP n, SEXP lag1, SEXP d, SEXP cov, SEXP k, SEXP nh, SEXP nvf, SEXP vectF, SEXP matGi, SEXP icont, SEXP mj3, SEXP mj4, SEXP mj6, SEXP mj7 )
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13,*i14,*i15,*i16,*i17;

    SEXP ans =  R_NilValue, id = R_NilValue, ir = R_NilValue, ij = R_NilValue, ik = R_NilValue, ngr = R_NilValue, gr =  R_NilValue, a1 = R_NilValue, a = R_NilValue, b = R_NilValue, vd = R_NilValue, iqm = R_NilValue, bm = R_NilValue, au = R_NilValue, zz = R_NilValue, v = R_NilValue, aic = R_NilValue;

    double *xgr, *xa1, *xa, *xb, *xvd, *xbm, *xau, *xzz, *xv, *xaic = NULL;
    int      *xid, *xir, *xij, *xik, *xngr, *xiqm = NULL;
    int      i, dd, kk, m4, m7;

    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(lag1);
    i3 = INTEGER_POINTER(d);
    d1 = NUMERIC_POINTER(cov);
    i4 = INTEGER_POINTER(k);
    i5 = INTEGER_POINTER(nh);
    i6 = INTEGER_POINTER(nvf);
    d2 = NUMERIC_POINTER(vectF);
    d3 = NUMERIC_POINTER(matGi);
    i7 = INTEGER_POINTER(icont);
    i14 = INTEGER_POINTER(mj3);
    i15 = INTEGER_POINTER(mj4);
    i16 = INTEGER_POINTER(mj6);
    i17 = INTEGER_POINTER(mj7);

    dd = *i3;
    kk = *i4;
    m4 = *i15;
    m7 = *i17;
    PROTECT(ans = allocVector(VECSXP, 16));
    SET_VECTOR_ELT(ans, 0, id = allocVector(INTSXP, kk));
    SET_VECTOR_ELT(ans, 1, ir = allocVector(INTSXP, kk));
    SET_VECTOR_ELT(ans, 2, ij = allocVector(INTSXP, dd));
    SET_VECTOR_ELT(ans, 3, ik = allocVector(INTSXP, dd));
    SET_VECTOR_ELT(ans, 4, ngr = allocVector(INTSXP, 1)); 
    SET_VECTOR_ELT(ans, 5, gr = allocVector(REALSXP, m4));
    SET_VECTOR_ELT(ans, 6, a1 = allocVector(REALSXP, kk*kk));
    SET_VECTOR_ELT(ans, 7, a = allocVector(REALSXP, kk*kk)); 
    SET_VECTOR_ELT(ans, 8, b = allocVector(REALSXP, kk*dd));
    SET_VECTOR_ELT(ans, 9, vd = allocVector(REALSXP, m4*m4));
    SET_VECTOR_ELT(ans, 10, iqm = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 11, bm = allocVector(REALSXP, dd*dd*m7)); 
    SET_VECTOR_ELT(ans, 12, au = allocVector(REALSXP, dd*dd*m7));
    SET_VECTOR_ELT(ans, 13, zz = allocVector(REALSXP, dd*dd*m7)); 
    SET_VECTOR_ELT(ans, 14, v = allocVector(REALSXP, dd*dd));
    SET_VECTOR_ELT(ans, 15, aic = allocVector(REALSXP, 1)); 

    i8 = INTEGER_POINTER(id);
    i9 = INTEGER_POINTER(ir);
    i10 = INTEGER_POINTER(ij);
    i11 = INTEGER_POINTER(ik);
    i12 = INTEGER_POINTER(ngr);
    d4 = NUMERIC_POINTER(gr);
    d5 = NUMERIC_POINTER(a1);
    d6 = NUMERIC_POINTER(a);
    d7 = NUMERIC_POINTER(b);
    d8 = NUMERIC_POINTER(vd);
    i13 = INTEGER_POINTER(iqm);
    d9 = NUMERIC_POINTER(bm);
    d10 = NUMERIC_POINTER(au);
    d11= NUMERIC_POINTER(zz);
    d12= NUMERIC_POINTER(v);
    d13 = NUMERIC_POINTER(aic);

    F77_CALL(markovf) (i1,i2,i3,d1,i4,i5,i6,d2,d3,i7,i8,i9,i10,i11,i12,d4,d5,d6,d7,d8,i13,d9,d10,d11,d12,d13,i14,i15,i16,i17);

    xid = INTEGER(id);
    xir = INTEGER(ir);
    xij = INTEGER(ij);
    xik = INTEGER(ik);
    xngr = INTEGER(ngr);
    xgr = REAL(gr);
    xa1 = REAL(a1);
    xa = REAL(a);
    xb = REAL(b);
    xvd = REAL(vd);
    xiqm = INTEGER(iqm);
    xbm = REAL(bm);
    xau = REAL(au);
    xzz = REAL(zz);
    xv = REAL(v);
    xaic = REAL(aic);

    for(i=0; i<kk; i++) xid[i] = i8[i];
    for(i=0; i<kk; i++) xir[i] = i9[i];
    for(i=0; i<dd; i++) xij[i] = i10[i];
    for(i=0; i<dd; i++) xik[i] = i11[i];
    *xngr = *i12;
    for(i=0; i<m4; i++) xgr[i] = d4[i];
    for(i=0; i<kk*kk; i++) xa1[i] = d5[i];
    for(i=0; i<kk*kk; i++) xa[i] = d6[i];
    for(i=0; i<kk*dd; i++) xb[i] = d7[i];
    for(i=0; i<m4*m4; i++) xvd[i] = d8[i];
    *xiqm = *i13;
    for(i=0; i<dd*dd*m7; i++) xbm[i] = d9[i];
    for(i=0; i<dd*dd*m7; i++) xau[i] = d10[i];
    for(i=0; i<dd*dd*m7; i++) xzz[i] = d11[i];
    for(i=0; i<dd*dd; i++) xv[i] = d12[i];
    *xaic = *d13;

    UNPROTECT(1);
    return ans; 
}



