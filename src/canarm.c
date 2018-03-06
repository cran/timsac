#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(canarmf) (int*,int*,double*,double*,int*,double*,double*,double*,int*,double*,int*,int*,int*,double*,double*,double*,double*,int*,double*,double*,int*,int*,double*,int*,double*,int*,int*);

SEXP CanarmC( SEXP n, SEXP morder1, SEXP autcv, SEXP l1, SEXP mmax, SEXP nmax )
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13;

    SEXP ans =  R_NilValue, arcoef = R_NilValue, v = R_NilValue, aic = R_NilValue, oaic = R_NilValue, mo = R_NilValue, parcor = R_NilValue;
    SEXP nc = R_NilValue, m1 = R_NilValue, m2 = R_NilValue, w = R_NilValue, z = R_NilValue, Rs = R_NilValue, chi = R_NilValue, ndt = R_NilValue, dic = R_NilValue, dicm = R_NilValue, po = R_NilValue;
    SEXP k = R_NilValue, b = R_NilValue, l = R_NilValue, a = R_NilValue;

    double *xarcoef, *xv, *xaic, *xoaic, *xparcor, *xw, *xz, *xRs, *xchi, *xdic, *xdicm, *xb, *xa = NULL;
    int      *xmo, *xnc, *xm1, *xm2, *xndt, *xpo, *xk, *xl = NULL;
    int      i, mm, mm2, mm3, nn;


    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(morder1);
    d1 = NUMERIC_POINTER(autcv);
    i3 = INTEGER_POINTER(l1);
    i12 = INTEGER_POINTER(mmax);
    i13 = INTEGER_POINTER(nmax);

    mm = *i12;
    mm2 = mm * mm;
    mm3 = mm2 * mm;
    nn = *i13;
    PROTECT(ans = allocVector(VECSXP, 21));
    SET_VECTOR_ELT(ans, 0, arcoef = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 1, v = allocVector(REALSXP, mm+1));
    SET_VECTOR_ELT(ans, 2, aic = allocVector(REALSXP, mm+1)); 
    SET_VECTOR_ELT(ans, 3, oaic = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 4, mo = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 5, parcor = allocVector(REALSXP, mm));
    SET_VECTOR_ELT(ans, 6, nc = allocVector(INTSXP, 1)); 
    SET_VECTOR_ELT(ans, 7, m1 = allocVector(INTSXP, mm));
    SET_VECTOR_ELT(ans, 8, m2 = allocVector(INTSXP, mm));
    SET_VECTOR_ELT(ans, 9, w = allocVector(REALSXP, mm3)); 
    SET_VECTOR_ELT(ans, 10, z = allocVector(REALSXP, mm2));
    SET_VECTOR_ELT(ans, 11, Rs = allocVector(REALSXP, mm2)); 
    SET_VECTOR_ELT(ans, 12, chi = allocVector(REALSXP, mm2)); 
    SET_VECTOR_ELT(ans, 13, ndt = allocVector(INTSXP, mm2)); 
    SET_VECTOR_ELT(ans, 14, dic = allocVector(REALSXP, mm2)); 
    SET_VECTOR_ELT(ans, 15, dicm = allocVector(REALSXP, mm)); 
    SET_VECTOR_ELT(ans, 16, po = allocVector(INTSXP, mm));
    SET_VECTOR_ELT(ans, 17, k = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 18, b = allocVector(REALSXP, mm)); 
    SET_VECTOR_ELT(ans, 19, l = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 20, a = allocVector(REALSXP, mm)); 

    d2 = NUMERIC_POINTER(arcoef);
    d3 = NUMERIC_POINTER(v);
    d4 = NUMERIC_POINTER(aic);
    d5 = NUMERIC_POINTER(oaic);
    i4 = INTEGER_POINTER(mo);
    d6 = NUMERIC_POINTER(parcor);
    i5 = INTEGER_POINTER(nc);
    i6 = INTEGER_POINTER(m1);
    i7 = INTEGER_POINTER(m2);
    d7 = NUMERIC_POINTER(w);
    d8 = NUMERIC_POINTER(z);
    d9 = NUMERIC_POINTER(Rs);
    d10 = NUMERIC_POINTER(chi);
    i8 = INTEGER_POINTER(ndt);
    d11 = NUMERIC_POINTER(dic);
    d12 = NUMERIC_POINTER(dicm);
    i9 = INTEGER_POINTER(po);
    i10 = INTEGER_POINTER(k);
    d13 = NUMERIC_POINTER(b);
    i11 = INTEGER_POINTER(l);
    d14 = NUMERIC_POINTER(a);

    F77_CALL(canarmf) (i1,i2,d1,d2,i3,d3,d4,d5,i4,d6,i5,i6,i7,d7,d8,d9,d10,i8,d11,d12,i9,i10,d13,i11,d14,i12,i13);

    xarcoef = REAL(arcoef);
    xv = REAL(v);
    xaic = REAL(aic);
    xoaic = REAL(oaic);
    xmo = INTEGER(mo);
    xparcor = REAL(parcor);
    xnc = INTEGER(nc);
    xm1 = INTEGER(m1);
    xm2 = INTEGER(m2);
    xw = REAL(w);
    xz = REAL(z);
    xRs = REAL(Rs);
    xchi = REAL(chi);
    xndt = INTEGER(ndt);
    xdic = REAL(dic);
    xdicm = REAL(dicm);
    xpo = INTEGER(po);
    xk = INTEGER(k);
    xb = REAL(b);
    xl = INTEGER(l);
    xa = REAL(a);

    for(i=0; i<nn; i++) xarcoef[i] = d2[i];
    for(i=0; i<mm+1; i++) xv[i] = d3[i];
    for(i=0; i<mm+1; i++) xaic[i] = d4[i];
    *xoaic = *d5;
    *xmo = *i4;
    for(i=0; i<mm; i++) xparcor[i] = d6[i];
    *xnc = *i5;
    for(i=0; i<mm; i++) xm1[i] = i6[i];
    for(i=0; i<mm; i++) xm2[i] = i7[i];
    for(i=0; i<mm3; i++) xw[i] = d7[i];
    for(i=0; i<mm2; i++) xz[i] = d8[i];
    for(i=0; i<mm2; i++) xRs[i] = d9[i];
    for(i=0; i<mm2; i++) xchi[i] = d10[i];
    for(i=0; i<mm2; i++) xndt[i] = i8[i];
    for(i=0; i<mm2; i++) xdic[i] = d11[i];
    for(i=0; i<mm; i++) xdicm[i] = d12[i];
    for(i=0; i<mm; i++) xpo[i] = i9[i];
    *xk = *i10;
    for(i=0; i<mm; i++) xb[i] = d13[i];
    *xl = *i11;
    for(i=0; i<mm; i++) xa[i] = d14[i];

    UNPROTECT(1);
    return ans; 
}

