#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(canocaf) (int*,int*,int*,int*,int*,double*,int*,double*,double*,int*,double*,double*,int*,int*,int*,double*,double*,double*,double*,int*,double*,double*,int*,double*,int*,int*,double*,int*,double*,int*,int*,int*);

SEXP canoca( SEXP ir, SEXP inw, SEXP n, SEXP lag1, SEXP ip0, SEXP cov, SEXP lmax, SEXP mj0, SEXP mj1 )
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13,*i14,*i15,*i16,*i17,*i18;

    SEXP ans =  R_NilValue, l = R_NilValue, aic = R_NilValue, oaic = R_NilValue, mo = R_NilValue, v = R_NilValue, ac = R_NilValue;
    SEXP nc = R_NilValue, m1 = R_NilValue, m2 = R_NilValue, w = R_NilValue, z = R_NilValue, Rs = R_NilValue, chi = R_NilValue, ndt = R_NilValue, dic = R_NilValue, dicm = R_NilValue, po = R_NilValue;
    SEXP f = R_NilValue, k = R_NilValue, nh = R_NilValue, g = R_NilValue, ivf = R_NilValue, vf = R_NilValue;

    double *xaic, *xoaic, *xv, *xac, *xw, *xz, *xRs, *xchi, *xdic, *xdicm, *xf, *xg, *xvf = NULL;
    int      *xl, *xmo, *xnc, *xm1, *xm2, *xndt, *xpo, *xk, *xnh, *xivf = NULL;
    int      i, d, nj0, nj1, nj12, nj13;

    i1 = INTEGER_POINTER(ir);
    i2 = INTEGER_POINTER(inw);
    i3 = INTEGER_POINTER(n);
    i4 = INTEGER_POINTER(lag1);
    i5 = INTEGER_POINTER(ip0);
    d1 = NUMERIC_POINTER(cov);
    i16 = INTEGER_POINTER(lmax);
    i17 = INTEGER_POINTER(mj0);
    i18 = INTEGER_POINTER(mj1);

    d = *i1;
    nj0 = *i17;
    nj1 = *i18;
    nj12 = nj1 * nj1;
    nj13 = nj1 * nj12;
    PROTECT(ans = allocVector(VECSXP, 23));
    SET_VECTOR_ELT(ans, 0, l = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 1, aic = allocVector(REALSXP, nj0));
    SET_VECTOR_ELT(ans, 2, oaic = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 3, mo = allocVector(INTSXP, 1)); 
    SET_VECTOR_ELT(ans, 4, v = allocVector(REALSXP, d*d)); 
    SET_VECTOR_ELT(ans, 5, ac = allocVector(REALSXP, nj0*d*d));
    SET_VECTOR_ELT(ans, 6, nc = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 7, m1 = allocVector(INTSXP, nj1));
    SET_VECTOR_ELT(ans, 8, m2 = allocVector(INTSXP, nj1));
    SET_VECTOR_ELT(ans, 9, w = allocVector(REALSXP, nj13));
    SET_VECTOR_ELT(ans, 10, z = allocVector(REALSXP, nj12)); 
    SET_VECTOR_ELT(ans, 11, Rs = allocVector(REALSXP, nj12));
    SET_VECTOR_ELT(ans, 12, chi = allocVector(REALSXP, nj12)); 
    SET_VECTOR_ELT(ans, 13, ndt = allocVector(INTSXP, nj12));
    SET_VECTOR_ELT(ans, 14, dic = allocVector(REALSXP, nj12)); 
    SET_VECTOR_ELT(ans, 15, dicm = allocVector(REALSXP, nj1)); 
    SET_VECTOR_ELT(ans, 16, po = allocVector(INTSXP, nj1)); 
    SET_VECTOR_ELT(ans, 17, f = allocVector(REALSXP, nj12)); 
    SET_VECTOR_ELT(ans, 18, k = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 19, nh = allocVector(INTSXP, nj1));
    SET_VECTOR_ELT(ans, 20, g = allocVector(REALSXP, nj1*d)); 
    SET_VECTOR_ELT(ans, 21, ivf = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 22, vf = allocVector(REALSXP, nj12)); 

    i6 = INTEGER_POINTER(l);
    d2 = NUMERIC_POINTER(aic);
    d3 = NUMERIC_POINTER(oaic);
    i7 = INTEGER_POINTER(mo);
    d4 = NUMERIC_POINTER(v);
    d5 = NUMERIC_POINTER(ac);
    i8 = INTEGER_POINTER(nc);
    i9 = INTEGER_POINTER(m1);
    i10 = INTEGER_POINTER(m2);
    d6 = NUMERIC_POINTER(w);
    d7 = NUMERIC_POINTER(z);
    d8 = NUMERIC_POINTER(Rs);
    d9 = NUMERIC_POINTER(chi);
    i11 = INTEGER_POINTER(ndt);
    d10 = NUMERIC_POINTER(dic);
    d11 = NUMERIC_POINTER(dicm);
    i12 = INTEGER_POINTER(po);
    d12 = NUMERIC_POINTER(f);
    i13 = INTEGER_POINTER(k);
    i14 = INTEGER_POINTER(nh);
    d13 = NUMERIC_POINTER(g);
    i15 = INTEGER_POINTER(ivf);
    d14 = NUMERIC_POINTER(vf);

    F77_CALL(canocaf) (i1,i2,i3,i4,i5,d1,i6,d2,d3,i7,d4,d5,i8,i9,i10,d6,d7,d8,d9,i11,d10,d11,i12,d12,i13,i14,d13,i15,d14,i16,i17,i18);

    xl = INTEGER(l);
    xaic = REAL(aic);
    xoaic = REAL(oaic);
    xmo = INTEGER(mo);
    xv = REAL(v);
    xac = REAL(ac);
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
    xf = REAL(f);
    xk = INTEGER(k);
    xnh = INTEGER(nh);
    xg = REAL(g);
    xivf = INTEGER(ivf);
    xvf = REAL(vf);

    *xl = *i6;
    for(i=0; i<nj0; i++) xaic[i] = d2[i];
    *xoaic = *d3;
    *xmo = *i7;
    for(i=0; i<d*d; i++) xv[i] = d4[i];
    for(i=0; i<nj0*d*d; i++) xac[i] = d5[i];
    *xnc = *i8;
    for(i=0; i<nj1; i++) xm1[i] = i9[i];
    for(i=0; i<nj1; i++) xm2[i] = i10[i];
    for(i=0; i<nj13; i++) xw[i] = d6[i];
    for(i=0; i<nj12; i++) xz[i] = d7[i];
    for(i=0; i<nj12; i++) xRs[i] = d8[i];
    for(i=0; i<nj12; i++) xchi[i] = d9[i];
    for(i=0; i<nj12; i++) xndt[i] = i11[i];
    for(i=0; i<nj12; i++) xdic[i] = d10[i];
    for(i=0; i<nj1; i++) xdicm[i] = d11[i];
    for(i=0; i<nj1; i++) xpo[i] = i12[i];
    for(i=0; i<nj12; i++) xf[i] = d12[i];
    *xk = *i13;
    for(i=0; i<nj1; i++) xnh[i] = i14[i];
    for(i=0; i<nj1*d; i++) xg[i] = d13[i];
    *xivf = *i15;
    for(i=0; i<nj12; i++) xvf[i] = d14[i];

    UNPROTECT(1);
    return ans; 
}


