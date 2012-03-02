#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(prdctrf) (int*,int*,int*,int*,int*,int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

SEXP prdctr( SEXP n, SEXP r, SEXP s, SEXP h, SEXP d, SEXP p, SEXP q, SEXP jsw, SEXP y, SEXP arcoef, SEXP macoef, SEXP impuls, SEXP v )
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15;
    int    *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;

    SEXP ans =  R_NilValue, yreal = R_NilValue, yori = R_NilValue, ypre = R_NilValue, ys = R_NilValue, z1 = R_NilValue, z2 = R_NilValue, z3 = R_NilValue, zz1 = R_NilValue, zz2 = R_NilValue, zz3 = R_NilValue;

    double *xyreal, *xyori, *xypre, *xys, *xz1, *xz2, *xz3, *xzz1, *xzz2, *xzz3 = NULL;
    int      i, nn, ns, nh, nd, k;

    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(r);
    i3 = INTEGER_POINTER(s);
    i4 = INTEGER_POINTER(h);
    i5 = INTEGER_POINTER(d);
    i6 = INTEGER_POINTER(p);
    i7 = INTEGER_POINTER(q);
    i8 = INTEGER_POINTER(jsw);
    d1 = NUMERIC_POINTER(y);
    d2 = NUMERIC_POINTER(arcoef);
    d3 = NUMERIC_POINTER(macoef);
    d4 = NUMERIC_POINTER(impuls);
    d5 = NUMERIC_POINTER(v);

    nn = *i1;
    ns = *i3;
    nh = *i4;
    nd = *i5;
    k = (ns+nh)*nd;
    PROTECT(ans = allocVector(VECSXP, 10));
    SET_VECTOR_ELT(ans, 0, yreal = allocVector(REALSXP, k));
    SET_VECTOR_ELT(ans, 1, yori = allocVector(REALSXP, (nh+1)*nd));
    SET_VECTOR_ELT(ans, 2, ypre = allocVector(REALSXP, k)); 
    SET_VECTOR_ELT(ans, 3, ys = allocVector(REALSXP, nn*nd));
    SET_VECTOR_ELT(ans, 4, z1 = allocVector(REALSXP, k));
    SET_VECTOR_ELT(ans, 5, z2 = allocVector(REALSXP, k)); 
    SET_VECTOR_ELT(ans, 6, z3 = allocVector(REALSXP, k)); 
    SET_VECTOR_ELT(ans, 7, zz1 = allocVector(REALSXP, k));
    SET_VECTOR_ELT(ans, 8, zz2 = allocVector(REALSXP, k));
    SET_VECTOR_ELT(ans, 9, zz3 = allocVector(REALSXP, k)); 

    d6 = NUMERIC_POINTER(yreal);
    d7 = NUMERIC_POINTER(yori);
    d8 = NUMERIC_POINTER(ypre);
    d9 = NUMERIC_POINTER(ys);
    d10 = NUMERIC_POINTER(z1);
    d11 = NUMERIC_POINTER(z2);
    d12 = NUMERIC_POINTER(z3);
    d13 = NUMERIC_POINTER(zz1);
    d14 = NUMERIC_POINTER(zz2);
    d15 = NUMERIC_POINTER(zz3);

    F77_CALL(prdctrf) (i1,i2,i3,i4,i5,i6,i7,i8,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15);

    xyreal = REAL(yreal);
    xyori = REAL(yori);
    xypre = REAL(ypre);
    xys = REAL(ys);
    xz1 = REAL(z1);
    xz2 = REAL(z2);
    xz3 = REAL(z3);
    xzz1 = REAL(zz1);
    xzz2 = REAL(zz2);
    xzz3 = REAL(zz3);

    for(i=0; i<k; i++) xyreal[i] = d6[i];
    for(i=0; i<(nh+1)*nd; i++) xyori[i] = d7[i];
    for(i=0; i<k; i++) xypre[i] = d8[i];
    for(i=0; i<nn*nd; i++) xys[i] = d9[i];
    for(i=0; i<k; i++) xz1[i] = d10[i];
    for(i=0; i<k; i++) xz2[i] = d11[i];
    for(i=0; i<k; i++) xz3[i] = d12[i];
    for(i=0; i<k; i++) xzz1[i] = d13[i];
    for(i=0; i<k; i++) xzz2[i] = d14[i];
    for(i=0; i<k; i++) xzz3[i] = d15[i];

    UNPROTECT(1);
    return ans; 
}

