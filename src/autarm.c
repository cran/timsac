#include <R.h>
#include <Rdefines.h>
#include "timsac.h"

extern void F77_NAME(autarmf) (int*,int*,double*,int*,int*,double*,int*,double*,int*,int*,double*,int*,double*,double*,double*,double*,double*,double*,int*,int*,int*,int*,int*);

SEXP AutarmC( SEXP n, SEXP morder1, SEXP autcv, SEXP inc, SEXP p1, SEXP arcoef1, SEXP q1, SEXP macoef1, SEXP lmax, SEXP mmax, SEXP nmax )
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
    int      *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13;

    SEXP ans =  R_NilValue, newn = R_NilValue, p = R_NilValue, a = R_NilValue, q = R_NilValue, b = R_NilValue;
    SEXP std =  R_NilValue, v = R_NilValue, gr = R_NilValue, aic = R_NilValue, aicm = R_NilValue, pbest = R_NilValue, qbest = R_NilValue;
    double *xa, *xb, *xstd, *xv, *xgr, *xaic, *xaicm = NULL;
    int      *xnewn, *xp, *xq, *xpbest, *xqbest = NULL;
    int      i, mm, nn, mn;

    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(morder1);
    d1 = NUMERIC_POINTER(autcv);
    i3 = INTEGER_POINTER(inc);
    i4 = INTEGER_POINTER(p1);
    d2 = NUMERIC_POINTER(arcoef1);
    i5 = INTEGER_POINTER(q1);
    d3 = NUMERIC_POINTER(macoef1);
    i11 = INTEGER_POINTER(lmax);
    i12 = INTEGER_POINTER(mmax);
    i13 = INTEGER_POINTER(nmax);

    mm = *i12;
    nn = *i13;
    mn = mm * nn;
    PROTECT(ans = allocVector(VECSXP, 12));
    SET_VECTOR_ELT(ans, 0, newn = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 1, p = allocVector(INTSXP, nn)); 
    SET_VECTOR_ELT(ans, 2, a = allocVector(REALSXP, mn));
    SET_VECTOR_ELT(ans, 3, q = allocVector(INTSXP, nn));
    SET_VECTOR_ELT(ans, 4, b = allocVector(REALSXP, mn));
    SET_VECTOR_ELT(ans, 5, std = allocVector(REALSXP, mn)); 
    SET_VECTOR_ELT(ans, 6, v = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 7, gr = allocVector(REALSXP, mn));
    SET_VECTOR_ELT(ans, 8, aic = allocVector(REALSXP, nn)); 
    SET_VECTOR_ELT(ans, 9, aicm = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 10, pbest = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 11, qbest = allocVector(INTSXP, 1)); 

    i6 = INTEGER_POINTER(newn);
    i7 = INTEGER_POINTER(p);
    d4 = NUMERIC_POINTER(a);
    i8 = INTEGER_POINTER(q);
    d5 = NUMERIC_POINTER(b);
    d6 = NUMERIC_POINTER(std);
    d7 = NUMERIC_POINTER(v);
    d8 = NUMERIC_POINTER(gr);
    d9 = NUMERIC_POINTER(aic);
    d10 = NUMERIC_POINTER(aicm);
    i9 = INTEGER_POINTER(pbest);
    i10 = INTEGER_POINTER(qbest);

    F77_CALL(autarmf)(i1,i2,d1,i3,i4,d2,i5,d3,i6,i7,d4,i8,d5,d6,d7,d8,d9,d10,i9,i10,i11,i12,i13);

    xnewn = INTEGER(newn);
    xp = INTEGER(p);
    xa = REAL(a);
    xq = INTEGER(q);
    xb = REAL(b);
    xstd = REAL(std);
    xv = REAL(v);
    xgr = REAL(gr);
    xaic = REAL(aic);
    xaicm = REAL(aicm);
    xpbest = INTEGER(pbest);
    xqbest = INTEGER(qbest);

    *xnewn = *i6;
    for(i=0; i<nn; i++) xp[i] = i7[i];
    for(i=0; i<mn; i++) xa[i] = d4[i];
    for(i=0; i<nn; i++) xq[i] = i8[i];
    for(i=0; i<mn; i++) xb[i] = d5[i];
    for(i=0; i<mn; i++) xstd[i] = d6[i];
    for(i=0; i<nn; i++) xv[i] = d7[i];
    for(i=0; i<mn; i++) xgr[i] = d8[i];
    for(i=0; i<nn; i++) xaic[i] = d9[i];
    for(i=0; i<nn; i++) xaicm[i] = d10[i];
    *xpbest = *i9;
    *xqbest = *i10;

    UNPROTECT(1);
    return ans; 
}
