/*
*  Timsac : Time Series Analysis and ControlProgram Package
*  Copyright (C) 2005    The Institute of Statistical Mathematics
*
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program; if not, write to the Free Software
*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
*
*
*  ismrp at grp.ism.ac.jp
*/

#include <R.h>
#include <Rinternals.h>
#include <libintl.h>

#define _(String) (String)

/* Fortran : */

/*** timsac72 ***/

void F77_NAME(auspecf)(int *n, int *lag1, double *acov, double *spec1,
              double *spec2, double *stat);

void F77_NAME(autcorf)(double *y, int *n, double *acov, double *acor, int *lag1,
              double *mean);

void F77_NAME(fftcorf)(int *ld, int *lag1, int *n, int *n2p, int *isw,
              double *x1, double *y1, double *acov, double *cov21,
              double *cov12, double *acor, double *cor21, double *cor12,
              double *mean);

void F77_NAME(fpeautf)(int *morder, int *n, double *sd, double *cxx,
              double *sig2, double *fpe, double *rfpe, double *parcor,
              double *chi2, double *ofpe, double *fpem, double *rfpem, int *mo,
			  double *sig2m,  double *a,  double *ao);

void F77_NAME(fpec7)(int *n, int *morder, int *ncon, int *ip, int *d, int *inw,
              double *cov, double *ccv, double *fpec, double *rfpec,
              double *aic, int *mo, double *fpecm, double *rfpecm, double *aicm,
              double *perr, double *arcoef);

void F77_NAME(mulcorf)(double *y, int *n, int *d, int *lag1, double *mean,
              double *cov, double *cor);

void F77_NAME(mulfrff)(int *nv, int *iovar, int *n, int *lag1, int *d,
              double *spec, Rcomplex *cospec, double *freqr, double *freqi,
              double *gain, double *phase, double *pcoh, double *errstat,
              double *mcoh);

void F77_NAME(mulnosf)(int *h, int *l, int *ip, double *sd, double *arcoef,
              double *nsd, double *drpc, double *irpc);

void F77_NAME(mulrspf)(int *h, int *l, int *d, int *k, double *cov,
              double *arcoef, double *macoef, Rcomplex *rspec, double *scoh);

void F77_NAME(mulspef)(int *n, int *d, int *lag1, int *lag3, double *cov,
              double *spec1, double *spec2, double *stat, double *coh1,
              double *coh2);

void F77_NAME(optdesf)(int *ncon, int *nman, int *ns, int *order, double *q,
              double *r, double *osd, double *a, double *b, double *gain);

void F77_NAME(optsimf)(int *ns, int *order, int *ir, int *il, double *trans,
              double *gamma, double *gain, double *noise, double *x, double *y,
              double *xmean, double *ymean, double *x2sum, double *y2sum,
              double *x2mean, double *y2mean, double *xvar, double *yvar);

void F77_NAME(raspecf)(int *h, int *l, int *k, double *var, double *arcoef,
              double *macoef, double *rspec);

void F77_NAME(sglfref)(int *invar, int *outvar, int *n, int *lag1, int *d,
              double *spec, double *spec1, double *spec2, double *cspec,
              double *qspec, double *gain, double *coh, double *freqr,
              double *freqi, double *err, double *phase);
	
void F77_NAME(wnoisef)(int *len, int *ir, double *perr, double *wn);


/*** timsac74 ***/

void F77_NAME(autarm)(int *n, int *morder1, double *autcv, int *inc, int *p1,
              double *arcoef1, int *q1, double *macoef1, int *newn, int *p,
              double *a,  int *q, double *b, double *std, double *v, double *gr,
              double *aic, double *aicm, int *pbest, int *qbest, int *lmax,
              int *mmax, int *nmax);
			  
void F77_NAME(bispecf)(int *n, int *lag, double *cv, double *tmnt,
              double *pspec1, double *pspec2, double *sig, double *ch,
              double *br, double *bi, double *rat);

void F77_NAME(canarmf)(int *n, int *morder1, double *autcv, double *arcoef,
              int *l1, double *v, double *aic, double *oaic, int *mo,
              double *parcor, int *nc, int *m1, int *m2, double *w, double *z,
              double *Rs, double *chi, int *ndt, double *dic, double *dicm,
              int *po, int *k, double *b, int *l, double *a, int *mmax,
              int *nmax);

void F77_NAME(canocaf)(int *ir, int *inw, int *n, int *lag1, int *ip0,
              double *cov, int *l, double *aic, double *oaic, int *mo,
              double *v, double *ac, int *nc, int *m1, int *m2, double *w,
              double *z, double *Rs, double *chi, int *ndt, double *dic,
              double *dicm, int *po, double *f, int *k, int *nh, double *g,
              int *ivf, double *vf, int *lmax, int *mj0, int *mj1);

void F77_NAME(covgenf)(int *lag, int *k, double *f, double *gain, double *acov,
              double *acor);

void F77_NAME(markovf)(int *n, int *lag1, int *d, double *cov, int *k, int *nh,
              int *nvf, double *vectF, double *matGi, int *icont, int *id,
              int *ir, int *ij, int *ik, int *ngr, double *gr, double *a1,
              double *a, double *b, double *vd, int *iqm, double *bm,
              double *au, double *zz, double *v, double *aic, int *mj3,
              int *mj4, int *mj6, int *mj7);

void F77_NAME(nonstf)(int *n, int *span, double *y, int *ns, int *morder,
              int *p, double *coef, double *v, double *aic, double *daic21,
			  double *daic, int *ks, int *ke, double *pspec);

void F77_NAME(prdctrf)(int *n, int *r, int *s, int *h, int *d, int *p, int *q,
              int *jsw, double *y, double *arcoef, double *macoef,
              double *impuls, double *v, double *yreal, double *yori,
              double *ypre, double *ys, double *z1, double *z2, double *z3,
              double *zz1, double *zz2, double *zz3);

void F77_NAME(simconf)(int *d, int *k, int *span, int *len, int *r,
              double *arcoef, double *impuls, double *v, double *weight,
			  double *bc, double *bd, double *g, double *av, double *si,
              double *s2);

void F77_NAME(thirmof)(int *n, int *lag, double *y, double *mean, double *acov,
              double *acor, double *mnt);


/*** timsac78 ***/

void F77_NAME(blocarf)(double *y, int *n, int *sorder, int *span, int *ns,
              double *mean, double *var, double *aic, double *bw, double *b,
              double *a, double *v, int *ks, int *ke, double *pxx);

void F77_NAME(blomarf)(double *y, int *n, int *d, double *calb, int *morder,
              int *span, int *ns, double *mean, double *var, double *bw,
              double *raic, double *a, double *e, double *aic, int *ks,
              int *ke, int *nns);

void F77_NAME(bsubstf)(double *y, int *n, int *mtype, int *lag, int *nreg,
              int *cstep, int *reg, int *tlag, double *ymean, double *yvar,
              int *m, double *aicm, double *vm, double *a1, double *v,
              double *aic, double *daic, double *aicb, double *vb, double *ek,
              double *a2, int *ind, double *c, double *c1, double *c2,
              double *b, double *eicmin, double *esum, double *npm,
              double *npmnreg, double *e, double *mean, double *var,
              double *skew, double *peak, double *cov, double *pxx);

void F77_NAME(exsarf)(double *y, int *n, int *morder, double *mean, double *var,
              double *v, double *aic, double *daic, int *m, double *aicm,
              double *sdm1, double *a1, double *sdm2, double *a2, int *ier);

void F77_NAME(mlocarf)(double *y, int *n, int *morder, int *span, int *cnst,
              int *ns, double *mean, double *var, double *a, int *mf,
              double *sdf, int *ks, int *ke, double *pxx, int *ld1, int *ld2,
              int *ms, double *sdms, double *aics, int *mp, double *sdmp,
              double *aicp);

void F77_NAME(mlomarf)(double *y, int *n, int *d, double *calb, int *morder,
              int *span, int *cnst, int *ns, double *mean, double *var,
              int *ld1, int *ld2, int *ms, double *aicm, int *mp, double *aicc,
              int *mf, double *aic, double *a, double *e, int *ks, int *ke,
              int *nns);

void F77_NAME(mulbarf)(double *y, int *n, int *d, double *calb, int *morder,
              double *mean, double *var, double *v, double *aic, double *daic,
              int *m, double *aicm, double *vm, double *w1, double *w2,
              double *a, double *b, double *g, double *h, double *e,
              double *aicb);

void F77_NAME(mulmarf)(double *y, int *n, int *d, double *calb, int *lag,
              double *mean, double *var, double *v, double *aic, double *daic,
              int *m, double *aicm, double *vm, int *npr, int *jnd, double *a,
              double *rv, double *aicf, double *ei, double *bi, double *matv,
              double *arcoef, int *morder, double *aics);

void F77_NAME(perarsf)(double *y, int *n, int *ni, int *lag, int *ksw,
              double *mean, double *var, int *np, int *jnd, double *a,
              double *aic, double *b, double *v, double *c, double *osd,
              int *mo);

void F77_NAME(unibarf)(double *y, int *n, int *arorder, double *mean,
              double *var, double *v, double *aic, double *daic, int *m,
              double *aicm, double *vm, double *pa, double *bw, double *sbw,
              double *pab, double *aicb, double *vb, double *pn, double *a,
              double *pxx);

void F77_NAME(unimarf)(double *y, int *n, int *morder, double *mean,
              double *var, double *v, double *aic, double *daic, int *m,
              double *aicm, double *vm, double *a);

void F77_NAME(xsarmaf)(double *y, int *n, int *p, int *q, double *p01,
              double *g1, double *tl1, double *p02, double *g2, double *alpar,
              double *alpma, double *tl2, double *sigma2);


/*** timsac84 ***/

void F77_NAME(bayseaf)(double *y, int *ndata, int *forecast, double *outlier,
              double *dmoi, double *trend, double *season, double *tdcmp,
              double *irreg, double *adjust, double *est, double *psds,
              double *psdt, double *avabic, int *ipara, double *para,
              double *arft, double *arfs, double *arfn, int *iart, int *iars,
              int *iarn);

void F77_NAME(decompf)(double *y, int *n, int *ipar, double *trend,
              double *seasonal, double *ar, double *trad, double *noise,
              double *para, int *miss, double *omax, int *ier);

void F77_NAME(spgrh)(double *y, int *n, int *lag1, int *ifpl1, int *mode,
              int *period, double *cxx, double *cn, double *xmean, double *sd,
              double *aic, double *parcor, double *pxx, int *ier);
