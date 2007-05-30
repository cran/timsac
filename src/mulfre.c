#include "timsac.h"

struct Rcomplex {
	double r;
	double i;
};

extern int mulfrf(long*, long*, long*, long*, long*, double*, struct Rcomplex*, double*, double*, double*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */ int mulfrf(i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9)

	double *d1,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
	long *i1,*i2,*i3,*i4,*i5;
	struct Rcomplex *d2;

{
	extern int mulfrff_(long*, long*, long*, long*, long*, double*, struct Rcomplex*, double*, double*, double*, double*, double*, double*, double*);

	mulfrff_(i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9);

	return 0;
}
