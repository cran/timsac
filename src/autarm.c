#include "timsac.h"
extern int autarm(long*,long*,double*,long*,long*,double*,long*,double*,long*,long*,double*,long*,double*,double*,double*,double*,double*,double*,long*,long*,char**,long*,long*,long*);

/* rtimsac74.dll subroutine */ int autarm(i1,i2,d1,i3,i4,d2,i5,d3,i6,i7,d4,i8,d5,d6,d7,d8,d9,d10,i9,i10,c1,i11,i12,i13)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13;
	char **c1;

{
	extern int autarmf_(long*,long*,double*,long*,long*,double*,long*,double*,long*,long*,double*,long*,double*,double*,double*,double*,double*,double*,long*,long*,long*,long*,long*,long*);

	long *ic;
	ic = (long*)*c1;

	autarmf_(i1,i2,d1,i3,i4,d2,i5,d3,i6,i7,d4,i8,d5,d6,d7,d8,d9,d10,i9,i10,ic,i11,i12,i13);

	return 0;
}

