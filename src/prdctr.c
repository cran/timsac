#include "timsac.h"
extern int prdctr(long*,long*,long*,long*,long*,long*,long*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,char**);

/* rtimsac74.dll subroutine */ int prdctr(i1,i2,i3,i4,i5,i6,i7,i8,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,c1)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;
	char **c1;

{
	extern int prdctrf_(long*,long*,long*,long*,long*,long*,long*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,long*);

	long *i9;
	i9 = (long*)*c1;

	prdctrf_(i1,i2,i3,i4,i5,i6,i7,i8,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,i9);

	return 0;
}
