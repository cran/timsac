#include "timsac.h"
extern int simcon(long*,long*,long*,long*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

/* rtimsac74.dll subroutine */ int simcon(i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
	long *i1,*i2,*i3,*i4,*i5;

{
	extern int simconf_(long*,long*,long*,long*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

	simconf_(i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10);

	return 0;
}
