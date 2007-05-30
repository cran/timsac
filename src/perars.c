#include "timsac.h"
extern int perars(double*,long*,long*,long*,long*,double*,double*,long*,long*,double*,double*,double*,double*,double*,double*,long*);

/* rtimsac78.dll subroutine */ int perars(d1,i1,i2,i3,i4,d2,d3,i5,i6,d4,d5,d6,d7,d8,d9,i7)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7;

{
	extern int perarsf_(double*,long*,long*,long*,long*,double*,double*,long*,long*,double*,double*,double*,double*,double*,double*,long*);

	perarsf_(d1,i1,i2,i3,i4,d2,d3,i5,i6,d4,d5,d6,d7,d8,d9,i7);

	return 0;
}
