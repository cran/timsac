#include "timsac.h"
extern int nonst(long*,long*,double*,long*,long*,long*,double*,double*,double*,double*,double*,long*,long*,double*,char**);

/* rtimsac74.dll subroutine */ int nonst(i1,i2,d1,i3,i4,i5,d2,d3,d4,d5,d6,i6,i7,d7,c1)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7;
	char **c1;

{
	extern int nonstf_(long*,long*,double*,long*,long*,long*,double*,double*,double*,double*,double*,long*,long*,double*,long*);

	long *i8;
	i8 = (long*)*c1;

	nonstf_(i1,i2,d1,i3,i4,i5,d2,d3,d4,d5,d6,i6,i7,d7,i8);

	return 0;
}
