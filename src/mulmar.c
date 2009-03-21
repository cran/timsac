#include "timsac.h"
#include <R_ext/RS.h>

extern int mulmar(double*,long*,long*,double*,long*,double*,double*,double*,double*,double*,long*,double*,double*,long*,long*,double*,double*,double*,double*,double*,double*,double*,long*,double*,char**);

/* rtimsac78.dll subroutine */ int mulmar(d1,i1,i2,d2,i3,d3,d4,d5,d6,d7,i4,d8,d9,i5,i6,d10,d11,d12,d13,d14,d15,d16,i7,d17,c1)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7;
	char **c1;

{
	extern int F77_NAME(mulmarf)  (double*,long*,long*,double*,long*,double*,double*,double*,double*,double*,long*,double*,double*,long*,long*,double*,double*,double*,double*,double*,double*,double*,long*,double*,long*);

	long *i8;
	i8 = (long*)*c1;

	F77_CALL(mulmarf)  (d1,i1,i2,d2,i3,d3,d4,d5,d6,d7,i4,d8,d9,i5,i6,d10,d11,d12,d13,d14,d15,d16,i7,d17,i8);

	return 0;
}
