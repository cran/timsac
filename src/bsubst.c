#include "timsac.h"
#include <R_ext/RS.h>

extern int bsubst(double*,long*,long*,long*,long*,long*,long*,long*,char**,double*,double*,double*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

/* rtimsac78.dll subroutine */ int bsubst(d1,i1,i2,i3,i4,i5,i6,i7,c1,d2,d3,d4,i9,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,i10,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26,d27,d28,d29)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17,*d18,*d19,*d20,*d21,*d22,*d23,*d24,*d25,*d26,*d27,*d28,*d29;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i9,*i10;
	char **c1;
{

	extern int F77_NAME(bsubstf) (double*,long*,long*,long*,long*,long*,long*,long*,long*,double*,double*,double*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

	long *i8;
	i8 = (long*)*c1;

	F77_CALL(bsubstf) (d1,i1,i2,i3,i4,i5,i6,i7,i8,d2,d3,d4,i9,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,i10,d15,d16,d17,d18,d19,d20,d21,d22,d23,d24,d25,d26,d27,d28,d29);

	return 0;

}
