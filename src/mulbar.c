#include "timsac.h"
#include <R_ext/RS.h>

extern int mulbar(double*,long*,long*,double*,long*,double*,double*,double*,double*,double*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

/* rtimsac78.dll subroutine */ int mulbar(d1,i1,i2,d2,i3,d3,d4,d5,d6,d7,i4,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17;
	long *i1,*i2,*i3,*i4;

{
	extern int F77_NAME(mulbarf) (double*,long*,long*,double*,long*,double*,double*,double*,double*,double*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

	F77_CALL(mulbarf) (d1,i1,i2,d2,i3,d3,d4,d5,d6,d7,i4,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17);

	return 0;
}
