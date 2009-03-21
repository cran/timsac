#include "timsac.h"
#include <R_ext/RS.h>

extern int optsim(long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */	int optsim(i1,i2,i3,i4,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14;
	long *i1,*i2,*i3,*i4;

{
	extern int F77_NAME(optsimf) (long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

	F77_CALL(optsimf) (i1,i2,i3,i4,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14);

	return 0;
}
