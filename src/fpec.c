#include "timsac.h"
#include <R_ext/RS.h>

extern int fpec7(long*, long*, long*, long*, long*, long*, double*, double*, double*, double*, double*, long*, double*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */ int fpec7(i1,i2,i3,i4,i5,i6,d1,d2,d3,d4,d5,i7,d6,d7,d8,d9,d10)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7;

{
	extern int F77_NAME(fpec7f) (long*, long*, long*, long*, long*, long*, double*, double*, double*, double*, double*, long*, double*, double*, double*, double*, double*);

	F77_CALL(fpec7f) (i1,i2,i3,i4,i5,i6,d1,d2,d3,d4,d5,i7,d6,d7,d8,d9,d10);

    return 0;
}
