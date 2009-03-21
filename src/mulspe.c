#include "timsac.h"
#include <R_ext/RS.h>

extern int mulspe(long*, long*, long*, long*, double*, double*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */ int mulspe(i1,i2,i3,i4,d1,d2,d3,d4,d5,d6)

	double *d1,*d2,*d3,*d4,*d5,*d6;
	long *i1,*i2,*i3,*i4;

{
	extern int F77_NAME(mulspef) (long*, long*, long*, long*, double*, double*, double*, double*, double*, double*);

	F77_CALL(mulspef) (i1,i2,i3,i4,d1,d2,d3,d4,d5,d6);

	return 0;
}
