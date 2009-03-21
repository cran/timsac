#include "timsac.h"
#include <R_ext/RS.h>

extern int auspec(long*, long*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */ int auspec(i1,i2,d1,d2,d3,d4)

	double *d1,*d2,*d3,*d4;
	long *i1,*i2;

{
	extern int F77_NAME(auspecf) (long*, long*, double*, double*, double*, double*);

	F77_CALL(auspecf) (i1,i2,d1,d2,d3,d4);

	return 0;
}

