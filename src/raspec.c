#include "timsac.h"
#include <R_ext/RS.h>

extern int raspec(long*, long*, long*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */ int raspec(i1,i2,i3,d1,d2,d3,d4)

	double *d1,*d2,*d3,*d4;
	long *i1,*i2,*i3;

{
	extern int F77_NAME(raspecf) (long*, long*, long*, double*, double*, double*, double*);

	F77_CALL(raspecf) (i1,i2,i3,d1,d2,d3,d4);

	return 0;
}
