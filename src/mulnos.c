#include "timsac.h"
#include <R_ext/RS.h>

extern int mulnos(long*, long*, long*, double*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */ int mulnos(i1,i2,i3,d1,d2,d3,d4,d5)

	double *d1,*d2,*d3,*d4,*d5;
	long *i1,*i2,*i3;

{
	extern int F77_NAME(mulnosf) (long*, long*, long*, double*, double*, double*, double*, double*);

	F77_CALL(mulnosf) (i1,i2,i3,d1,d2,d3,d4,d5);

	return 0;
}
