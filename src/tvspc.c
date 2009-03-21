#include "timsac.h"
#include <R_ext/RS.h>

extern int tvspc(long*, long*, long*, long*, long*, double*, double*, double*, double*);

int tvspc(i1,i2,i3,i4,i5,d1,d2,d3,d4)

	double *d1,*d2,*d3,*d4;
	long *i1,*i2,*i3,*i4,*i5;

{
	extern int F77_NAME(tvspcf) (long*, long*, long*, long*, long*, double*, double*, double*, double*);

	F77_CALL(tvspcf) (i1,i2,i3,i4,i5,d1,d2,d3,d4);

	return 0;
}
