#include "timsac.h"
#include <R_ext/RS.h>

extern int wnoise(long*, long*, double*, double*);

/* rtimsac72.dll subroutine */	int wnoise(i1,i2,d1,d2)

	double *d1,*d2;
	long *i1,*i2;

{
	extern int F77_NAME(wnoisef) (long*, long*, double*, double*);

	F77_CALL(wnoisef) (i1,i2,d1,d2);

	return 0;
}
