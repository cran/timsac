#include "timsac.h"
#include <R_ext/RS.h>

extern int thirmo(long*,long*,double*,double*,double*,double*,double*);

/* rtimsac74.dll subroutine */ int thirmo(i1,i2,d1,d2,d3,d4,d5)

	double *d1,*d2,*d3,*d4,*d5;
	long *i1,*i2;

{
	extern int F77_NAME(thirmof) (long*,long*,double*,double*,double*,double*,double*);

	F77_CALL(thirmof) (i1,i2,d1,d2,d3,d4,d5);

	return 0;
}
