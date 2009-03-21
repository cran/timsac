#include "timsac.h"
#include <R_ext/RS.h>

extern int bispec(long*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

/* rtimsac74.dll subroutine */ int bispec(i1,i2,d1,d2,d3,d4,d5,d6,d7,d8,d9)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
	long *i1,*i2;

{
	extern int F77_NAME(bispecf)(long*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

	F77_CALL(bispecf) (i1,i2,d1,d2,d3,d4,d5,d6,d7,d8,d9);

	return 0;
}
