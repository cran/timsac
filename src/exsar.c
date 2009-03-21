#include "timsac.h"
#include <R_ext/RS.h>

extern int exsar(double*,long*,long*,double*,double*,double*,double*,double*,long*,double*,double*,double*,double*,double*,char**);

/* rtimsac78.dll subroutine */ int exsar(d1,i1,i2,d2,d3,d4,d5,d6,i3,d7,d8,d9,d10,d11,c1)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11;
	long *i1,*i2,*i3;
	char **c1;

{
	extern int F77_NAME(exsarf)(double*,long*,long*,double*,double*,double*,double*,double*,long*,double*,double*,double*,double*,double*,long*);

	long *i4;
	i4 = (long*)*c1;

	F77_CALL(exsarf) (d1,i1,i2,d2,d3,d4,d5,d6,i3,d7,d8,d9,d10,d11,i4);

	return 0;
}
