#include "timsac.h"
#include <R_ext/RS.h>

extern int xsarma(double*,long*,long*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

/* rtimsac78.dll subroutine */ int xsarma(d1,i1,i2,i3,d2,d3,d4,d5,d6,d7,d8,d9,d10)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
	long *i1,*i2,*i3;

{
	extern int F77_NAME(xsarmaf) (double*,long*,long*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

	F77_CALL(xsarmaf) (d1,i1,i2,i3,d2,d3,d4,d5,d6,d7,d8,d9,d10);

	return 0;
}
