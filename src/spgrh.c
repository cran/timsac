#include "timsac.h"
#include <R_ext/RS.h>

extern int spgrh(double*, long*, long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, long*);

int spgrh(d1,i1,i2,i3,i4,i5,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,i6)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11;
	long *i1,*i2,*i3,*i4,*i5,*i6;

{
	extern int F77_NAME(spgrhf) (double*, long*, long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, long*);

	F77_CALL(spgrhf) (d1,i1,i2,i3,i4,i5,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,i6);

	return 0;
}

