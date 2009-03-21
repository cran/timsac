#include "timsac.h"
#include <R_ext/RS.h>

extern int decomp(double*, long*, long*, double*, double*, double*, double*, double*, double*, long*, double*);

/* rdecomp subroutine */ int decomp(d1,i1,i2,d2,d3,d4,d5,d6,d7,i3,d8)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8;
	long *i1,*i2,*i3;

{
	extern int F77_NAME(decompf) (double*, long*, long*, double*, double*, double*, double*, double*, double*, long*, double*);
  
	F77_CALL(decompf) (d1,i1,i2,d2,d3,d4,d5,d6,d7,i3,d8);

	return 0;
}
