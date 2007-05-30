#include "timsac.h"
extern int optdes(long*, long*, long*, long*, double*, double*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */	int optdes(i1,i2,i3,i4,d1,d2,d3,d4,d5,d6)

	double *d1,*d2,*d3,*d4,*d5,*d6;
	long *i1,*i2,*i3,*i4;

{
	extern int optdesf_(long*, long*, long*, long*, double*, double*, double*, double*, double*, double*);

	optdesf_(i1,i2,i3,i4,d1,d2,d3,d4,d5,d6);

	return 0;
}
