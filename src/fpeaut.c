#include "timsac.h"
extern int fpeaut(long*, long*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, long*, double*, double*, double*);

/* rtimsac72.dll subroutine */ int fpeaut(i1,i2,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,i3,d11,d12,d13)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13;
	long *i1,*i2,*i3;

{
	extern int fpeautf_(long*, long*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, long*, double*, double*, double*);

	fpeautf_(i1,i2,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,i3,d11,d12,d13);

    return 0;
}
