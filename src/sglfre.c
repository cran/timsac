#include "timsac.h"
extern int sglfre(long*, long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */  int sglfre(i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11;
	long *i1,*i2,*i3,*i4,*i5;

{
	extern int sglfref_(long*, long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

	sglfref_(i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11);

	return 0;
}
