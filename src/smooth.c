#include "timsac.h"
extern int smooth(double*, long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, long*, long*, double*, double*, long*, long*, long*, double*, double*, double*, double*);

int smooth(d1,i1,i2,i3,i4,d2,d3,d4,d5,d6,d7,d8,i5,i6,d9,d10,i7,i8,i9,d11,d12,d13,d14)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9;

{
	extern int smoothf_(double*, long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, long*, long*, double*, double*, long*, long*, long*, double*, double*, double*, double*);

	smoothf_(d1,i1,i2,i3,i4,d2,d3,d4,d5,d6,d7,d8,i5,i6,d9,d10,i7,i8,i9,d11,d12,d13,d14);

	return 0;
}
