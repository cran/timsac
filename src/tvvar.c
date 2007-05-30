#include "timsac.h"
extern int tvvar(double*, long*, long*, double*, long*, double*, double*, double*, double*, long*, double*, double*, double*, double*, double*, double*);

int tvvar(d1,i1,i2,d2,i3,d3,d4,d5,d6,i4,d7,d8,d9,d10,d11,d12)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12;
	long *i1,*i2,*i3,*i4;

{
	extern int tvvarf_(double*, long*, long*, double*, long*, double*, double*, double*, double*, long*, double*, double*, double*, double*, double*, double*);

	tvvarf_(d1,i1,i2,d2,i3,d3,d4,d5,d6,i4,d7,d8,d9,d10,d11,d12);

	return 0;
}
