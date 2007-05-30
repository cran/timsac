#include "timsac.h"
extern int tvar(double*, long*, long*, long*, long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, double*);

int tvar(d1,i1,i2,i3,i4,i5,i6,i7,d2,d3,d4,d5,d6,d7,d8,d9)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7;

{
	extern int tvarf_(double*, long*, long*, long*, long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, double*);

	tvarf_(d1,i1,i2,i3,i4,i5,i6,i7,d2,d3,d4,d5,d6,d7,d8,d9);

	return 0;
}

