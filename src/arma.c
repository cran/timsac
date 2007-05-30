#include "timsac.h"
extern int arma(long*, long*, double*, double*, double*, long*, long*, long*, double*, double*, double*, double*, double*, double*, long*, long*);

int arma(i1,i2,d1,d2,d3,i3,i4,i5,d4,d5,d6,d7,d8,d9,i6,i7)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7;

{
	extern int armaf_(long*, long*, double*, double*, double*, long*, long*, long*, double*, double*, double*, double*, double*, double*, long*, long*);

	armaf_(i1,i2,d1,d2,d3,i3,i4,i5,d4,d5,d6,d7,d8,d9,i6,i7);

	return 0;
}

