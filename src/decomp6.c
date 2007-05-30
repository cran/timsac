#include "timsac.h"
extern int decomp(double*, long*, long*, double*, double*, double*, double*, double*, double*, long*, double*);

/* rdecomp subroutine */ int decomp(d1,i1,i2,d2,d3,d4,d5,d6,d7,i3,d8)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8;
	long *i1,*i2,*i3;

{
	extern int decompf_(double*, long*, long*, double*, double*, double*, double*, double*, double*, long*, double*);
  
	decompf_(d1,i1,i2,d2,d3,d4,d5,d6,d7,i3,d8);

	return 0;
}

