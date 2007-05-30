#include "timsac.h"
extern int fftcor(long*, long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */ int fftcor(i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
	long *i1,*i2,*i3,*i4,*i5;

{
	extern int fftcorf_(long*, long*, long*, long*, long*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

	fftcorf_(i1,i2,i3,i4,i5,d1,d2,d3,d4,d5,d6,d7,d8,d9);

	return 0;
}
