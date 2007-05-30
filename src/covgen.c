#include "timsac.h"
extern int covgen(long*, long*, double*, double*, double*, double*);

/* rtimsac74.dll subroutine */ int covgen(i1, i2, d1, d2, d3, d4)

	double *d1,*d2,*d3,*d4;
	long *i1,*i2;

{
	extern int covgenf_(long*,long*,double*,double*,double*,double*);

	covgenf_(i1,i2,d1,d2,d3,d4);

	return 0;
}
