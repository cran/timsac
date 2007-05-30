#include "timsac.h"
extern int mulcor(double*, long*, long*, long*, double* ,double*, double*);

/* rtimsac72.dll subroutine */ int mulcor (d1,i1,i2,i3,d2,d3,d4)

	double *d1,*d2,*d3,*d4;
	long *i1,*i2,*i3;

{
	extern int mulcorf_(double*, long*, long*, long*, double*, double*, double*);

	mulcorf_(d1,i1,i2,i3,d2,d3,d4);

	return 0;
}
