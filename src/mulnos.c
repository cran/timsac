#include "timsac.h"
extern int mulnos(long*, long*, long*, double*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */ int mulnos(i1,i2,i3,d1,d2,d3,d4,d5)

	double *d1,*d2,*d3,*d4,*d5;
	long *i1,*i2,*i3;

{
	extern int mulnosf_(long*, long*, long*, double*, double*, double*, double*, double*);

	mulnosf_(i1,i2,i3,d1,d2,d3,d4,d5);

	return 0;
}
