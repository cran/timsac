#include "timsac.h"
extern int auspec(long*, long*, double*, double*, double*, double*);

/* rtimsac72.dll subroutine */ int auspec(i1,i2,d1,d2,d3,d4)

	double *d1,*d2,*d3,*d4;
	long *i1,*i2;

{
	extern int auspecf_(long*, long*, double*, double*, double*, double*);

	auspecf_(i1,i2,d1,d2,d3,d4);

	return 0;
}
