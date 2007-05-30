#include "timsac.h"
extern int autcor(double*, long*, double*, double*, long*, double*);

/* rtimsac72.dll subroutine */ int autcor(d1,i1,d2,d3,i2,d4)

	double *d1,*d2,*d3,*d4;
	long *i1,*i2;

{
	extern int autcorf_(double*, long*, double*, double*, long*, double*);

	autcorf_(d1,i1,d2,d3,i2,d4);

	return 0;
}
