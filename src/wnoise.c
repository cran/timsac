#include "timsac.h"
extern int wnoise(long*, long*, double*, double*);

/* rtimsac72.dll subroutine */	int wnoise(i1,i2,d1,d2)

	double *d1,*d2;
	long *i1,*i2;

{
	extern int wnoisef_(long*, long*, double*, double*);

	wnoisef_(i1,i2,d1,d2);

	return 0;
}
