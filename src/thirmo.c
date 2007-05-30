#include "timsac.h"
extern int thirmo(long*,long*,double*,double*,double*,double*,double*);

/* rtimsac74.dll subroutine */ int thirmo(i1,i2,d1,d2,d3,d4,d5)

	double *d1,*d2,*d3,*d4,*d5;
	long *i1,*i2;

{
	extern int thirmof_(long*,long*,double*,double*,double*,double*,double*);

	thirmof_(i1,i2,d1,d2,d3,d4,d5);

	return 0;
}
