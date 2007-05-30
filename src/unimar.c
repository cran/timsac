#include "timsac.h"
extern int unimar(double*,long*,long*,double*,double*,double*,double*,double*,long*,double*,double*,double*,char**);

/* rtimsac78.dll subroutine */ int unimar(d1,i1,i2,d2,d3,d4,d5,d6,i3,d7,d8,d9,c1)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
	long *i1,*i2,*i3;
	char **c1;

{
	extern int unimarf_(double*,long*,long*,double*,double*,double*,double*,double*,long*,double*,double*,double*,long*);
	
	long *i4;
	i4 = (long*)*c1;

	unimarf_(d1,i1,i2,d2,d3,d4,d5,d6,i3,d7,d8,d9,i4);

	return 0;
}
