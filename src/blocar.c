#include "timsac.h"
extern int blocar(double*,long*,long*,long*,long*,double*,double*,double*,double*,double*,double*,double*,long*,long*,double*);

/* rtimsac78.dll subroutine */ int blocar(d1,i1,i2,i3,i4,d2,d3,d4,d5,d6,d7,d8,i5,i6,d9)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
	long *i1,*i2,*i3,*i4,*i5,*i6;

{
	extern int blocarf_(double*,long*,long*,long*,long*,double*,double*,double*,double*,double*,double*,double*,long*,long*,double*);

	blocarf_(d1,i1,i2,i3,i4,d2,d3,d4,d5,d6,d7,d8,i5,i6,d9);

	return 0;
}
