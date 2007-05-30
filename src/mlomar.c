#include "timsac.h"
extern int mlomar(double*,long*,long*,double*,long*,long*,long*,long*,double*,double*,long*,long*,long*,double*,long*,double*,long*,double*,double*,double*,long*,long*);

/* rtimsac78.dll subroutine */ int mlomar(d1,i1,i2,d2,i3,i4,i5,i6,d3,d4,i7,i8,i9,d5,i10,d6,i11,d7,d8,d9,i12,i13)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13;

{

	extern int mlomarf_(double*,long*,long*,double*,long*,long*,long*,long*,double*,double*,long*,long*,long*,double*,long*,double*,long*,double*,double*,double*,long*,long*);

	mlomarf_(d1,i1,i2,d2,i3,i4,i5,i6,d3,d4,i7,i8,i9,d5,i10,d6,i11,d7,d8,d9,i12,i13);

	return 0;

}
