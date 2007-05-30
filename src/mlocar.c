#include "timsac.h"
extern int mlocar(double*,long*,long*,long*,long*,long*,double*,double*,double*,long*,double*,long*,long*,double*,long*,long*,long*,double*,double*,long*,double*,double*);

/* rtimsac78.dll subroutine */ int mlocar(d1,i1,i2,i3,i4,i5,d2,d3,d4,i6,d5,i7,i8,d6,i9,i10,i11,d7,d8,i12,d9,d10)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12;

{
	extern int mlocarf_(double*,long*,long*,long*,long*,long*,double*,double*,double*,long*,double*,long*,long*,double*,long*,long*,long*,double*,double*,long*,double*,double*);

	mlocarf_(d1,i1,i2,i3,i4,i5,d2,d3,d4,i6,d5,i7,i8,d6,i9,i10,i11,d7,d8,i12,d9,d10);

	return 0;
}
