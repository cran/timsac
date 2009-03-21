#include "timsac.h"
#include <R_ext/RS.h>

extern int canarm(long*,long*,double*,double*,long*,double*,double*,double*,long*,double*,long*,long*,long*,double*,double*,double*,double*,long*,double*,double*,long*,long*,double*,long*,double*,char**,long*,long*);

/* rtimsac74.dll subroutine */ int canarm(i1,i2,d1,d2,i3,d3,d4,d5,i4,d6,i5,i6,i7,d7,d8,d9,d10,i8,d11,d12,i9,i10,d13,i11,d14,c1,i12,i13)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13;
	char **c1;

{
	extern int F77_NAME(canarmf) (long*,long*,double*,double*,long*,double*,double*,double*,long*,double*,long*,long*,long*,double*,double*,double*,double*,long*,double*,double*,long*,long*,double*,long*,double*,long*,long*,long*);

	long *ic;
	ic = (long*)*c1;

	F77_CALL(canarmf) (i1,i2,d1,d2,i3,d3,d4,d5,i4,d6,i5,i6,i7,d7,d8,d9,d10,i8,d11,d12,i9,i10,d13,i11,d14,ic,i12,i13);

	return 0;
}

