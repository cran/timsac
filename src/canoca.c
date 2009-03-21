#include "timsac.h"
#include <R_ext/RS.h>

extern int canoca(long*,long*,long*,long*,long*,double*,long*,double*,double*,long*,double*,double*,long*,long*,long*,double*,double*,double*,double*,long*,double*,double*,long*,double*,long*,long*,double*,long*,double*,long*,long*,long*);

/* rtimsac74.dll subroutine */ int canoca(i1,i2,i3,i4,i5,d1,i6,d2,d3,i7,d4,d5,i8,i9,i10,d6,d7,d8,d9,i11,d10,d11,i12,d12,i13,i14,d13,i15,d14,i16,i17,i18)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13,*i14,*i15,*i16,*i17,*i18;

{
	extern int F77_NAME(canocaf) (long*,long*,long*,long*,long*,double*,long*,double*,double*,long*,double*,double*,long*,long*,long*,double*,double*,double*,double*,long*,double*,double*,long*,double*,long*,long*,double*,long*,double*,long*,long*,long*);

	F77_CALL(canocaf) (i1,i2,i3,i4,i5,d1,i6,d2,d3,i7,d4,d5,i8,i9,i10,d6,d7,d8,d9,i11,d10,d11,i12,d12,i13,i14,d13,i15,d14,i16,i17,i18);

	return 0;
}
