#include "timsac.h"
extern int markov(long*,long*,long*,double*,long*,long*,long*,double*,double*,long*,long*,long*,long*,long*,long*,double*,double*,double*,double*,double*,long*,double*,double*,double*,double*,double*,char**,long*,long*,long*,long*);

/* rtimsac74.dll subroutine */ int markov(i1,i2,i3,d1,i4,i5,i6,d2,d3,i7,i8,i9,i10,i11,i12,d4,d5,d6,d7,d8,i13,d9,d10,d11,d12,d13,c1,i14,i15,i16,i17)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13,*i14,*i15,*i16,*i17;
	char **c1;

{
	extern int markovf_(long*,long*,long*,double*,long*,long*,long*,double*,double*,long*,long*,long*,long*,long*,long*,double*,double*,double*,double*,double*,long*,double*,double*,double*,double*,double*,long*,long*,long*,long*,long*);

	long *i18;
	i18 = (long*)*c1;

	markovf_(i1,i2,i3,d1,i4,i5,i6,d2,d3,i7,i8,i9,i10,i11,i12,d4,d5,d6,d7,d8,i13,d9,d10,d11,d12,d13,i18,i14,i15,i16,i17);

	return 0;
}
