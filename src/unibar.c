#include "timsac.h"
extern int unibar(double*,long*,long*,double*,double*,double*,double*,double*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

/* rtimsac78.dll subroutine */ int unibar(d1,i1,i2,d2,d3,d4,d5,d6,i3,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16,*d17;
	long *i1,*i2,*i3;

{
	extern int unibarf_(double*,long*,long*,double*,double*,double*,double*,double*,long*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

	unibarf_(d1,i1,i2,d2,d3,d4,d5,d6,i3,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17);

	return 0;
}
