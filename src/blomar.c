#include <R.h>
#include "timsac.h"
//#include <R_ext/RS.h>

extern int blomar(double*,long*,long*,double*,long*,long*,long*,double*,double*,double*,double*,double*,double*,double*,long*,long*);

/* rtimsac78.dll subroutine */ int blomar(d1,i1,i2,d2,i3,i4,i5,d3,d4,d5,d6,d7,d8,d9,i6,i7)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
	long *i1,*i2,*i3,*i4,*i5,*i6,*i7;

{

	extern int F77_NAME(blomarf) (double*,long*,long*,double*,long*,long*,long*,double*,double*,double*,double*,double*,double*,double*,long*,long*);

	F77_CALL(blomarf) (d1,i1,i2,d2,i3,i4,i5,d3,d4,d5,d6,d7,d8,d9,i6,i7);

	return 0;

}


