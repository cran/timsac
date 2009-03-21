#include "timsac.h"
#include <R_ext/RS.h>

extern int ngsmth(double*, long*, long*, double*, double*, long*, double*, double*, long*, double*, float*, double*, long*, long*, long*, long*);

int ngsmth(d1,i1,i2,d2,d3,i3,d4,d5,i4,d6,f1,d7,i5,i6,i7,i8)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7;
	float  *f1;
	long   *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;

{
	extern int F77_NAME(ngsmthf) (double*, long*, long*, double*, double*, long*, double*, double*, long*, double*, float*, double*, long*, long*, long*, long*);

	F77_CALL(ngsmthf) (d1,i1,i2,d2,d3,i3,d4,d5,i4,d6,f1,d7,i5,i6,i7,i8);

	return 0;
}
