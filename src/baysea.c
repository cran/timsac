#include <R_ext/RS.h>
#include "timsac.h"

extern int baysea(double*, long*, long*, double*, double*, double*, double*, double*, double*, double*,
				  double*, double*, double*, double*, long*, double*, double*, double*, double*, 
				  long*, long*, long*);

int baysea(d1,i1,i2,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,i3,d13,d14,d15,d16,i4,i5,i6)

	double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15,*d16;
	long *i1,*i2,*i3,*i4,*i5,*i6;

{
	extern int F77_NAME(bayseaf) (double*, long*, long*, double*, double*, double*, double*, double*, double*, double*,
								  double*, double*, double*, double*, long*, double*, double*, double*, double*,
								  long*, long*, long*);

	F77_CALL(bayseaf) (d1,i1,i2,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,i3,d13,d14,d15,d16,i4,i5,i6);

	return 0;
}

