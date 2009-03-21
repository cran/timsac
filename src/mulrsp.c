#include "timsac.h"
#include <R_ext/RS.h>

struct Rcomplex {
	double r;
	double i;
};

extern int mulrsp(long*, long*, long*, long*, double*, double*, double*, struct Rcomplex*, double*);

/* rtimsac72.dll subroutine */ int mulrsp(i1,i2,i3,i4,d1,d2,d3,d4,d5)

	double *d1,*d2,*d3,*d5;
	long *i1,*i2,*i3,*i4;
	struct Rcomplex *d4;

{
	extern int F77_NAME(mulrspf) (long*, long*, long*, long*, double*, double*, double*, struct Rcomplex*, double*);

	F77_CALL(mulrspf) (i1,i2,i3,i4,d1,d2,d3,d4,d5);

    return 0;
}
