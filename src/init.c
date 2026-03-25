#include "regF77.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* .Fortran calls */

static const R_FortranMethodDef FortEntries[] = {
    {"auspecf",   (DL_FUNC) &F77_NAME(auspecf),    6},
    {"autcorf",   (DL_FUNC) &F77_NAME(autcorf),    6},
    {"fftcorf",   (DL_FUNC) &F77_NAME(fftcorf),   14},
    {"fpeautf",   (DL_FUNC) &F77_NAME(fpeautf),   16},
    {"fpec7",    (DL_FUNC) &F77_NAME(fpec7),    17},	
    {"mulcorf",   (DL_FUNC) &F77_NAME(mulcorf),    7},
    {"mulfrff",   (DL_FUNC) &F77_NAME(mulfrff),   14},
    {"mulnosf",   (DL_FUNC) &F77_NAME(mulnosf),    8},
    {"mulrspf",   (DL_FUNC) &F77_NAME(mulrspf),    9},
    {"mulspef",   (DL_FUNC) &F77_NAME(mulspef),   10},
    {"optdesf",   (DL_FUNC) &F77_NAME(optdesf),   10},
    {"optsimf",   (DL_FUNC) &F77_NAME(optsimf),   18},
    {"raspecf",   (DL_FUNC) &F77_NAME(raspecf),    7},
    {"sglfref",   (DL_FUNC) &F77_NAME(sglfref),   16},
    {"wnoisef",   (DL_FUNC) &F77_NAME(wnoisef),    4},
	
    {"autarm",   (DL_FUNC) &F77_NAME(autarm),   23},
    {"bispecf",   (DL_FUNC) &F77_NAME(bispecf),   11},
    {"canarmf",   (DL_FUNC) &F77_NAME(canarmf),   27},
    {"canocaf",   (DL_FUNC) &F77_NAME(canocaf),   32},
    {"covgenf",   (DL_FUNC) &F77_NAME(covgenf),    6},
    {"markovf",   (DL_FUNC) &F77_NAME(markovf),   30},
    {"nonstf",    (DL_FUNC) &F77_NAME(nonstf),    14},
    {"prdctrf",   (DL_FUNC) &F77_NAME(prdctrf),   23},
    {"simconf",   (DL_FUNC) &F77_NAME(simconf),   15},
    {"thirmof",   (DL_FUNC) &F77_NAME(thirmof),    7},
	
    {"blocarf",   (DL_FUNC) &F77_NAME(blocarf),   15},
    {"blomarf",   (DL_FUNC) &F77_NAME(blomarf),   17},
    {"bsubstf",   (DL_FUNC) &F77_NAME(bsubstf),   37},
	{"exsarf",    (DL_FUNC) &F77_NAME(exsarf),    15},
    {"mlocarf",   (DL_FUNC) &F77_NAME(mlocarf),   22},
    {"mlomarf",   (DL_FUNC) &F77_NAME(mlomarf),   23},
    {"mulbarf",   (DL_FUNC) &F77_NAME(mulbarf),   21},
    {"mulmarf",   (DL_FUNC) &F77_NAME(mulmarf),   24},
    {"perarsf",   (DL_FUNC) &F77_NAME(perarsf),   16},
    {"unibarf",   (DL_FUNC) &F77_NAME(unibarf),   20},
    {"unimarf",   (DL_FUNC) &F77_NAME(unimarf),   12},
    {"xsarmaf",   (DL_FUNC) &F77_NAME(xsarmaf),   13},
	
    {"bayseaf",   (DL_FUNC) &F77_NAME(bayseaf),   22},
    {"decompf",   (DL_FUNC) &F77_NAME(decompf),   12},
    {"spgrh",    (DL_FUNC) &F77_NAME(spgrh),    14},
    {NULL, NULL, 0}
};

void attribute_visible R_init_timsac(DllInfo *dll) 
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
