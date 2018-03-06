#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP AuspecC(SEXP, SEXP, SEXP);
extern SEXP AutarmC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP AutcorC(SEXP, SEXP, SEXP);
extern SEXP BayseaC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BispecC(SEXP, SEXP, SEXP, SEXP);
extern SEXP BlocarC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BlomarC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BsubstC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CanarmC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CanocaC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CovgenC(SEXP, SEXP, SEXP, SEXP);
extern SEXP DecompC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ExsarC(SEXP, SEXP, SEXP);
extern SEXP FftcorC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FpeautC(SEXP, SEXP, SEXP, SEXP);
extern SEXP Fpec7C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MarkovC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MlocarC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MlomarC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MulbarC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MulcorC(SEXP, SEXP, SEXP, SEXP);
extern SEXP MulfrfC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MulmarC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MulnosC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MulrspC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MulspeC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NonstC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP OptdesC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP OptsimC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PerarsC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PrdctrC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RaspecC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP SglfreC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP SimconC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP SpgrhC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ThirmoC(SEXP, SEXP, SEXP);
extern SEXP UnibarC(SEXP, SEXP, SEXP);
extern SEXP UnimarC(SEXP, SEXP, SEXP);
extern SEXP WnoiseC(SEXP, SEXP, SEXP);
extern SEXP XsarmaC(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"AuspecC", (DL_FUNC) &AuspecC,  3},
    {"AutarmC", (DL_FUNC) &AutarmC, 11},
    {"AutcorC", (DL_FUNC) &AutcorC,  3},
    {"BayseaC", (DL_FUNC) &BayseaC, 11},
    {"BispecC", (DL_FUNC) &BispecC,  4},
    {"BlocarC", (DL_FUNC) &BlocarC,  5},
    {"BlomarC", (DL_FUNC) &BlomarC,  7},
    {"BsubstC", (DL_FUNC) &BsubstC,  8},
    {"CanarmC", (DL_FUNC) &CanarmC,  6},
    {"CanocaC", (DL_FUNC) &CanocaC,  9},
    {"CovgenC", (DL_FUNC) &CovgenC,  4},
    {"DecompC", (DL_FUNC) &DecompC,  5},
    {"ExsarC",  (DL_FUNC) &ExsarC,   3},
    {"FftcorC", (DL_FUNC) &FftcorC,  7},
    {"FpeautC", (DL_FUNC) &FpeautC,  4},
    {"Fpec7C",  (DL_FUNC) &Fpec7C,   7},
    {"MarkovC", (DL_FUNC) &MarkovC, 14},
    {"MlocarC", (DL_FUNC) &MlocarC,  6},
    {"MlomarC", (DL_FUNC) &MlomarC,  8},
    {"MulbarC", (DL_FUNC) &MulbarC,  5},
    {"MulcorC", (DL_FUNC) &MulcorC,  4},
    {"MulfrfC", (DL_FUNC) &MulfrfC,  6},
    {"MulmarC", (DL_FUNC) &MulmarC,  5},
    {"MulnosC", (DL_FUNC) &MulnosC,  5},
    {"MulrspC", (DL_FUNC) &MulrspC,  7},
    {"MulspeC", (DL_FUNC) &MulspeC,  5},
    {"NonstC",  (DL_FUNC) &NonstC,   5},
    {"OptdesC", (DL_FUNC) &OptdesC,  9},
    {"OptsimC", (DL_FUNC) &OptsimC,  8},
    {"PerarsC", (DL_FUNC) &PerarsC,  5},
    {"PrdctrC", (DL_FUNC) &PrdctrC, 13},
    {"RaspecC", (DL_FUNC) &RaspecC,  6},
    {"SglfreC", (DL_FUNC) &SglfreC,  6},
    {"SimconC", (DL_FUNC) &SimconC,  9},
    {"SpgrhC",  (DL_FUNC) &SpgrhC,   6},
    {"ThirmoC", (DL_FUNC) &ThirmoC,  3},
    {"UnibarC", (DL_FUNC) &UnibarC,  3},
    {"UnimarC", (DL_FUNC) &UnimarC,  3},
    {"WnoiseC", (DL_FUNC) &WnoiseC,  3},
    {"XsarmaC", (DL_FUNC) &XsarmaC,  5},
    {NULL, NULL, 0}
};

void R_init_timsac(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

