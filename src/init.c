/* 
   Generated from: tools::package_native_routine_registration_skeleton("tweedie", character_only = FALSE)
*/
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(twcdf)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(twpdf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"twcdf", (DL_FUNC) &F77_NAME(twcdf),  9},
  {"twpdf", (DL_FUNC) &F77_NAME(twpdf), 10},
  {NULL, NULL, 0}
};

void R_init_tweedie(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
