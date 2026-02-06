#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// DECLARATION: This declares the function twcomputation as a standard C function.
// Since the Fortran routine uses BIND(C), it exports this exact symbol name
// with the C calling convention (no hidden arguments).
extern void twcomputation(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *); 
// Arguments are declared as void* to match R's .C interface for simple vectors.

/* .C calls */
static const R_CMethodDef CEntries[] = {
  // REGISTRATION: The name MUST be exactly "twcomputation"
  {"twcomputation", (DL_FUNC) &twcomputation, 11},
  {NULL, NULL, 0}
};

// Fortran entries table MUST be NULL or empty if you removed the F77_NAME definition
static const R_FortranMethodDef FortranEntries[] = {
  {NULL, NULL, 0}
};

void R_init_tweedie(DllInfo *dll)
{
  // Register only the CEntries table
  R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
