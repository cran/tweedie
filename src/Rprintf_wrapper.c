#include <R.h>

// Print a double with a label
void rprintf_double(const char *label, double *x) {
    Rprintf("%s %f\n", label, *x);
}

// Print an integer with a label
void rprintf_int(const char *label, int *x) {
    Rprintf("%s %d\n", label, *x);
}
