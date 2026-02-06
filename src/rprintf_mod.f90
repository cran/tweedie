MODULE rprintf_mod
  ! For printing at the R console from within FORTRAN
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE rprintf_double(label, x) BIND(C, NAME="rprintf_double")
      USE ISO_C_BINDING
      CHARACTER(KIND=C_CHAR), DIMENSION(*) :: label
      REAL(C_DOUBLE), INTENT(IN) :: x
    END SUBROUTINE rprintf_double

    SUBROUTINE rprintf_int(label, x) BIND(C, NAME="rprintf_int")
      USE ISO_C_BINDING
      CHARACTER(KIND=C_CHAR), DIMENSION(*) :: label
      INTEGER(C_INT), INTENT(IN) :: x
    END SUBROUTINE rprintf_int
  END INTERFACE

END MODULE rprintf_mod
