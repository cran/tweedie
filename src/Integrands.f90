MODULE Integrands_MOD
  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  USE Calcs_Imag, ONLY: evaluateImk
  USE Calcs_Real
  IMPLICIT NONE
  
  REAL(KIND=C_DOUBLE) :: Imk, Rek, lambda
  LOGICAL(C_BOOL)     :: error
  
CONTAINS

  FUNCTION Integrands(i, t) RESULT(integrand_result)
    ! Function to return the integrand values

    USE tweedie_params_mod
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
      
    IMPLICIT NONE
    
    INTEGER(C_INT), INTENT(IN)        :: i
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t                ! The internal variable for integration

    REAL(KIND=C_DOUBLE)               :: integrand_result ! The result of the function
     
    
    ! Grab the relevant scalar values for this iteration:
    current_y    = Cy(i)
    current_mu   = Cmu(i)


    ! Check for when t = 0 (t is very close to zero)
    IF (ABS(t) .LT. 1.0E-14_C_DOUBLE) THEN
      ! This should ideally be handled by the integrator (limits), 
      ! but returning the analytic limit is safest.
      integrand_result = current_mu - current_y
  
      RETURN
    ELSE
      CALL evaluateRek(t, Rek)
      CALL evaluateImk(t, Imk, error)
      IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, t, 1)
      
      IF (Cpdf) THEN
        IF (CpSmall) THEN
          CALL evaluateLambda(lambda)
          integrand_result = DEXP( Rek ) * DCOS( Imk ) - DEXP( -lambda ) * DCOS(t * current_y )
        ELSE
          integrand_result = DEXP( Rek ) * DCOS( Imk )
        END IF
      ELSE
        integrand_result = DEXP( Rek ) * DSIN( Imk ) / t
      END IF
    END IF
    
  END FUNCTION Integrands

END MODULE Integrands_MOD
