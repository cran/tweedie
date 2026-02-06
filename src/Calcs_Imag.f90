
MODULE Calcs_Imag

  USE tweedie_params_mod
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE R_interfaces
  
  IMPLICIT NONE
  
  PUBLIC :: evaluateImk, evaluateImkM, evaluateImkd, evaluateImkdd, evaluateImkdZero
  
CONTAINS

  SUBROUTINE evaluateImk(t, Imk, error) 
    ! Evaluate Im k(t)
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: Imk
    LOGICAL(C_BOOL), INTENT(OUT)        :: error
    
    REAL(KIND=C_DOUBLE)   :: tanArg, omega, front, alpha, pi

  
    ! Initialisation
    error = .FALSE.
    pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
    front = current_mu ** (2.0_C_DOUBLE - Cp) / ( current_phi * (2.0_C_DOUBLE - Cp))
    tanArg = (1.0_C_DOUBLE - Cp) * t * current_phi  / (current_mu ** (1.0_C_DOUBLE - Cp) )
    omega = DATAN( tanArg )
  
    IF ((omega .GT. 0.0_C_DOUBLE ) .OR.    &    
        (omega .LT. (-pi/2.0_C_DOUBLE)) ) THEN
      ! CALL DBLEPR("ERROR (evaluateImk): Omega out of range:", -1, omega, 1)
      ! CALL DBLEPR("ERROR (evaluateImk): t:", -1, t, 1)
      ! CALL DBLEPR("         p:", -1, Cp, 1)
      ! CALL DBLEPR("         y:", -1, current_y, 1)
      ! CALL DBLEPR("        mu:", -1, current_mu, 1)
      ! CALL DBLEPR("       phi:", -1, current_phi, 1)

      error = .TRUE.
      RETURN
    END IF
    alpha = (2.0_C_DOUBLE - Cp)/(1.0_C_DOUBLE - Cp)
  
    Imk = front *   &
          DSIN(omega * alpha)/(DCOS(omega) ** alpha) - t * current_y
  
    RETURN
  
  END SUBROUTINE evaluateImk
  
  
  

  SUBROUTINE evaluateImkM(t, f, df, m)
    ! Evaluate  Im k(t) - m*pi - pi/2  (PDF) or  Im k(t) - m*pi  (CDF) for 
    ! finding the zeros of the integrand, for the PDf and CDF
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: f, df
    INTEGER(C_INT), INTENT(IN)        :: m
  
    REAL(KIND=C_DOUBLE)               :: pi, Imk_val
    LOGICAL(C_BOOL)                   :: error
  
    pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
    CALL evaluateImk(t, Imk_val, error)
    IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, t, 1)

    ! The expression depends on whether we are working with the PDF or the CDF.
    ! The PDF uses cos Im k(t) in the integrand; the CDF has sin Im k(t) in the integrand.
    ! Thus, the PDF has integrand zeros at Im k(t) = pi/2 + m pi/y;
    !       the CDF has integrand zeros at Im k(t) =        m pi/y.
    IF (Cpdf) THEN
      f = Imk_val - REAL(m, KIND=C_DOUBLE) * pi - pi/2.0_C_DOUBLE
    ELSE
      f = Imk_val - REAL(m, KIND=C_DOUBLE) * pi
    END IF
    CALL evaluateImkd(t, df)
    
  END SUBROUTINE evaluateImkM




  SUBROUTINE evaluateImkd(t, Imkd) 
    ! Evaluate Im k'(t)
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imkd  ! The result of the calculation
    
    REAL(KIND=C_DOUBLE) :: omega, pindex

  
    pindex = 1.0_C_DOUBLE / (1.0_C_DOUBLE - Cp)
    omega = DATAN( ( (1.0_C_DOUBLE - Cp) * t * current_phi) / (current_mu ** (1.0_C_DOUBLE - Cp) ) )
    
    Imkd = current_mu * (DCOS(omega * pindex) / (DCOS(omega) ** pindex)) - current_y
  
  END SUBROUTINE evaluateImkd





  SUBROUTINE evaluateImkdd(t, Imkdd)
    ! Evaluate Im k''(t)
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)   :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)  :: Imkdd
    
    REAL(KIND=C_DOUBLE)    :: front, omega, pindex


    pindex = Cp / (1.0_C_DOUBLE - Cp)
    front = -current_phi * current_mu ** (Cp/(1.0_C_DOUBLE - Cp))
    omega = DATAN( ( (1.0_C_DOUBLE - Cp) * t * current_phi) / (current_mu ** (1.0_C_DOUBLE - Cp) ) )
  
    Imkdd = front * (DSIN(omega * pindex) / (DCOS(omega) ** pindex) )
  
  END SUBROUTINE evaluateImkdd
  

  
  
  SUBROUTINE evaluateImkdZero(t, f, df) 
  ! Evaluate Im k'(t)  and  Im k''(t)  for solving for Kmax (i.e., Im k'(t) = 0)
    
    USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
  
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: t
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: f, df
  
    REAL(KIND=C_DOUBLE)  :: Imkd, Imkdd
  
    CALL evaluateImkd( t, Imkd)
    CALL evaluateImkdd(t, Imkdd)
  
    f  = Imkd
    df = Imkdd
  
  END SUBROUTINE evaluateImkdZero



END MODULE Calcs_Imag

