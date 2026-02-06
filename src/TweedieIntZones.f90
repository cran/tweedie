
MODULE TweedieIntZones
  
  USE tweedie_params_mod
  USE Calcs_K
  USE Calcs_Imag
  USE Calcs_Real
  USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_BOOL
  USE R_interfaces

  IMPLICIT NONE

CONTAINS 
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE checkStopPreAcc(tmax, zeroL, &
                             stop_PreAccelerate, converged_Pre)
    ! Determine if it is OK to stop pre-accelerating, and start using acceleration
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN) :: tmax
    REAL(KIND=C_DOUBLE), INTENT(IN) :: zeroL
     LOGICAL(C_BOOL), INTENT(OUT)    :: stop_PreAccelerate, converged_Pre
    
    ! Local vars
    INTEGER(C_INT)        :: nmax
    REAL(KIND=C_DOUBLE)   :: MM, Rek, Rekd, tstop, pi
    
    ! NOTE: 
    ! If  stop_PreAccelerate  is  .TRUE.  it means to stop pre-accelerating
    ! If  converged_Pre  is  .TRUE.  it means to convergence has been identified during pre-accelerating
    ! These are NOT necessarily the same.
    !
    ! If convergence is detected (converged_Pre is .TRUE.), then  stop_PreAccelerate  is also .TRUE.
    !
    ! But sometimes, pre-acceleration should stop (stop_PreAccelerate  is  .TRUE.)  even if
    ! convergence is NOT detected; for example:
    !  - the maximum number of iterations has been reached (so stop_PreAccelerate is .TRUE, 
    !    but  converged_Pre  is  .FALSE.); or
    ! - things are so well behaved that we can go straight to acceleration (stop_PreAccelerate is 
    !    .TRUE.  but  converged_Pre  is still  .FALSE.)
    
    
    pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
    

    ! Initialise
    stop_PreAccelerate = .FALSE.
    converged_Pre = .FALSE.

    ! Stop condition for pre-acceleration.
    ! Ensure that we have passed the peak of Im k(t), so that acceleration can be used
    
    IF ( Cp .GT. 2.0_C_DOUBLE) THEN
      IF (current_y .GT. current_mu) THEN
        ! Stop pre-accelerating, and go straight to acceleration 
        stop_PreAccelerate = .TRUE.
        converged_Pre = .FALSE.
      ELSE
        IF (zeroL .GT. tmax) THEN
          ! past tmax, we can usually stop pre-accelerating, and go straight to acceleration 
          stop_PreAccelerate = .TRUE.
        END IF
      END IF
    ELSE
      IF (current_y .GT. current_mu) THEN
        ! Stop pre-accelerating, and go straight to acceleration 
        stop_PreAccelerate = .TRUE.
      ELSE
        MM = 1.0_C_DOUBLE / (2.0_C_DOUBLE  * (Cp - 1.0_C_DOUBLE))
        IF ( MM .LE. 1.0_C_DOUBLE ) THEN
          stop_PreAccelerate = .TRUE.
        ELSE
          ! Check when t is larger than the last turning point of Re k(t) 
          ! Or when exp{Re(k)/t} is so small that it makes no difference... but 
          ! care is needed: Re k(t) is not necessarily convex here
           IF ( ABS( MM - FLOOR(MM) ) < 1.0E-09_C_DOUBLE ) then
              nmax = FLOOR(MM) - 1_C_INT
           ELSE
              nmax = FLOOR(MM)
           END IF
           tstop = current_mu**(1.0_C_DOUBLE - Cp) / ((1.0_C_DOUBLE - Cp) * current_phi) *   & 
                   DTAN( DBLE(nmax) * pi * (1.0_C_DOUBLE - Cp) )
           IF (zeroL .GT. tstop) stop_PreAccelerate = .TRUE.
        END IF
      END IF
    END IF
    
    ! Sometimes this takes forever to flag  stop_PreAccelerate  as .TRUE.
    ! so also check if exp{Re k(t)/t} is very small
    CALL evaluateRek( zeroL, Rek)
    CALL evaluateRekd(zeroL, Rekd)
    
    IF (zeroL .GT. 0.0_C_DOUBLE) THEN
      IF ( ( (DEXP(Rek)/zeroL) .LT. 1.0E-07_C_DOUBLE) .AND.          & 
          (Rekd .LT. 0.0_C_DOUBLE) ) then
        stop_PreAccelerate = .TRUE.
        converged_Pre = .TRUE.
        
        IF (Cverbose) THEN 
          ! CALL DBLEPR("checkStopPreAcc: Pre-accelerating stopping (a). Rek(t) small:", -1, (DEXP(Rek)/zeroL), 1)
        END IF
      END IF
    END IF

    IF (zeroL .GT. 0.0_C_DOUBLE) THEN
      IF ( ( (DEXP(Rek)/zeroL) .LT. 1.0E-15_C_DOUBLE)  ) THEN
        stop_PreAccelerate = .TRUE.
        converged_Pre = .TRUE.

        IF (Cverbose) THEN 
          ! CALL DBLEPR("checkStopPreAcc: Pre-accelerating stopping (b). Rek(t) small:", -1, (DEXP(Rek)/zeroL), 1)
        END IF
      END IF
    END IF
    
    ! If converged, then always stop pre-accelerating
    IF (converged_Pre) stop_PreAccelerate = .TRUE.

  END SUBROUTINE checkStopPreAcc
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE updateTM(i, tmax, mmax, left_Of_Max, &
                      m, zeroL, zeroR)
                  
    ! Update the values of  t  and  m  to the values needed
    ! for the next integration region.
    ! The updated value of  m  that is returned corresponds to the
    ! updated value of  zeroR.
    
    IMPLICIT NONE
  
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: tmax
    REAL(KIND=C_DOUBLE), INTENT(INOUT)  :: zeroR
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: zeroL
    INTEGER(C_INT), INTENT(IN)          :: i, mmax
    INTEGER(C_INT), INTENT(INOUT)       :: m
    LOGICAL(C_BOOL), INTENT(INOUT)      :: left_Of_Max
    
    ! Local vars
    INTEGER(C_INT)        :: mOld
    REAL(KIND=C_DOUBLE)   :: pi
    REAL(KIND=C_DOUBLE)   :: zeroBoundL, zeroBoundR, zeroStartPoint
    REAL(KIND=C_DOUBLE)   :: current_y, current_mu, current_phi

    
    pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)
    
    ! Grab the relevant scalar values for this iteration:
    current_y    = Cy(i)
    current_mu   = Cmu(i)
    current_phi  = Cphi(i)


    ! NEXT LEFT-SIDE ZERO: 
    ! Move the previous right-side zero to be the next left zero
    zeroL = zeroR
    ! Move the current m to mOld
    mOld = m
    
    ! NEXT RIGHT-SIDE ZERO
    ! Now work out the value of the next right zero
    
    ! - First: update to the next value of m
    CALL advanceM(m, mmax, mOld, left_Of_Max)

    ! - Secondly, find the next zero, corresponding to this value of m,
    !   which becomes  zeroR.
    !   So first find some bounds on this next zero.
    IF ( left_Of_Max ) THEN
      ! left_Of_Max is  TRUE  if the just-found value of m is to the left of mmax
      !
      ! If the value of m for the right-side zero is still to the left of tmax,
      ! leftOfMax is TRUE, and the upper bound is tmax. 
      zeroBoundR = tmax
      zeroBoundL = zeroR
    ELSE 
      ! If the value of zeroL is to the right of tmax,
      ! leftOfMax is FALSE, and the lower bound is tmax. 
      zeroBoundL = tmax
      zeroBoundR = zeroBoundL * 20.0_C_DOUBLE
    END IF
    
    ! With these bounds, we can now find the right-side zero,  zeroR
    ! Find a reasonable starting point for the algorithm:
    zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE

    ! Improve the starting point (sometimes very useful):
    CALL improveKZeroBounds(m, left_Of_Max, zeroStartPoint, &
                            zeroBoundL, zeroBoundR)
    zeroStartPoint = (zeroBoundL + zeroBoundR)/2.0_C_DOUBLE

    ! Now find the zero, within the bounds, with this starting point
    CALL findExactZeros(m, zeroBoundL, zeroBoundR, &
                        zeroStartPoint, zeroR, left_Of_Max)
    ! The zero just found  (zeroR)  is the right-side zero


    ! RETURNING: m, zeroL, zeroR
    
  END SUBROUTINE updateTM
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  SUBROUTINE findInitialZeroR(mfirst, left_Of_Max, tmax, &
                              zeroR)
  
    USE tweedie_params_mod
    
    IMPLICIT NONE
    
    REAL(KIND=C_DOUBLE), INTENT(OUT)    :: zeroR
    INTEGER(C_INT), INTENT(IN)          :: mfirst
    LOGICAL(C_BOOL), INTENT(INOUT)      :: left_Of_Max
    REAL(KIND=C_DOUBLE), INTENT(IN)     :: tmax

    ! Local vars
    REAL(KIND=C_DOUBLE)                 :: t_Start_Point, zeroBoundL, zeroBoundR
    REAL(KIND=C_DOUBLE)                 :: pi, TMP
    LOGICAL(C_BOOL)                     :: error


    pi = 4.0_C_DOUBLE * DATAN(1.0_C_DOUBLE)


    ! Initialisation
    t_Start_Point = 0.0_C_DOUBLE
    zeroBoundL = 0.0_C_DOUBLE
    zeroBoundR = 0.0_C_DOUBLE

    ! Find starting point for the first zero
    IF (left_Of_Max) THEN
      t_Start_Point = pi / current_y  
      zeroBoundL = 0_C_DOUBLE
      zeroBoundR = tmax   ! WAS: t_Start_Point * 2.0_C_DOUBLE
    ELSE
      ! Searching to the right of tmax
      t_Start_Point = tmax + pi / current_y  
      zeroBoundL = tmax
      zeroBoundR = t_Start_Point * 2.0_C_DOUBLE
    END IF

    IF ( (t_Start_Point .GT. zeroBoundR) .OR. (t_Start_Point .LT. zeroBoundL) ) Then
      t_Start_Point = (zeroBoundL + zeroBoundR) / 2.0_c_DOUBLE
    END IF
  
    ! Find the zero
    CALL findExactZeros(mfirst, zeroBoundL, zeroBoundR, t_Start_Point, zeroR, left_Of_Max)
    ! findExactZeros may change the value of  left_Of_Max

    CALL evaluateImk(zeroR, TMP, error)
    IF (error) CALL DBLEPR("ERROR: integrand zero =", -1, zeroR, 1)

  END SUBROUTINE findInitialZeroR

END MODULE TweedieIntZones

