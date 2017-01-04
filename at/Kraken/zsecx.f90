SUBROUTINE ZSECX( x2, Tolerance, Iteration, MaxIteration, ErrorMessage ) 

  ! Secant method
  ! Input:
  !   x2 is an initial guess
  !   Tolerance is an error bound for the returned root
  !   MaxIteration is the maximum allowable number of iterations
  !
  ! Output:
  !   x2 is the root
  !   ErrorMessage is empty unless there was a problem.

  IMPLICIT NONE
  INTEGER,          INTENT( IN    ) :: MaxIteration
  INTEGER,          INTENT( OUT   ) :: Iteration
  INTEGER                           :: IPower0, IPower1
  REAL ( KIND=8 ),  INTENT( IN    ) :: Tolerance
  REAL ( KIND=8 ),  INTENT( INOUT ) :: x2
  REAL ( KIND=8 )                   :: x0, x1, shift, f0, f1, cNum, cDen
  CHARACTER (LEN=80), INTENT( OUT ) :: ErrorMessage 

  ! Secant method                                                     

  ErrorMessage = ' ' 
  x1 = x2 + 10.0D0 * Tolerance

  CALL FUNCT( x1, F1, IPower1 )
  ! WRITE( *, * )
  ! WRITE( *, FMT="( 2G24.16, I5 )" ) SQRT( x1 ), F1, IPower1

  DO Iteration = 1, MaxIteration 
     x0      = x1
     F0      = F1 
     IPower0 = IPower1 
     x1      = x2 

     CALL FUNCT( x1, F1, IPower1 )

!!$     IF ( F1 == 0.0D0 ) THEN 
!!$        shift = 0.0D0
!!$     ELSE 
!!$        shift = ( x1 - x0 ) / ( 1.0D0 - F0 / F1 * 10.0D0 ** ( IPower0 - IPower1 ) )
!!$     ENDIF

     ! ugly stuff to block overflows by forcing shift to be bounded
     CNum = f1 * ( x1 - x0 )
     CDen = f1 - f0 * 10.0D0 ** ( IPower0 - IPower1 )

     IF ( ABS( CNum ) >= ABS( CDen * x1 ) ) THEN 
        shift = 0.1D0 * Tolerance
     ELSE
        shift = CNum / CDen 
     ENDIF

     x2 = x1 - shift
     ! WRITE( *, FMT="( 3G24.16, I5, G24.16 )" ) SQRT( x1 ), x1, F1, IPower1, x2

     ! The following convergence test is very restrictive
     ! Found it was necessary to check all 3 points (x0, x1, x2) for certain problems ...
     IF ( ABS( x2 - x1 ) + ABS( x2 - x0 ) < Tolerance ) THEN
        ! WRITE( *, * ) 'converged', x0, x1, x2, Tolerance
        RETURN
     END IF

  END DO

  ErrorMessage = ' *** FAILURE TO CONVERGE IN SECANT' 

END SUBROUTINE ZSECX
