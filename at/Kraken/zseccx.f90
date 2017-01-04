SUBROUTINE ZSECCX( x2, Tolerance, Iteration, MaxIteration, ErrorMessage )

  ! Secant method
  ! Input:
  ! x2 is an initial guess
  ! Tolerance is an error bound for the returned root
  ! MaxIteration is the maximum allowable number of iterations
  !
  ! Output:
  ! x2 is the root
  ! Iteration is the number of iterations used
  ! ErrorMessage is empty unless there was a problem.

  IMPLICIT NONE
  INTEGER, INTENT(  IN )               :: MaxIteration
  INTEGER, INTENT( OUT )               :: Iteration
  INTEGER                              :: IPower0, IPower1
  REAL    ( KIND=8 ), INTENT( IN    )  :: Tolerance
  COMPLEX ( KIND=8 ), INTENT( INOUT )  :: x2
  COMPLEX ( KIND=8 )                   :: x0, x1, shift, f0, f1, CNum, CDen
  CHARACTER (LEN=80), INTENT( OUT    ) :: ErrorMessage 

  ErrorMessage = ' ' 
  IF ( Tolerance <= 0.0D0 ) THEN 
     ErrorMessage = 'Non-positive tolerance specified' 
     STOP 
  ENDIF

  x1     = x2 + 100.0D0 * Tolerance 
  CALL FUNCT( x1, f1, IPower1 ) 

  ! WRITE( *, * )
  ! WRITE( *, "( 4G17.9, I5 )" ) SQRT( x1 ), f1, IPower1

  DO Iteration = 1, MaxIteration 
     x0      = x1
     f0      = f1 
     IPower0 = IPower1 

     x1 = x2 
     CALL FUNCT( x1, f1, IPower1 )

     ! ugly stuff to block overflows by forcing shift to be bounded
     CNum = f1 * ( x1 - x0 )
     CDen = f1 - f0 * 10.0D0 ** ( IPower0 - IPower1 )

     IF ( ABS( CNum ) >= ABS( CDen * x1 ) ) THEN 
        shift = 0.1D0 * Tolerance
     ELSE
        shift = CNum / CDen 
     ENDIF

     x2 = x1 - shift 
     ! WRITE( *, "( 6G24.16, I5 )" ) SQRT( x1 ), x1, f1, IPower1

     ! The following convergence test is very restrictive
     ! Found it was necessary to check all 3 points (x0, x1, x2) for certain problems ...
     IF ( ABS( x2 - x1 ) + ABS( x2 - x0 ) < Tolerance ) then
        ! write( *, * ) 'converged', x0, x1, x2, Tolerance
        RETURN
     END IF
  ENDDO

  ErrorMessage = ' *** FAILURE TO CONVERGE IN SECANT'

END SUBROUTINE ZSECCX
