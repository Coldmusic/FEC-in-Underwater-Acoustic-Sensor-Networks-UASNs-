SUBROUTINE SUBTABd( x, Nx )

  ! If x(3) = -999.9 then subtabulation is performed
  ! i.e., a vector is generated with Nx points in [ x(1), x(2) ]
  ! If x(2) = -999.9 then x(1) is repeated

  IMPLICIT NONE
  INTEGER, INTENT( IN ) :: Nx
  INTEGER               :: I
  REAL (KIND=8)         :: x( Nx ), deltax

  IF ( Nx >= 3 ) THEN
     IF ( x( 3 ) == -999.9 ) THEN   ! testing for equality here is dangerous
        IF ( x( 2 ) == -999.9 ) x( 2 ) = x( 1 )
        deltax      = ( x( 2 ) - x( 1 ) ) / ( Nx - 1 )
        x( 1 : Nx ) = x( 1 ) + [ ( I, I = 0, Nx - 1 ) ] * deltax
     END IF
  END IF

END SUBROUTINE SUBTABd
