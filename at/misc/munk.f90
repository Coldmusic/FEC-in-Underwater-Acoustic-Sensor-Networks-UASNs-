SUBROUTINE ANALYT( Depth, cP, cS, rho, Medium, N1, Freq, ATUNIT, Task )

  !     Munk profile

  !     Returns
  !        cS, cP, rho at depths i*h i = 1, N
  !        Depths of interfaces

  IMPLICIT NONE
  REAL (KIND=8), PARAMETER   :: eps = 0.00737
  INTEGER           :: i, Depth, Freq, Medium, N, N1
  REAL     (KIND=8) :: h, x, z
  COMPLEX  (KIND=8) :: cP( N1 ), cS( N1 )
  REAL     (KIND=8) :: rho( N1 )
  CHARACTER (LEN=1) :: ATUNIT
  CHARACTER (LEN=8) :: Task

  !     Initialization is only important for providing a print-out of
  !     the SSP for the user.  Unfortunately it overruns the array
  !     rho( 1 ) so we suppress it:

  IF ( Task( 1 : 4 ) == 'INIT' ) RETURN

  N = N1 - 1

  SELECT CASE ( Medium )

  CASE ( 1 )   ! THE OCEAN
     h = 5000.0 / N
     DO i = 1, N1
        z = ( i - 1 ) * h
        x = 2.0 * ( z - 1300.0 ) / 1300.0
        cP(  i ) = 1500.0 * ( 1.0 + eps * ( x - 1.0 + EXP( -x ) ) )
        cS(  i ) = 0.0
        rho( i ) = 1.0
     END DO
     RETURN

  CASE ( 2 )   ! THE FLUID HALF-SPACE
     cP(  1 ) = 1551.91
     cS(  1 ) = 0.0
     rho( 1 ) = 1.0E20
     RETURN

  CASE( 9 )  ! AN ELASTIC LAYER
     h = 1000.0 / N
     z = 5000.0

     DO i = 1,N+1
        cP( i ) = 4700.0 + ( z - 5000.0 ) / 10.0
        cS( i ) = 2000.0 + ( z - 5000.0 ) / 10.0
        cP( i ) = 4700.0
        cS( i ) = 2000.0
        rho( i ) = 2.0
        z = z + h
     END DO
  END SELECT

END SUBROUTINE ANALYT
