SUBROUTINE EVAL( C, phi, Nz, R, Nr, rr, k, M, Option, P )

  ! Given modes and wavenumbers, compute pressure field
  ! Normalized to pressure of point source at 1 meter
  !
  ! Option = X     Cartesian   (x, z) coordinates
  ! Option = R     Cylindrical (r, z) coordinates

  IMPLICIT NONE
  REAL,    PARAMETER :: pi = 3.1415926
  COMPLEX, PARAMETER :: i = ( 0.0, 1.0 )
  INTEGER, PARAMETER :: MaxM = 20000, MinExp = -100
  INTEGER            :: M, ir, iz, Nr, Nz
  REAL               :: rr( Nz ), r( Nr )
  COMPLEX            :: C( M ), phi( MaxM, Nz ), k( MaxM ), P( Nz, Nr ), T, Hank( M ), ik( M ), const( M ), Cmat( M, Nz ), factor
  CHARACTER (LEN=50) :: Option

  ! If no modes, return vanishing pressure
  IF ( M <= 0 ) THEN
     P = 0.0
     RETURN
  END IF

  ! Initialization
  factor = i * SQRT( 2.0 * pi ) * EXP( i * pi / 4.0 )

  IF ( Option(1:1) == 'X' ) THEN
     const( 1 : M ) = factor * C( 1 : M ) / k( 1 : M )
  ELSE
     const( 1 : M ) = factor * C( 1 : M ) / SQRT( k( 1 : M ) )
  ENDIF

  ik( 1 : M ) = -i * k( 1 : M )   ! use e{i(wt-kr)} form
  IF ( Option( 4 : 4 ) == 'I' ) ik = REAL( ik )   ! Incoherent case

  DO iz = 1, Nz
     Cmat( :, iz ) = const( : ) * phi( 1 : M, iz ) * EXP( ik( : ) * rr( iz ) )
  END DO

  ! Loop over range

  DO ir = 1, Nr
     ! eliminate underflows (can raise CPU time)
     !WHERE (  REAL( ik * r( ir ) ) > MinExp )
     Hank = EXP( ik * r( ir ) )
     !ELSEWHERE
     !   Hank = 0.0
     !END WHERE

     ! Loop over depth
     DO iz = 1, Nz
        IF ( Option( 4 : 4 ) /= 'I' )  THEN         ! coherent   case
           T =       SUM(    Cmat( :, iz ) * Hank( : ) )
        ELSE                                    ! incoherent case
           T = SQRT( SUM(  ( Cmat( :, iz ) * Hank( : ) ) ** 2 ) )
        ENDIF
        ! Cylindrical or cartesian coordinates?
        IF ( Option( 1 : 1 ) == 'R' .AND. ABS( r( ir ) + rr( iz ) ) > TINY( R( 1 ) ) ) T = T / SQRT( r( ir ) + rr( iz ) )
        P( iz, ir ) = T
     END DO

  END DO   ! next range step

END SUBROUTINE EVAL
