SUBROUTINE EVALAD( FileRoot, rProf, NProf, phiS, phiR, rd, Nrd, r, Nr, M, Opt, P )

  ! Computes pressure field using adiabatic mode theory
  ! Normalized to pressure at 1 meter

  ! Opt:
  !   X     Cartesian   (x, z) coordinates
  !   T     Translationally invariant ocean
  !   R     Cylindrical (r, z) coordinates
  !   S     Scaled cylindrical coordinates ( 1/r fall-off removed )

  IMPLICIT NONE
  INTEGER, PARAMETER :: MaxM = 20000
  COMPLEX, PARAMETER :: i = ( 0.0, 1.0 )
  REAL,    PARAMETER :: pi = 3.1415926
  INTEGER            :: iProf, ir, ird, NProf, Nrd, Nr, M, M1
  REAL               :: rd( * ), rProf( * ), r( * ), freq, rLeft, rMid, w, wMid
  COMPLEX            :: phiS( * ), phiL( MaxM, Nrd ), phiR( MaxM, * ), const( MaxM ), Hank( MaxM ), &
                        SUM, kInt( MaxM ), phiINT( MaxM ), P( Nrd, * ), kL( MaxM ), kR( MaxM ), kMid( MaxM )
  COMPLEX   (KIND=8) :: sumk( MaxM ), sumkinv( MaxM )
  CHARACTER          :: Opt*( * )
  CHARACTER (LEN=80) :: Title, FileRoot

  ! Initialization
  rProf( NProf + 1 ) = HUGE( rProf( NProf ) )
  iProf              = 1

  ! Receiver depths at left  of segment
  CALL GetModes( FileRoot, iProf    , MaxM, rd, Nrd, 'N', kL, phiL, M1, freq, Title )
  M = MIN( M, M1 )

  ! Receiver depths at right of segment
  CALL GetModes( FileRoot, iProf + 1, MaxM, rd, Nrd, 'N', kR, phiR, M1, freq, Title )
  M = MIN( M, M1 )

  const(   1 : M ) = i * SQRT( 2.0 * pi ) * EXP( i * pi / 4.0 ) * phiS( 1 : M )
  sumk(    1 : M ) = 0.0
  sumkinv( 1 : M ) = 0.0
  IF ( Opt( 1 : 1 ) == 'T' ) const( 1 : M ) = const( 1 : M ) / SQRT( kL( 1 : M ) )

  ! March forward in range
  DO ir = 1, Nr
     ! Crossing into new range segment?
     DO WHILE ( r( ir ) > 1000.0 * rProf( iProf + 1 ) )
        CALL NewSegment( FileRoot, r, ir, rProf, iProf, NProf, kL, kR, phiL, phiR, M, MaxM, rd, Nrd, sumk, sumkinv )
     END DO

     ! Compute proportional distance, W, and interpolate
     IF ( ir > 1 ) THEN
        rLeft = MAX( r( ir - 1 ), 1000.0 * rProf( iProf ) )
     ELSE
        rLeft = 1000.0 * rProf( iProf )
     ENDIF

     rMid = 0.5 * ( r( ir ) + rLeft )
     W    = ( r( ir ) / 1000.0 - rProf( iProf ) ) / ( rProf( iProf + 1 ) - rProf( iProf ) )
     wMid = ( rMid    / 1000.0 - rProf( iProf ) ) / ( rProf( iProf + 1 ) - rProf( iProf ) )

     kInt(    1 : M ) = kL(      1 : M ) + W    * ( kR( 1 : M ) - kL( 1 : M ) )
     kMid(    1 : M ) = kL(      1 : M ) + wMid * ( kR( 1 : M ) - kL( 1 : M ) )
     sumk(    1 : M ) = sumk(    1 : M ) + kMid(1 : M) * ( r( ir ) - rLeft )
     sumkinv( 1 : M ) = sumkinv( 1 : M ) + ( r( ir ) - rLeft ) / kMid( 1 : M )

     IF ( Opt( 4 : 4 ) /= 'I' ) THEN   ! coherent   case
        Hank( 1 : M ) = const( 1 : M ) * CMPLX( EXP(       -i * sumk( 1 : M ) ) )
     ELSE                              ! incoherent case
        Hank( 1 : M ) = const( 1 : M ) * CMPLX( EXP( REAL( -i * sumk( 1 : M ) ) ) )
     ENDIF

     SELECT CASE( Opt( 1 : 1 ) )
     CASE ( 'R' )                             ! Cylindrical coords.
        IF ( r( ir ) == 0.0 ) THEN
           Hank( 1 : M ) = 0.0
        ELSE
           Hank( 1 : M ) = Hank( 1 : M ) / SQRT( kInt( 1 : M ) * r( ir ) )
        ENDIF
     CASE ( 'X' )                             ! Cartesian coords.
        Hank( 1 : M ) = Hank( 1 : M ) / kInt (1 : M )
     CASE ( 'T' )                             ! Translationally invariant
        Hank( 1 : M ) = Hank( 1 : M ) / CMPLX( SQRT( kInt( 1 : M ) * sumkinv( 1 : M ) ) )
     CASE DEFAULT                             ! Scaled cylindrical coords.
        Hank( 1 : M ) = Hank( 1 : M ) / SQRT( kInt( 1 : M ) )
     END SELECT

     ! For each rcvr, add up modal contributions
     DO ird = 1, Nrd
        phiINT( 1 : M ) = phiL( 1 : M, ird ) + W * ( phiR( 1 : M, ird ) - phiL( 1 : M, ird ) )
        IF ( Opt( 4 : 4 ) /= 'I' )  THEN   ! coherent   case
           P( ird, ir ) =       SUM(       phiINT( 1 : M ) * Hank( 1 : M ) )
        ELSE                               ! incoherent case
           P( ird, ir ) = SQRT( SUM( ABS(  phiINT( 1 : M ) * Hank( 1 : M ) ) ** 2 ) )
        ENDIF
     END DO

  END DO   ! Next range step

END SUBROUTINE EVALAD

!**********************************************************************!

SUBROUTINE NewSegment( FileRoot, r, ir, rProf, iProf, NProf, kL, kR, phiL, phiR, M, MaxM, rd, Nrd, sumk, sumkinv )

  ! Treats the crossing into a new range segment

  IMPLICIT NONE
  INTEGER            :: ir, iProf, M, M1, MaxM, NProf, Nrd
  REAL               :: rd( * ), rProf( * ), r( * ), freq, rLeft, rMid, wMid
  COMPLEX            :: phiL( MaxM, * ), phiR( MaxM, * ), kL( * ), kR( * ), kMid( M )
  COMPLEX   (KIND=8) :: sumk( * ), sumkinv( * )
  CHARACTER (LEN=80) :: Title, FileRoot

  ! Do phase integral up to the new range
  IF ( ir > 1 ) THEN
     rLeft = MAX( r( ir - 1 ), 1000.0 * rProf( iProf ) )
  ELSE
     rLeft = 1000.0 * rProf( iProf )
  ENDIF

  rMid             = 0.5 * ( 1000.0 * rProf( iProf + 1 ) + rLeft )
  wMid             = ( rMid / 1000.0 - rProf( iProf ) ) / ( rProf( iProf + 1 ) - rProf( iProf ) )
  kMid(    1 : M ) = kL(      1 : M ) + wMid  * ( kR( 1 : M ) - kL( 1 : M ) )
  sumk(    1 : M ) = sumk(    1 : M ) + kMid( 1 : M ) * ( 1000.0 * rProf( iProf + 1 ) - rLeft )
  sumkinv( 1 : M ) = sumkinv( 1 : M ) + ( 1000.0 * rProf( iProf + 1 ) - rLeft ) / kMid( 1 : M )

  ! Copy right modes to left
  kL(   1 : M )        = kR(   1 : M )
  phiL( 1 : M, 1 : Nrd ) = phiR( 1 : M, 1 : Nrd ) 

  ! Read in the new right mode set
  iProf = iProf + 1

  IF ( iProf + 1  <= NProf ) THEN
     CALL GetModes( FileRoot, iProf + 1, MaxM, rd, Nrd, 'N', kR, phiR, M1, freq, Title )
     M = MIN( M, M1 )
     WRITE( *, * ) 'New profile read', r( ir ), iProf + 1, rProf( iProf + 1 ), ' #modes=', M
  ENDIF

END SUBROUTINE NewSegment
