PROGRAM SCOOTER

  ! Finite-element wavenumber integration program

  ! Copyright (C) 2009 Michael B. Porter

  ! This program is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.

  ! You should have received a copy of the GNU General Public License
  ! along with this program.  If not, see <http://www.gnu.org/licenses/>.

  ! Initial version developed at the Naval Research Laboratory in 1985.

  USE SdRdRMod
  USE scomod
  IMPLICIT NONE
  INTEGER                       :: ik, NPoints, NTotal1, IAllocStat
  REAL                          :: kMin, kMax, deltak, TStart, TEnd
  REAL,             ALLOCATABLE :: k( : )
  REAL    (KIND=8)              :: Freq, xS, YS
  REAL    (KIND=8), ALLOCATABLE :: rho( : )
  COMPLEX (KIND=8), ALLOCATABLE :: B1( : ), B2( : ), B3( : ), B4( : )
  COMPLEX,          ALLOCATABLE :: Green( :, :, : )     ! G( Nsd, Nrd, Nk )
  CHARACTER  (LEN=80)           :: Title, FileRoot
  CHARACTER  (LEN=10)           :: PlotType

  CALL CPU_TIME( Tstart )

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  CALL GETPAR( FileRoot, Title, Freq )

  IF ( NMedia > 1 ) THEN
     IF ( ANY( sigma( 2 : NMedia ) /= 0.0 ) ) CALL ERROUT( PRTFile, 'F', 'SCOOTER', 'Rough interfaces not allowed' )
  ENDIF

  ! Set up vector of wavenumber samples
  kMin = REAL( SQRT( omega2 ) / CHigh )
  kMAX = REAL( SQRT( omega2 ) / CLow  )
  Nk   = INT( 1000.0 * RMax * ( kMax - kMin ) / REAL( pi ) )

  WRITE( PRTFile, * ) 'Nk = ', Nk
  ALLOCATE( k( Nk ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'SCOOTER', 'Insufficient memory to allocate k( Nk ) vector' )

  ! Set-up the vector of k-space points
  Deltak = ( kMax - kMin ) / ( Nk - 1 )
  Atten  = Deltak
  k      = kMin + [ ( Ik, Ik = 0, Nk - 1 ) ] * Deltak

  ALLOCATE ( Green( Nsd, Nrd, Nk ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) &
       CALL ERROUT( PRTFile, 'F', 'SCOOTER', 'Insufficient memory to allocate Green''s function matrix; reduce Rmax, Nsd, or Nrd' )

  h( 1 : NMedia ) = ( Depth( 2 : NMedia + 1 ) - Depth( 1 : NMedia ) ) / N( 1 : NMedia )  ! vector of mesh widths
  NPoints          = SUM( N( 1 : NMedia ) ) + NMedia                                  ! number of solution points

  ALLOCATE ( B1( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints ), rho( NPoints ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'SCOOTER', 'Insufficient memory to allocate B1, B2, B3, B4 vectors' )

  ! Write header for Green's function file
  PlotType   = 'Green'
  ALLOCATE( theta( 1 ) )
  theta( 1 ) = 0.0
  Ntheta     = 1   ! dummy bearing angle
  XS         = 0.0
  YS         = 0.0
  CALL WriteHeader( TRIM( FileRoot ) // '.grn', Title, theta, Ntheta, sd, Nsd, rd, Nrd, &
                    k, Nk, REAL( Freq ), REAL( Atten ), PlotType, REAL( xS ), REAL( YS ) )

  CALL INIT( B1, B2, B3, B4, rho, NPoints )    ! Initialize matrices

  NTotal1 = SUM( N( FirstAcoustic : LastAcoustic ) ) + 1       ! size of matrix for acoustic part
  CALL KERNEL( B1, B2, B3, B4, rho, NPoints, k, Green, NTotal1 )

  CALL CPU_TIME( Tend )
  WRITE( PRTFile, "(' CPU TIME: ', G15.5, 's')" ) Tend - Tstart

END PROGRAM SCOOTER
!**********************************************************************!
SUBROUTINE GetPar( FileRoot, Title, Freq )

  !     Read in the ENVFile data

  USE scomod
  USE SdRdRMod
  USE RefCoMod
  IMPLICIT NONE
  INTEGER, PARAMETER :: iProf = 1
  CHARACTER (LEN=80) :: Title, FileRoot
  REAL               :: zMin, zMAX
  REAL      (KIND=8) :: Freq

  Title = 'SCOOTER- '

  CALL READIN( FileRoot, Title, Freq, MaxMedium, NMedia, TopOpt, HSTop, N, sigma, Depth, BotOpt, HSBot, ENVFile, PRTFile )

  READ(  ENVFile, *    ) CLow, CHigh                 ! Spectral limits
  WRITE( PRTFile, "( /, ' cLow = ', G10.5, 'm/s      cHigh = ', G10.5, 'm/s' )" ) cLow, cHigh

  IF ( CLow <= 0.0 .OR. CHigh <= 0.0 .OR. CLow >= CHigh ) &
       CALL ERROUT( PRTFile, 'F', 'GETPAR', 'Need phase speeds CLow, CHigh > 0 and CLow < CHigh'  )

  READ(  ENVFile, * ) RMax                           ! Maximum range for calculations
  WRITE( PRTFile, * ) 'RMax = ', RMax
  IF ( RMax <= 0.0 ) CALL ERROUT( PRTFile, 'F', 'GETPAR', 'RMax must be positive'  )

  zMin = REAL( Depth( 1 ) )
  zMAX = REAL( Depth( NMedia + 1 ) )
  CALL SDRD( ENVFile, PRTFile, zMin, zMAX )           ! Read source/receiver depths

  CLOSE ( ENVFile )
  omega2 = ( 2.0 * pi * Freq ) ** 2

  CALL ReadRC( FileRoot, BotOpt( 1 : 1 ), TopOpt( 2 : 2 ), PRTFile )   ! Optionally read in Bot, Top reflection coefficients

END SUBROUTINE GetPar
!**********************************************************************!
SUBROUTINE INIT( B1, B2, B3, B4, rho, NPoints )

  ! Initializes arrays defining difference equations

  USE scomod
  IMPLICIT NONE
  INTEGER           :: II, J, Medium, N1, NPoints
  REAL     (KIND=8) :: rho( NPoints ), cMinV, two_h, Freq
  COMPLEX  (KIND=8) :: cp( NPoints ), cs( NPoints ), cp2, cs2
  COMPLEX  (KIND=8) :: B1( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints )
  CHARACTER (LEN=8) :: Task

  cMin          = 1.0E6
  FirstAcoustic = 0
  Loc( 1 )      = 0

  DO Medium = 1, NMedia  ! Loop over media
     IF ( Medium /= 1 ) Loc( Medium ) = Loc( Medium - 1 ) + N( Medium - 1 ) + 1
     N1   = N(   Medium ) + 1
     II   = Loc( Medium ) + 1
     Task = 'TAB'
     CALL Profil( Depth, cp( II ), cs( II ), rho( II ), Medium, N1, Freq, TopOpt( 1 : 1 ), TopOpt( 3 : 4 ), Task, ENVFile, PRTFile )

     IF ( cs( II ) == ( 0.0, 0.0 ) ) THEN   ! Case of an acoustic medium
        Material( Medium )  =  'ACOUSTIC'
        IF ( FirstAcoustic == 0 ) FirstAcoustic = Medium
        LastAcoustic = Medium

        cMinV = MINVAL( DBLE( cp( II:II + N( Medium ) ) ) )
        cMin  = MIN( cMin, cMinV )
        B1( II : II + N( Medium ) ) = omega2 / cp( II:II + N( Medium ) ) ** 2

     ELSE                                  ! Case of an elastic medium
        Material( Medium ) = 'ELASTIC'
        two_h         = 2.0 * h( Medium )

        DO j = II, II + N( Medium )
           cMin = MIN( DBLE( cs( j ) ), cMin )

           cp2 = cp( j ) ** 2
           cs2 = cs( j ) ** 2

           B1(  j ) = two_h / ( rho( j ) * cs2 )
           B2(  j ) = two_h / ( rho( j ) * cp2 )
           B3(  j ) = 4.0 * two_h * rho( j ) * cs2 * ( cp2 - cs2 ) / cp2
           B4(  j ) = two_h * ( cp2 - 2.0 * cs2 ) / cp2
           rho( j ) = two_h * omega2 * rho( j )
        END DO

     ENDIF
  END DO   ! next Medium

END SUBROUTINE INIT
!**********************************************************************!
SUBROUTINE BCimpedance( B1, B2, B3, B4, rho, NPoints, x, BotTop, HS, F, G, IPower )

  !     Compute Boundary Condition Impedance
  !     Same subroutine as in KRAKENC except
  !        PEKRT    is replaced by SQRT
  !        COMC     is replaced by COMSCO
  !        cInside  is related to B1 differently

  USE scomod
  USE RefCoMod
  IMPLICIT NONE
  COMPLEX (KIND=8), PARAMETER :: zero = (0.0D0, 0.0D0 )
  INTEGER, INTENT( OUT ) :: IPower
  INTEGER                :: Ibot, Itop, Medium, NPoints
  REAL     (KIND=8)      :: rho( NPoints ), c0, RadDeg, rhoInside, omega
  COMPLEX  (KIND=8)      :: x, kx, kz, Twersky, gammaS2, gammaP2, gammaS, gammaP, mu, yV( 5 ), RCmplx, &
                            B1( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints ), cInside
  COMPLEX  (KIND=8), INTENT( OUT ) :: F, G
  CHARACTER (LEN=3)      :: BotTop
  TYPE( HSInfo )         :: HS
  TYPE( ReflectionCoef ) :: RInt

  IPower = 0

  ! Get rho, C just INSide the boundary
  SELECT CASE ( BotTop )
  CASE ( 'TOP' )
     Itop      = 1
     rhoInside = rho( Itop )
     cInside   = SQRT( omega2 / B1( Itop ) )
  CASE ( 'BOT' )
     Ibot      = Loc( LastAcoustic ) + N( LastAcoustic ) + 1
     rhoInside = rho( Ibot )
     cInside   = SQRT( omega2 / B1( Ibot ) )
  END SELECT

  ! impedance for different bottom types

  SELECT CASE ( HS%BC )
  CASE ( 'V' )                   ! Vacuum
     F       = 1.0
     G       = -i * SQRT( omega2 / cInside ** 2 - x ) * sigma( 1 ) ** 2
     yV      = CMPLX( [ F, G, zero, zero, zero ] )
  CASE (  'S', 'H', 'T', 'I'  )  ! Vacuum with Twersky scatter model
     omega   = SQRT( omega2 )
     kx      = SQRT( x )
     F       = 1.0
     c0      = REAL( cInside )
     G       = Twersky( omega, HS, kx, rhoInside, c0 )
     G       = G / ( i * omega * rhoInside )
     yV      = CMPLX( [ F, G, zero, zero, zero ] )
  CASE ( 'R' )                    ! Rigid
     F       = 0.0
     G       = 1.0
     yV      = CMPLX( [ F, G, zero, zero, zero ] )
  CASE ( 'A' )                    ! Acousto-elastic half-space
     IF ( REAL( HS%cs ) > 0.0 ) THEN
        gammaS2 = x - omega2 / HS%cs ** 2
        gammaP2 = x - omega2 / HS%cp ** 2
        gammaS  = SQRT( gammaS2 )
        gammaP  = SQRT( gammaP2 )
        mu      = HS%rho * HS%cs ** 2

        yV( 1 ) = ( gammaS * gammaP - x ) / mu
        yV( 2 ) = ( ( gammaS2 + x ) ** 2 - 4.0 * gammaS * gammaP * x ) * mu
        yV( 3 ) = 2.0 * gammaS * gammaP - gammaS2 - x
        yV( 4 ) = gammaP * ( x - gammaS2 )
        yV( 5 ) = gammaS * ( gammaS2 - x )

        F = omega2 * yV( 4 )
        G = yV( 2 )
     ELSE
        gammaP = SQRT( x - omega2 / HS%cp ** 2 )
        F    = 1.0
        G    = HS%rho / gammaP
     ENDIF
  CASE ( 'F' )                    ! Tabulated reflection coefficient
     ! Compute the grazing angle Theta
     kx         = SQRT( x )
     kz         = SQRT( omega2 / cInside ** 2 - kx ** 2 )
     RadDeg     = 180.0 / pi
     RInt%theta = RadDeg * ATAN2( REAL( kz ), REAL( kx ) )

     ! Evaluate R( TheInt )
     SELECT CASE ( BotTop )
     CASE ( 'TOP' )
        CALL RefCo( RInt, RTop, NTopPts, PRTFile )
     CASE ( 'BOT' )
        CALL RefCo( RInt, RBot, NBotPts, PRTFile )
     END SELECT

     ! Convert R( Theta ) to (f,g) in Robin BC
     RCmplx = RInt%R * EXP( i * RInt%phi )
     F      = 1.0
     G      = ( 1.0 + RCmplx ) / ( i * kz * ( 1.0 - RCmplx ) )
  CASE ( 'P' )                    ! Precalculated reflection coef
     CALL IRCINT( x, F, G, IPower, xTab, FTab, GTab, ITab, NkTab )
  END SELECT

  IF ( BotTop == 'TOP' ) G = -G    ! A top BC has the sign flipped relative to a bottom BC

  ! Shoot through elastic layers
  SELECT CASE ( BotTop )
  CASE ( 'TOP' )
     IF ( FirstAcoustic > 1 ) THEN
        DO Medium = 1, FirstAcoustic - 1            ! Shooting down from top
           CALL ElasDn( B1, B2, B3, B4, rho, NPoints, x, yV, IPower, Medium )
        END DO
        F = omega2 * yV( 4 )
        G = yV( 2 )
     ENDIF
  CASE ( 'BOT' )
     IF ( LastAcoustic < NMedia ) THEN
        DO Medium = NMedia, LastAcoustic + 1, -1    ! Shooting up from bottom
           CALL ElasUp( B1, B2, B3, B4, rho, NPoints, x, yV, IPower, Medium )
        END DO
        F = omega2 * yV( 4 )
        G = yV( 2 )
     ENDIF
  END SELECT

END SUBROUTINE Bcimpedance
!**********************************************************************!
SUBROUTINE ElasUp( B1, B2, B3, B4, rho, NPoints, x, yV, IPower, Medium )

  ! Propagates through an elastic layer using compound matrix formulation

  USE scomod
  IMPLICIT NONE
  INTEGER          :: IPowerR, IPowerF, II, J, NPoints, IPower, Medium
  REAL    (KIND=8) :: rho( NPoints ), two_h, Roof, Floor
  COMPLEX (KIND=8) :: X, xV( 5 ), yV( 5 ), zV( 5 ), two_x, xB3, four_h_x
  COMPLEX (KIND=8) :: B1( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints )
  PARAMETER ( Roof = 1.0E5, Floor = 1.0E-5, IPowerR = 5, IPowerF = -5 )

  ! Euler's method for first step
  two_x    = 2.0 * x
  two_h    = 2.0 * h( Medium )
  four_h_x = 4.0 * h( Medium ) * x
  j      = Loc( Medium ) + N( Medium ) + 1

  xB3 = x * B3( j ) - rho( j )

  zV( 1 ) = yV( 1 ) - 0.5 * (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
  zV( 2 ) = yV( 2 ) - 0.5 * ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
  zV( 3 ) = yV( 3 ) - 0.5 * (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
  zV( 4 ) = yV( 4 ) - 0.5 * (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
  zV( 5 ) = yV( 5 ) - 0.5 * (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

  ! Modified midpoint method
  DO II = N( Medium ), 1, -1
     j = j - 1

     xV = yV
     yV = zV

     xB3 = x * B3( j ) - rho( j )

     zV( 1 ) = xV( 1 ) - (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
     zV( 2 ) = xV( 2 ) - ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
     zV( 3 ) = xV( 3 ) - (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
     zV( 4 ) = xV( 4 ) - (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
     zV( 5 ) = xV( 5 ) - (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

     ! Scale if necessary
     IF ( II /= 1 ) THEN
        IF      ( ABS( DBLE( zV( 2 ) ) ) < Floor ) THEN
           zV     = Roof * zV
           yV     = Roof * yV
           IPower = IPower - IPowerR
        ELSE IF ( ABS( DBLE( zV( 2 ) ) ) > Roof ) THEN
           zV     = Floor * zV
           yV     = Floor * yV
           IPower = IPower - IPowerF
        ENDIF

     ENDIF
  END DO

  yV = ( xV + 2.0 * yV + zV ) / 4.0   ! Apply the standard filter at the terminal point

END SUBROUTINE ElasUp
!**********************************************************************!
SUBROUTINE ElasDn( B1, B2, B3, B4, rho, NPoints, x, yV, IPower, Medium )

  ! Propagates through an elastic layer using compound matrix formulation

  USE scomod
  IMPLICIT NONE
  INTEGER          :: IPowerR, IPowerF, II, J, NPoints, IPower, Medium
  REAL    (KIND=8) :: rho( NPoints ), two_h, Roof, Floor
  COMPLEX (KIND=8) :: x, xV( 5 ), yV( 5 ), zV( 5 ), two_x, xB3, four_h_x
  COMPLEX (KIND=8) :: B1( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints )
  PARAMETER ( Roof = 1.0E5, Floor = 1.0E-5, IPowerR = 5, IPowerF = -5 )

  ! Euler's method for first step

  two_x    = 2.0 * x
  two_h    = 2.0 * h( Medium )
  four_h_x = 4.0 * h( Medium ) * x
  j        = Loc( Medium ) + 1

  xB3 = x * B3( j ) - rho( j )

  zV( 1 ) = yV( 1 ) + 0.5 * (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
  zV( 2 ) = yV( 2 ) + 0.5 * ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
  zV( 3 ) = yV( 3 ) + 0.5 * (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
  zV( 4 ) = yV( 4 ) + 0.5 * (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
  zV( 5 ) = yV( 5 ) + 0.5 * (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

  ! Modified midpoint method
  DO II = 1, N( Medium )
     j = j + 1

     xV = yV
     yV = zV

     xB3 = x * B3( j ) - rho( j )

     zV( 1 ) = xV( 1 ) + (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
     zV( 2 ) = xV( 2 ) + ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
     zV( 3 ) = xV( 3 ) + (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
     zV( 4 ) = xV( 4 ) + (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
     zV( 5 ) = xV( 5 ) + (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

     ! Scale if necessary
     IF ( II /= N( Medium ) ) THEN
        IF     ( ABS( DBLE( zV( 2 ) ) ) < Floor ) THEN
           zV     = Roof * zV
           yV     = Roof * yV
           IPower = IPower - IPowerR
        ELSE IF ( ABS( DBLE( zV( 2 ) ) ) > Roof ) THEN
           zV     = Floor * zV
           yV     = Floor * yV
           IPower = IPower - IPowerF
        ENDIF
     ENDIF
  END DO

  yV = ( xV + 2.0 * yV + zV ) / 4.0   ! Apply the standard filter at the terminal point


END SUBROUTINE ElasDn
!**********************************************************************!
SUBROUTINE Kernel( B1, B2, B3, B4, rho, NPoints, k, Green, NTotal1 )

  ! Solve system for a sequence of k-values

  USE SdRdRMod
  USE scomod
  IMPLICIT NONE
  INTEGER          :: II, Ik, Is, Ir, j, l, Medium, NPoints, NTotal1
  REAL             :: z( NTotal1 ), rhoElement( NPoints ), k( Nk )
  REAL    (KIND=8) :: rho( NPoints ), rhoh
  COMPLEX          :: Green( Nsd, Nrd, Nk )
  COMPLEX (KIND=8) :: BElement, DF( NTotal1 ), EF( NTotal1 )
  COMPLEX (KIND=8) :: B1( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints ), x

  ! Tabulate z coordinates
  z( 1 ) = REAL( Depth( FirstAcoustic ) )
  j      = 2

  DO Medium = FirstAcoustic, LastAcoustic
     z( j : j + N( Medium ) - 1 ) = REAL( Depth( Medium ) + [ ( II * h( Medium ), II = 1, N( Medium ) ) ] )
     j = j + N( Medium )
  END DO

  ! Compute weights for source/rcvr depth interpolation
  CALL WEIGHT( z, NTotal1, sd, Nsd, WS, Isd )
  CALL WEIGHT( z, NTotal1, rd, Nrd, WR, Ird )

  ! Assemble matrix
  j       = 1
  l       = Loc( FirstAcoustic ) + 1
  DF( 1 ) = 0.0

  DO Medium = FirstAcoustic, LastAcoustic
     DO II = 1, N( Medium )
        rhoElement( l ) = REAL( rho( l ) + rho( l + 1 ) ) / 2.0
        rhoH            = rhoElement( l ) * h( Medium )
        BElement        = h( Medium ) * ( ( B1( l ) + B1( l + 1 ) ) / 2.0 ) / ( 12.0 * rhoElement( l ) )

        DF( j     ) = DF( j ) - 1.0 / rhoH + 5.0 * BElement
        DF( j + 1 ) =         - 1.0 / rhoH + 5.0 * BElement
        EF( j + 1 ) =           1.0 / rhoH +       BElement

        j = j + 1
        l = l + 1
     END DO
     l = l + 1
  END DO

  DO Ik = 1, Nk   ! Step through each point in k-space
     x = ( k( Ik ) + i * Atten ) ** 2
     CALL Solve( B1, B2, B3, B4, rho, NPoints, NTotal1, x, Green, Ik, DF, EF, rhoElement )  ! Solve for G(k)
  END DO

  ! Write Green's function to file
  DO IS = 1, Nsd
     DO IR = 1, Nrd
        WRITE( GRNFile, REC = 7 + ( IS - 1 ) * Nrd + IR ) Green( IS, IR, : )
     END DO
  END DO

  CLOSE( GRNFile )

END SUBROUTINE Kernel

!**********************************************************************!

SUBROUTINE Solve( B1, B2, B3, B4, rho, NPoints, NTotal1, x, Green, Ik, DF, EF, rhoElement )

  ! Set up the linear system and solve

  USE SdRdRMod
  USE scomod
  IMPLICIT NONE
  INTEGER           :: NPoints, j, l, Medium, II, ik, Is, IPower, NTotal1
  REAL              :: rhoElement( * )
  REAL     (KIND=8) :: rho( NPoints ), rhoSd
  COMPLEX           :: Green( Nsd, Nrd, Nk )
  COMPLEX  (KIND=8) :: d( NTotal1 ), e( NTotal1 ), RV1( NTotal1 ), RV2( NTotal1 ), RV4( NTotal1 )
  COMPLEX  (KIND=8) :: DF( * ), EF( * ), BElement, xT, x, F, G
  COMPLEX  (KIND=8) :: B1( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints )

  ! Complete assembly of matrix by adding in x
  j = 1
  l = Loc( FirstAcoustic ) + 1
  d( 1 ) = DF( 1 )

  DO Medium = FirstAcoustic, LastAcoustic
     xT = -h( Medium ) * x / 12.0

     DO II = 1, N( Medium )
        BElement = xT / rhoElement( l )
        d( j     ) = d(  j     ) + 5.0 * BElement
        d( j + 1 ) = DF( j + 1 ) + 5.0 * BElement
        e( j + 1 ) = EF( j + 1 ) +       BElement
        j = j + 1
        l = l + 1
     END DO

     l = l + 1
  END DO

  ! Corner elt requires top impedance

  CALL BCImpedance( B1, B2, B3, B4, rho, NPoints, x, 'TOP', HSTop, F, G, IPower )
  IF ( G == 0.0 ) THEN
     d( 1 ) = 1.0D0
     e( 2 ) = 0.0D0
  ELSE
     d( 1 ) = d( 1 ) + F / G
  ENDIF

  ! Corner elt requires bottom impedance

  CALL BCImpedance( B1, B2, B3, B4, rho, NPoints, x, 'BOT', HSBot, F, G, IPower )
  IF ( G == 0.0 ) THEN
     d( NTotal1 ) = 1.0D0
     e( NTotal1 ) = 0.0D0
  ELSE
     d( NTotal1 ) =  d( NTotal1 ) - F / G
  ENDIF

  CALL FACTOR( NTotal1, d, e, RV1, RV2, RV4 )   !     * Do LU decomposition *

  DO IS = 1, Nsd   ! Loop over all source positions

     ! Set up RHS in D (previously used for diagonal)
     rhosd       = 1.0    ! assumes rho( zs ) = 1
     D           = 0.0
     II          = Isd( IS )

     d( II     ) = 2.0 * ( 1.0 - WS( IS ) ) / rhosd
     d( II + 1 ) = 2.0 *     WS( IS )       / rhosd

     CALL BackSb( NTotal1, RV1, RV2, RV4, d )    ! Solve the system
     Green( IS, :, Ik ) = CMPLX( d( Ird ) + WR * ( d( Ird + 1 ) - d( Ird ) ) )   ! extract the solution at the rcvr depths
  END DO

END SUBROUTINE Solve
