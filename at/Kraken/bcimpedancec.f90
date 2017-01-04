SUBROUTINE BCImpedance( x, BotTop, HS, F, G, IPower )

  ! Compute Boundary Condition Impedance

  USE krakcmod
  USE RefCoMod

  IMPLICIT NONE
  COMPLEX (KIND=8), PARAMETER :: zero = (0.0D0, 0.0D0 )
  INTEGER                :: Itop, Ibot, Medium
  INTEGER, INTENT( OUT ) :: IPower
  REAL          (KIND=8) :: omega, c0, rhoInside, RadDeg
  COMPLEX       (KIND=8) :: x, Kx, kz, Twersky, gammaS2, gammaP2, gammaS, gammaP, mu, &
                            yV( 5 ), PekerisRoot, RCmplx, CInside
  COMPLEX (KIND=8), INTENT( OUT ) :: F, G
  CHARACTER      (LEN=3) :: BotTop
  TYPE( HSInfo )         :: HS
  TYPE(ReflectionCoef)   :: RInt

  IPower = 0

  ! Get rho, c just INSide the boundary
  ! (There is at least one acoustic layer in the problem, except
  !  in the case where BOUNCE is used to get the refl. coef. for
  !  a purely elastic stack of layers.

  SELECT CASE ( BotTop )
  CASE ( 'TOP' )
     IF ( FirstAcoustic > 0 ) THEN
        Itop       = Loc( FirstAcoustic ) + N( FirstAcoustic ) + 1
        rhoInside  = rho( Itop )
        CInside    = PekerisRoot( omega2 * h( FirstAcoustic ) **2 / ( 2.0D0 + B1( Itop ) ) )
     ENDIF
  CASE ( 'BOT' )
     IF ( LastAcoustic > 0 ) THEN
        Ibot       = Loc( LastAcoustic ) + N( LastAcoustic ) + 1
        rhoInside  = rho( Ibot )
        CInside    = PekerisRoot( omega2 * h( LastAcoustic ) **2 / ( 2.0D0 + B1( Ibot ) ) )
     ENDIF
  END SELECT

  SELECT CASE ( HS%BC )
  CASE ( 'V' )                   ! Vacuum
     F       = 1.0D0
     !G = -i * PekerisRoot( omega2 / CInside ** 2 - x ) * SIGMA( 1 ) ** 2
     G       = 0.0D0
     yV      = CMPLX( [ F, G, zero, zero, zero ] )
  CASE (  'S', 'H', 'T', 'I'  )  ! Vacuum with Twersky scatter model
     omega   = SQRT( omega2 )
     kx      = SQRT( x )
     F       = 1.0D0
     C0      = REAL( CInside )
     G       = Twersky( omega, HS, kx, rhoInside, C0 )
     G       = G / ( i * omega * rhoInside )
     yV      = CMPLX( [ F, G, zero, zero, zero ] )
  CASE ( 'R' )                    ! Rigid
     F       = 0.0D0
     G       = 1.0D0
     yV      = CMPLX( [ F, G, zero, zero, zero ] )
  CASE ( 'A' )                    ! Acousto-elastic half-space
     IF ( REAL( HS%CS ) > 0.0 ) THEN
        gammaS2 = x - omega2 / HS%CS ** 2
        gammaP2 = x - omega2 / HS%CP ** 2
        gammaS  = PekerisRoot( gammaS2 )
        gammaP  = PekerisRoot( gammaP2 )

        mu = HS%rho * HS%CS ** 2

        yV( 1 ) = ( gammaS * gammaP - x ) / mu
        yV( 2 ) = ( ( gammaS2 + x ) ** 2 - 4.0D0 * gammaS * gammaP * x ) * mu
        yV( 3 ) = 2.0D0 * gammaS * gammaP - gammaS2 - x
        yV( 4 ) = gammaP * ( x - gammaS2 )
        yV( 5 ) = gammaS * ( gammaS2 - x )

        F = omega2 * yV( 4 )
        G = yV( 2 )
     ELSE
        gammaP = PekerisRoot( x - omega2 / HS%CP ** 2 )
        F = 1.0D0
        G = HS%rho / gammaP
     ENDIF
  CASE ( 'F' )                    ! Tabulated reflection coefficient
     ! Compute the grazing angle THETA
     kx       = SQRT( x )
     kz       = SQRT( omega2 / CInside ** 2 - x )
     RadDeg   = 180.0D0 / pi
     RInt%theta = RadDeg * DATAN2( DBLE( kz ), DBLE( kx ) )

     ! Evaluate R( ThetaInt )
     IF ( BotTop == 'TOP' ) THEN
        CALL REFCO( RInt, RTop, NTopPts, PRTFile )
     ELSE
        CALL REFCO( RInt, RBot, NBotPts, PRTFile )
     ENDIF

     ! Convert R(THETA) to (f,g) in Robin BC
     RCmplx = RInt%R * EXP( i * RInt%phi )
     F      = 1.0D0
     G      = ( 1.0D0 + RCmplx ) / ( i * kz * ( 1.0D0 - RCmplx ) )

  CASE ( 'P' )                    ! Precalculated reflection coef
     CALL IRCINT( x, F, G, IPower, xTab, FTab, GTab, ITab, NkTab )
  END SELECT

  IF ( BotTop == 'TOP' ) G = -G    ! A top BC has the sign flipped relative to a bottom BC

  ! Shoot through elastic layers
  SELECT CASE ( BotTop )
  CASE ( 'TOP' )
     IF ( FirstAcoustic > 1 ) THEN   ! Shoot down from top

        DO Medium = 1, FirstAcoustic - 1
           CALL ELASDN( x, yV, IPower, Medium )
        END DO

        F = omega2 * yV( 4 )
        G = yV( 2 )
     ENDIF
  CASE ( 'BOT' )
     IF ( LastAcoustic < NMedia ) THEN   ! Shoot up from bottom

        DO Medium = NMedia, LastAcoustic + 1, -1
           CALL ELASUP( x, yV, IPower, Medium )
        END DO

        F = omega2 * yV( 4 )
        G = yV( 2 )
     ENDIF
  END SELECT

END SUBROUTINE BCImpedance
!**********************************************************************C
SUBROUTINE ELASUP( x, yV, IPower, Medium )

  ! Propagates through an elastic layer using compound matrix formulation

  USE krakcmod

  IMPLICIT NONE
  INTEGER,       PARAMETER :: IPowerR = 50, IPowerF = -50
  REAL (KIND=8), PARAMETER :: Roof = 1.0D+50, Floor = 1.0D-50
  INTEGER                  :: IPower, Medium, II, j
  REAL    (KIND=8)         :: two_h
  COMPLEX (KIND=8)         :: x, xV( 5 ), yV( 5 ), zV( 5 ), two_x, xB3, four_h_x


  ! Euler's method for first step

  two_x    = 2.0D0 * x
  two_h    = 2.0D0 * h( Medium )
  four_h_x = 4.0D0 * h( Medium ) * x
  j        = Loc( Medium ) + N( Medium ) + 1
  xB3      = x * B3( j ) - rho( j )

  zV( 1 ) = yV( 1 ) - 0.5D0 * (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
  zV( 2 ) = yV( 2 ) - 0.5D0 * ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
  zV( 3 ) = yV( 3 ) - 0.5D0 * (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
  zV( 4 ) = yV( 4 ) - 0.5D0 * (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
  zV( 5 ) = yV( 5 ) - 0.5D0 * (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

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
        IF ( ABS( DBLE( zV( 2 ) ) ) < Floor ) THEN
           zV     = Roof * zV
           yV     = Roof * yV
           IPower = IPower - IPowerR
        ENDIF

        IF ( ABS( DBLE( zV( 2 ) ) ) > Roof  ) THEN
           zV     = Floor * zV
           yV     = Floor * yV
           IPower = IPower - IPowerF
        ENDIF
     ENDIF
  END DO

  yV = ( xV + 2.0D0 * yV + zV ) / 4.0D0   ! Apply the standard filter at the terminal point

END SUBROUTINE ELASUP
!**********************************************************************C
SUBROUTINE ELASDN( x, yV, IPower, Medium )

  ! Propagates through an elastic layer using compound matrix formulation

  USE krakcmod

  IMPLICIT NONE
  INTEGER,       PARAMETER :: IPowerR = 50, IPowerF = -50
  REAL (KIND=8), PARAMETER :: Roof = 1.0D+50, Floor = 1.0D-50
  INTEGER                  :: IPower, Medium, II, j
  REAL            (KIND=8) :: two_h
  COMPLEX         (KIND=8) :: x, xV( 5 ), yV( 5 ), zV( 5 ), two_x, xB3, four_h_x

  ! Euler's method for first step

  two_x    = 2.0D0 * x
  two_h    = 2.0D0 * h( Medium )
  four_h_x = 4.0D0 * h( Medium ) * x
  j        = Loc( Medium ) + 1
  xB3      = x * B3( j ) - rho( j )

  zV( 1 ) = yV( 1 ) + 0.5D0 * (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
  zV( 2 ) = yV( 2 ) + 0.5D0 * ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
  zV( 3 ) = yV( 3 ) + 0.5D0 * (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
  zV( 4 ) = yV( 4 ) + 0.5D0 * (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
  zV( 5 ) = yV( 5 ) + 0.5D0 * (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

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
        IF ( ABS( DBLE( zV( 2 ) ) ) < Floor ) THEN
           zV     = Roof * zV
           yV     = Roof * yV
           IPower = IPower - IPowerR
        ENDIF

        IF ( ABS( DBLE( zV( 2 ) ) ) > Roof  ) THEN
           zV     = Floor * zV
           yV     = Floor * yV
           IPower = IPower - IPowerF
        ENDIF
     ENDIF
  END DO

  yV = ( xV + 2.0D0 * yV + zV ) / 4.0D0   ! Apply the standard filter at the terminal point

END SUBROUTINE ELASDN
