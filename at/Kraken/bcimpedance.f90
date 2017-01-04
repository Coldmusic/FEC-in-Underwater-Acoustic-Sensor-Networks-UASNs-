SUBROUTINE BCImpedance( x, BotTop, HS, F, G, IPower )

  ! Compute Boundary Condition Impedance

  USE krakmod

  IMPLICIT NONE
  INTEGER, INTENT( OUT ) :: IPower
  INTEGER                :: Medium
  REAL (KIND=8)          :: yV( 5 ), x, gammaS, gammaP, mu
  REAL (KIND=8), INTENT( OUT ) :: F, G
  COMPLEX (KIND=8)       :: gammaS2, gammaP2
  CHARACTER      (LEN=3) :: BotTop
  TYPE( HSInfo )         :: HS

  IPower = 0

  SELECT CASE( HS%BC )
  CASE( 'V', 'S', 'H', 'T', 'I' )   ! Vacuum or Twersky
     F  = 1.0D0
     G  = 0.0D0
     yV = [ F, G, 0.D0, 0.D0, 0.D0 ]
  CASE ( 'R' )                      ! Rigid
     F  = 0.0D0
     G  = 1.0D0
     yV = [ F, G, 0.D0, 0.D0, 0.D0 ]
  CASE ( 'A' )                      ! Acousto-elastic half-space
     IF ( REAL( HS%cS ) > 0.0 ) THEN
        gammaS2 = x - omega2 / DBLE( HS%cS ) ** 2
        gammaP2 = x - omega2 / DBLE( HS%cP ) ** 2
        gammaS  = DBLE( SQRT( gammaS2 ) )
        gammaP  = DBLE( SQRT( gammaP2 ) )
        mu      = HS%rho * DBLE( HS%cS ) ** 2

        yV( 1 ) = ( gammaS * gammaP - x ) / mu
        yV( 2 ) = ( ( gammaS2 + x ) ** 2 - 4.0D0 * gammaS * gammaP * x ) * mu
        yV( 3 ) = 2.0D0 * gammaS * gammaP - gammaS2 - x
        yV( 4 ) = gammaP * ( x - gammaS2 )
        yV( 5 ) = gammaS * ( gammaS2 - x )

        F = omega2 * yV( 4 )
        G = yV( 2 )
        IF ( G > 0.0D0 ) ModeCount = ModeCount + 1
     ELSE
        gammaP2 = x - omega2 / DBLE( HS%cP ) ** 2
        gammaP  = DBLE( SQRT( x - omega2 / HS%cP ** 2 ) )
        F       = 1.0D0

        IF ( gammaP /= 0.0D0 ) THEN
           G    = HS%rho / gammaP
        ELSE
           G    = 0.0D0
        END IF
     ENDIF
  CASE ( 'F' )                      ! Tabulated reflection coefficient
     CALL ERROUT( PrtFile, 'F', 'KRAKEN - BCIMP', 'KRAKEN does not allow a tabulated reflection coef.' )
  CASE ( 'P' )                      ! Precalculated reflection coef
     CALL ERROUT( PrtFile, 'F', 'KRAKEN - BCIMP', 'KRAKEN does not allow a precalculated reflection coef.' )
  END SELECT

  IF ( BotTop == 'TOP' ) G = -G

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
!**********************************************************************!
SUBROUTINE ELASUP( x, yV, IPower, Medium )

  ! Propagates through an elastic layer using compound matrix formulation

  USE krakmod

  IMPLICIT NONE
  INTEGER,       PARAMETER :: IPowerR = 50, IPowerF = -50
  REAL (KIND=8), PARAMETER :: Roof = 1.0D50, Floor = 1.0D-50
  INTEGER                  :: IPower, Medium, II, j
  REAL (KIND=8)            :: x, two_h, two_x, four_h_x, xB3
  REAL (KIND=8)            :: xV( 5 ), yV( 5 ), zV( 5 )


  ! Euler's method for first step

  two_x    = 2.0D0 * x
  two_h    = 2.0D0 * h( Medium )
  four_h_x = 4.0D0 * h( Medium ) * x
  j        = Loc( Medium ) + N( Medium ) + 1
  xB3      = x * B3( j ) - rho( j )

  zV( 1 ) = yV( 1 ) - 0.5D0 * (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
  zV( 2 ) = yV( 2 ) - 0.5D0 * ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
  zV( 3 ) = yV( 3 ) - 0.5D0 * (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
  zV( 4 ) = yV( 4 ) - 0.5D0 * (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) -two_x * B4( j ) * yV( 3 ) )
  zV( 5 ) = yV( 5 ) - 0.5D0 * (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -       four_h_x * yV( 3 ) )

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

     !         P1 = yV( 2 ) * ( yV( 4 ) - yV( 5 ) )
     !         P0 = xV( 2 ) * ( xV( 4 ) - xV( 5 ) )
     !         IF ( ( P0 > 0.0 ) .AND. ( P1 < 0.0 ) ) ModeCount = ModeCount+1
     !         IF ( ( P0 < 0.0 ) .AND. ( P1 > 0.0 ) ) ModeCount = ModeCount+1

     ! Scale if necessary
     IF ( II /= 1 ) THEN
        IF ( ABS( zV( 2 ) ) < Floor ) THEN
           zV     = Roof * zV
           yV     = Roof * yV
           IPower = IPower - IPowerR
        ENDIF

        IF ( ABS( zV( 2 ) ) > Roof  ) THEN
           zV     = Floor * zV
           yV     = Floor * yV
           IPower = IPower - IPowerF
        ENDIF
     ENDIF
  END DO

  yV = ( xV + 2.0D0 * yV + zV ) / 4.0D0   ! Apply the standard filter at the terminal point

END SUBROUTINE ELASUP
!**********************************************************************!
SUBROUTINE ELASDN( x, yV, IPower, Medium )

  ! Propagates through an elastic layer using compound matrix formulation

  USE krakmod

  IMPLICIT NONE
  INTEGER,       PARAMETER :: IPowerR = 50, IPowerF = -50
  REAL (KIND=8), PARAMETER :: Roof = 1.0D50, Floor = 1.0D-50
  INTEGER                  :: IPower, Medium, II, j
  REAL (KIND=8)            :: x, two_h, two_x, four_h_x, xB3
  REAL (KIND=8)            :: xV( 5 ), yV( 5 ), zV( 5 )

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
        IF ( ABS( zV( 2 ) ) < Floor ) THEN
           zV     = Roof * zV
           yV     = Roof * yV
           IPower = IPower - IPowerR
        ENDIF

        IF ( ABS( zV( 2 ) ) > Roof  ) THEN
           zV     = Floor * zV
           yV     = Floor * yV
           IPower = IPower - IPowerF
        ENDIF
     ENDIF
  END DO

  yV = ( xV + 2.0D0 * yV + zV ) / 4.0D0   ! Apply the standard filter at the terminal point

END SUBROUTINE ELASDN
