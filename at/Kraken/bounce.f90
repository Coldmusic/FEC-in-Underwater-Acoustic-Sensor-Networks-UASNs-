PROGRAM BOUNCE

  ! Program for computing the reflection coefficient, R(theta)
  ! Michael B. Porter

  USE krakcmod
  USE RefCoMod
  IMPLICIT NONE
  INTEGER            :: IAllocStat
  REAL (KIND=8)      :: RKMin, RKmax
  CHARACTER (LEN=80) :: FileRoot

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  ! *** Read in environmental info ***

  TITLE = 'BOUNCE-'
  iProf = 1

  CALL READIN( FileRoot, TITLE, Freq, MaxMedium, NMedia, &
       TopOpt, HSTop, NG, Sigma, Depth,     &
       BotOpt, HSBot, ENVFile, PRTFile )

  READ(  ENVFile, *    ) Clow, Chigh        ! Spectral limits
  WRITE( PRTFile, "( /, ' cLow = ', G10.5, 'm/s      cHigh = ', G10.5, 'm/s' )" ) cLow, cHigh

  IF ( Clow == 0.0 ) CALL ERROUT( PRTFile, 'F', 'BOUNCE', 'Clow must be greater than zero' )

  READ(  ENVFile, * ) RMax                  ! Maximum range for calculations
  WRITE( PRTFile, * ) 'RMax = ', RMax
  CLOSE( ENVFile )

  omega2 = ( 2.0 * pi * Freq ) ** 2

  CALL READRC( FileRoot, BotOpt( 1 : 1 ), TopOpt( 1 : 1 ), PRTFile ) ! Optionally read in bottom reflection coefficient data

  ! *** Write internal reflection coefficient data ***

  WRITE( PRTFile, * ) 'Writing internal refl. coeff. table'
  ! Compute NkTab

  RKmin = SQRT( omega2 ) / Chigh
  RKmax = SQRT( omega2 ) / Clow
  IF ( Chigh > 1.0E6 ) RKmin = 0.0

  NkTab = INT( 1000.0 * RMax * ( RKmax - RKmin ) / ( 2.0 * pi ) )
  WRITE( PRTFile, * ) 'NkTab = ', NkTab

  ALLOCATE( XTab( NkTab ), FTab( NkTab ), GTab( NkTab), ITab( NkTab ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BOUNCE', 'Too many points in reflection coefficient' )

  CALL COMPR( FileRoot, RKmin, RKmax )

END PROGRAM BOUNCE

!**********************************************************************C

SUBROUTINE INIT

  ! Initializes arrays defining difference equations

  USE krakcmod
  IMPLICIT NONE
  LOGICAL           :: ElasticFlag = .FALSE.
  INTEGER           :: iAllocStat, ii, j, Med, N1, NPts
  REAL (KIND=8)     :: TWOH
  COMPLEX  (KIND=8) :: CP2, CS2
  COMPLEX (KIND=8), ALLOCATABLE :: CP( : ), CS( : )
  CHARACTER (LEN=8) :: TASK

  Cmin          = 1.0E6
  FirstAcoustic = 0
  LOC( 1 )      = 0
  NPTS          = SUM( N( 1 : NMedia ) ) + NMedia

  ALLOCATE ( B1( NPTS ), B2( NPTS ), B3( NPTS ), B4( NPTS ), rho( NPTS ), CP( NPTS ), CS( NPTS ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) &
       CALL ERROUT( PRTFile, 'F', 'BOUNCE - INIT', 'Insufficient memory: Reduce mesh.' )

  ! *** Loop over media ***

  DO Med = 1, NMedia
     IF ( Med /= 1 ) LOC( Med ) = LOC( Med - 1 ) + N( Med - 1 ) + 1
     N1 = N(   Med ) + 1
     II = LOC( Med ) + 1

     ! *** PROFIL reads in the data for a medium ***

     TASK = 'TAB'
     CALL PROFIL( Depth, CP( II ), CS( II ), rho( II ), Med, N1, Freq, TopOpt( 1 : 1 ), TopOpt( 3 : 4 ), TASK, ENVFile, PRTFile )

     ! *** Load diagonals of the finite-difference equations ***

     IF ( CS( II ) == ( 0.0, 0.0 ) ) THEN ! Case of an acoustic medium ---

        Material( Med ) = 'ACOUSTIC'
        IF ( FirstAcoustic == 0 ) FirstAcoustic = Med
        LastAcoustic = Med

        Cmin = MIN( MINVAL( DBLE( CP( II : II + N( Med ) ) ) ), Cmin )
        B1( II : II + N( Med ) ) = -2.0 +  h( Med ) ** 2 * omega2 / CP( II : II + N( Med ) ) ** 2

     ELSE                                ! Case of an elastic medium ---

        IF ( Sigma( Med ) /= 0.0 ) THEN
           WRITE( PRTFile, * ) 'Rough elastic interface not allowed'
           WRITE( PRTFile, * ) 'PROGRAM ABORTING'
           STOP 'ERROR IN BOUNCE: Rough elastic interface not allowed'
        ENDIF

        Material( Med ) = 'ELASTIC'
        ElasticFlag = .TRUE.
        TWOH   = 2.0 * h( Med )

        DO J = II, II + N( Med )
           Cmin = MIN( DBLE( CS( J ) ), Cmin )

           CP2 = CP( J ) ** 2
           CS2 = CS( J ) ** 2

           B1( J )  = TWOH / ( rho( J ) * CS2 )
           B2( J )  = TWOH / ( rho( J ) * CP2 )
           B3( J )  = 4.0 * TWOH * rho( J ) * CS2 * ( CP2 - CS2 ) / CP2
           B4( J )  = TWOH * ( CP2 - 2.0 * CS2 ) / CP2
           rho( J ) = TWOH * omega2 * rho( J )
        END DO
     ENDIF
  END DO   ! Next Med

  ! *** Bottom properties ***

  IF ( HSBot%BC( 1 : 1 ) == 'A' ) THEN
     IF ( HSBot%CS /= ( 0.0, 0.0 ) ) THEN ! Elastic bottom:
        ElasticFlag = .TRUE.
        Cmin = MIN( Cmin, DBLE(  HSBot%CS ) )
     ELSE                            ! Acoustic bottom:
        Cmin = MIN( Cmin, DBLE(  HSBot%CP ) )
     ENDIF
  ENDIF

  ! *** Top properties ***

  IF ( HSTop%BC( 1 : 1 ) == 'A' ) THEN
     IF (  HSTop%CS /= ( 0.0, 0.0 ) ) THEN   ! Elastic top:
        ElasticFlag = .TRUE.
        Cmin = MIN( Cmin, DBLE(  HSTop%CS ) )
     ELSE                              ! Acoustic top:
        Cmin = MIN( Cmin, DBLE(  HSTop%CP ) )
     ENDIF
  ENDIF

  IF ( ElasticFlag ) Cmin = 0.9 * Cmin
  Clow = MAX( Clow, 0.99 * Cmin )

END SUBROUTINE INIT
!**********************************************************************C
SUBROUTINE COMPR( FileRoot, RKmin, RKmax )

  ! Computes the reflection coefficient for k in [RKmin, RKmax]

  USE krakcmod
  USE RefCoMod
  IMPLICIT NONE
  REAL (KIND=8), PARAMETER :: RadDeg = 180.0 / pi
  INTEGER                  :: ik, IPow, itheta, Loops
  REAL      (KIND=8)       :: K0, C0, Deltak, RK, RKMin, RKMax
  COMPLEX   (KIND=8)       :: X, F, G
  REAL     (KIND=8), ALLOCATABLE :: Kx( : ), Kz( : ), Theta( : ), R( : ), phase( : )
  COMPLEX  (KIND=8), ALLOCATABLE :: RCmplx( : )
  CHARACTER (LEN=80)       :: FileRoot

  ALLOCATE( Kx( NkTab ), Kz( NkTab ), RCmplx( NkTab) , R( NkTab ), Theta( NkTab ), phase( NkTab ) )

  N( 1 : NMedia ) = NG( 1 : NMedia )
  h( 1 : NMedia ) = ( Depth( 2 : NMedia + 1 ) - Depth( 1 : NMedia ) ) / N( 1 : NMedia )

  HV( 1 ) = h( 1 )
  CALL INIT

  DeltaK = ( RKmax - RKmin ) / ( NkTab - 1 )

  DO ik = 1, NkTab
     RK = RKmin + ( ik - 1 ) * DeltaK
     X  = RK ** 2

     CALL BCImpedance( X, 'BOT', HSBot, F, G, IPow )  ! Bottom impedance
     CALL ACOUST( X, F, G, IPow  )  ! Shoot through acoustic layers
     XTab( ik ) = DBLE( X )
     FTab( ik ) = F
     GTab( ik ) = G
     ITab( ik ) = IPow
  END DO

  IF ( HSTop%BC( 1 : 1 ) == 'A' ) THEN
     C0 = DBLE( HSTop%CP )                    ! use upper halfspace speed for reference if a halfspace was specified
  ELSE
     C0 = 1500
  END IF

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Reference sound speed = ', C0

  K0    = SQRT( omega2 ) / C0 ! free-space wavenumber
  Kx    = SQRT( XTab )        ! horizontal wavenumber
  R     = 0.0
  phase = 0.0
  Kz    = 0

  WHERE( K0 > Kx )
     Kz = SQRT( K0 ** 2 - Kx ** 2 ) ! vertical   wavenumber
  END WHERE
  Theta  = RadDeg * ATAN2( Kz, Kx )    ! angle of incidence
  RCmplx =  - ( FTab - i * Kz * GTab ) / ( FTab + i * Kz * GTab )   ! complex reflection coef.
  R      = ABS( RCmplx )
  phase  = RadDeg * ATAN2( AIMAG( RCmplx ), REAL( RCmplx ) )

  ! unwrap the phase by counting loops in the complex plane
  Loops = 0
  DO itheta = NkTab - 1, 1, -1
     IF ( AIMAG( RCmplx( itheta ) ) > 0 .AND. AIMAG( RCmplx( itheta + 1 ) ) < 0 .AND. &
                                               REAL( RCmplx( itheta + 1 ) ) < 0 ) Loops = Loops + 1
     IF ( AIMAG( RCmplx( itheta ) ) < 0 .AND. AIMAG( RCmplx( itheta + 1 ) ) > 0 .AND. &
                                               REAL( RCmplx( itheta + 1 ) ) < 0 ) Loops = Loops - 1
     phase( itheta ) = phase( itheta ) - Loops * 360
  END DO

  OPEN( FILE = TRIM( FileRoot ) // '.irc', UNIT = IRCFile, STATUS = 'UNKNOWN' )  ! Internal Reflection Coef. format
  WRITE( IRCFile, * ) '''', TITLE, '''', Freq
  WRITE( IRCFile, * ) NkTab
  WRITE( IRCFile, FMT = "( 5G15.7, I5 )" ) ( XTab( ik ), FTab( ik ), GTab( ik ), ITab( ik ), ik = 1, NkTab )

  OPEN( FILE = TRIM( FileRoot ) // '.brc', UNIT = BRCFile, STATUS = 'UNKNOWN' )  ! Bottom Reflection Coef. format
  WRITE( BRCFile, * ) NkTab
  DO ik = NkTab, 1, -1
     WRITE( BRCFile, * ) Theta( ik ), R( ik ), phase( ik )
  END DO

END SUBROUTINE COMPR
!**********************************************************************C
SUBROUTINE ACOUST( X, F, G, IPow )

  ! Shoot through acoustic layers

  USE krakcmod
  IMPLICIT NONE
  INTEGER,       PARAMETER :: IPowF = -5
  REAL (KIND=8), PARAMETER :: ROOF = 1.0E5, FLOOR = 1.0E-5
  INTEGER                  :: II, IPow, Med
  REAL (KIND=8)            :: rhoM
  COMPLEX (KIND=8)         :: X, F, G, P0, P1, P2, H2K2

  IF ( FirstAcoustic == 0 ) RETURN

  ! *** Loop over successive acoustic media ***

  DO Med = LastAcoustic, FirstAcoustic, -1
     H2K2 = h( Med ) ** 2 * X
     II   = LOC( Med ) + N( Med ) + 1
     rhoM = rho( II )
     P1   = -2.0 * G
     P2   = ( B1( II ) - H2K2 ) * G - 2.0 * h( Med ) * F * rhoM

     ! *** Shoot through a single medium ***

     DO II = LOC( Med ) + N( Med ), LOC( Med ) + 1, -1

        P0 = P1
        P1 = P2
        P2 = ( H2K2 - B1( II ) ) * P1 - P0

        DO WHILE ( ABS( DBLE( P2 ) ) > ROOF )
           P0   = FLOOR * P0
           P1   = FLOOR * P1
           P2   = FLOOR * P2
           IPow = IPow - IPowF
        END DO
     END DO

     ! F = P'/rho and G = -P since FP+GP'/rho = 0
     F = -( P2 - P0 ) / ( 2.0 * h( Med ) ) / rhoM
     G = -P1
  END DO

END SUBROUTINE ACOUST
