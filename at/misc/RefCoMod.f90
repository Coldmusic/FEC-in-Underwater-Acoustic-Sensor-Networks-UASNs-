MODULE RefCoMod

  ! reflection coefficient data

  SAVE
  INTEGER, PARAMETER            :: BRCFile = 31, TRCFile = 32, IRCFile = 12
  INTEGER                       :: NBotPts, NTopPts, NkTab
  INTEGER,          ALLOCATABLE :: ITab( : )
  REAL    (KIND=8), ALLOCATABLE :: XTab( : )
  COMPLEX (KIND=8), ALLOCATABLE :: FTab( : ), GTab( : )

  TYPE ReflectionCoef
      REAL(KIND=8) :: theta, R, phi
  END TYPE

  TYPE(ReflectionCoef), ALLOCATABLE :: RBot( : ), RTop( : )

CONTAINS

  SUBROUTINE READRC( FileRoot, BotRC, TopRC, PRTFile )

    ! Optionally read in bottom reflection coefficient

    IMPLICIT NONE
    REAL ( KIND = 8 ), PARAMETER :: pi = 3.14159265, DegRad = pi / 180.0
    INTEGER               PRTFile, I, IK, IOStat, iAllocStat
    REAL  ( KIND = 8 ) :: Freq
    CHARACTER (LEN=1 ) :: BotRC, TopRC
    CHARACTER (LEN=80) :: Title2, FileRoot

    IF ( BotRC == 'F' ) THEN
       WRITE( PRTFile, * ) 'Using tabulated bottom reflection coef.'
       OPEN( FilE = TRIM( FileRoot ) // '.brc', UNIT = BRCFile, STATUS = 'OLD', IOSTAT = IOStat )
       IF ( IOStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'GETPAR', 'Unable to open Bottom Reflection Coefficient file' )

       READ(  BRCFile, * ) NBotPTS
       WRITE( PRTFile, * ) 'Number of points in bottom reflection coefficient = ', NBotPTS

       IF ( ALLOCATED( RBot ) ) DEALLOCATE( RBot )
       ALLOCATE( RBot( NBotPTS ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
            CALL ERROUT( PRTFile, 'F', 'READRD', 'Insufficient memory for bot. refl. coef.: reduce # points'  )

       READ(  BRCFile, * ) ( RBot( I ), I = 1, NBotPTS )
       CLOSE( BRCFile )
       RBot%phi = DegRad * RBot%phi   ! convert to radians

    ELSE   ! should allocate something anyway, since variable is passed
       ALLOCATE(  RBot( 1 ), Stat = IAllocStat )
    ENDIF

    ! Optionally read in top reflection coefficient

    IF ( TopRC == 'F' ) THEN
       WRITE( PRTFile, * ) 'Using tabulated top    reflection coef.'
       OPEN( FILE = TRIM( FileRoot ) // '.trc', UNIT = TRCFile, STATUS = 'OLD', IOSTAT = IOStat )
       IF ( IOStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'READRC', 'Unable to open Top Reflection Coefficient file' )

       READ(  TRCFile, * ) NTopPTS
       WRITE( PRTFile, * ) 'Number of points in top reflection coefficient = ', NTopPts

       IF ( ALLOCATED( RTop ) ) DEALLOCATE( RTop )
       ALLOCATE( RTop( NTopPTS ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'READRC', 'Insufficient memory for top refl. coef.: reduce # points'  )

       READ(  TRCFile, * ) ( RTop( I ), I = 1, NTopPTS )
       CLOSE( TRCFile )
       RTop%phi = DegRad *  RTop%phi   ! convert to radians
    ELSE   ! should allocate something anyway, since variable is passed
       ALLOCATE( RTop( 1 ), Stat = IAllocStat )
    ENDIF

    ! Optionally read in internal reflection coefficient data

    IF ( BotRC == 'P' ) THEN
       WRITE( PRTFile, * ) 'Reading precalculated refl. coeff. table'
       OPEN( FILE = TRIM( FileRoot ) // '.irc', UNIT = IRCFile, STATUS = 'OLD', IOSTAT = IOStat )
       IF ( IOsTAT /= 0 ) CALL ERROUT( PRTFile, 'F', 'READRC', 'Unable to open Internal Reflection Coefficient file' )

       READ(  IRCFile, * ) Title2, Freq
       READ(  IRCFile, * ) NkTab
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of points in internal reflection coefficient = ', NkTab

       IF ( ALLOCATED( XTab ) )  DEALLOCATE( XTab, FTab, GTab, ITab )
       ALLOCATE( XTab( NkTab ), FTab( NkTab ), GTab( NkTab ), ITab( NkTab ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'GETPAR', 'Too many points in reflection coefficient' )

       READ( IRCFile, FMT = "( 5G15.7, I5 )" ) ( XTab( IK ), FTab( IK ), GTab( IK ), ITab( IK ), IK = 1, NkTab )
       CLOSE( IRCFile )
    ENDIF

  END SUBROUTINE READRC

  SUBROUTINE REFCO( RInt, R, NPts, PRTFile )

    ! Given an angle ThetaInt, REFCO returns the magnitude and
    ! phase of the reflection coefficient (RInt%R, RInt%phi).

    ! Uses linear interpolation using the two nearest abscissas
    ! Assumes phi has been unwrapped so that it varies smoothly.
    ! I tried modifying it to allow a complex angle of incidence but
    ! stopped when I realized I needed to fuss with a complex atan2 routine

    IMPLICIT NONE
    INTEGER PRTFile, NPts, ILeft, IRight, IMid
    REAL ( KIND=8 ) :: alpha, Thetaintr
    TYPE(ReflectionCoef), INTENT( IN    ) :: R( NPts )
    TYPE(ReflectionCoef), INTENT( INOUT ) :: RInt

    ILeft  = 1
    IRight = NPts

    thetaIntr = REAL( RInt%Theta )   ! This should be unnecessary? probably used when I was doing complex angles

    ! Three cases: ThetaInt left, in, or right of tabulated interval

    IF     ( thetaIntr < R( ILeft  )%theta ) THEN
       !IRight = 2
       RInt%R   = 0.0     ! R( ILeft  )%R
       RInt%phi = 0.0     ! R( ILeft  )%phi
       CALL ERROUT( PRTFile, 'W', 'REFCO', 'Refl. Coef. being set to 0 outside tabulated domain' )
       WRITE( PRTFile, * ) 'angle = ', thetaintr, 'lower limit = ', R( ILeft)%theta

    ELSE IF( thetaIntr > R( IRight )%theta ) THEN
       !ILeft = NPts - 1
       RInt%R   = 0.0     ! R( IRight )%R
       RInt%phi = 0.0     ! R( IRight )%phi
       CALL ERROUT( PRTFile, 'W', 'REFCO', 'Refl. Coef. being set to 0 outside tabulated domain' )
       WRITE( PRTFile, * ) 'angle = ', ThetaIntr, 'upper limit = ', R( IRight )%theta

    ELSE
       ! Search for bracketting abscissas: Log2( NPts ) stabs required for a bracket

       DO WHILE ( ILeft /= IRight - 1 )
          IMid = ( ILeft + IRight ) / 2
          IF ( R( IMid )%theta > thetaIntr ) THEN
             IRight = IMid
          ELSE
             ILeft  = IMid
          ENDIF
       END DO

       ! Linear interpolation for reflection coef

       Alpha    = ( RInt%theta - R( ILeft )%theta ) / ( R( IRight )%theta - R( ILeft )%theta )
       RInt%R   = ( 1 - Alpha ) * R( ILeft )%R   + Alpha * R( IRight )%R
       RInt%phi = ( 1 - Alpha ) * R( ILeft )%phi + Alpha * R( IRight )%phi

    ENDIF

  END SUBROUTINE REFCO

  SUBROUTINE IRCINT( X, F, G, IPower, XTab, FTab, GTab, ITab, NkTab )

    ! Internal reflection coefficient interpolator.
    ! Returns F, G, IPower for given X using tabulated values.
    ! Uses polynomial interpolation using the two nearest abscissas.
    ! The order of the polynomial is set by N below

    IMPLICIT NONE

    INTEGER N
    PARAMETER ( N = 3 )

    INTEGER               I, J, IDel, ILeft, IMid, IRight, NAct, NkTab, IPower
    INTEGER               ITab( NkTab )
    REAL    ( KIND=8 ) :: XTab( NkTab ), Ten = 10.0, Xr
    COMPLEX ( KIND=8 ) :: FTab( NkTab ), GTab( NkTab ), XT( N ), FT( N ), GT( N ), PZ, X, F, G

    Xr     = DBLE( X )        ! taking the real part
    ILeft  = 1
    IRight = NkTab

    ! Three cases: X left, in, or right of tabulated interval

    IF (     Xr < XTab( ILeft  ) ) THEN
       F      = FTab( ILeft  )
       G      = GTab( ILeft  )
       IPower = ITab( ILeft  )
    ELSE IF( Xr > XTab( IRight ) ) THEN
       F      = FTab( IRight )
       G      = GTab( IRight )
       IPower = ITab( IRight )
    ELSE

       ! Search for bracketting abscissas:
       ! Log base 2 (NPTS) stabs required for a bracket

       DO WHILE ( ILeft /= IRight-1 )
          IMid = ( ILeft + IRight )/2
          IF ( XTab( IMid ) > Xr ) THEN
             IRight = IMid
          ELSE
             ILeft = IMid
          ENDIF
       END DO

       ! Extract the subset for interpolation and scale

       ILeft  = MAX( ILeft  - ( N - 2 ) / 2, 1 )
       IRight = MIN( IRight + ( N - 1 ) / 2, NkTab )

       NAct = IRight - ILeft + 1
       DO I = 1, NAct
          J       = I + ILeft - 1
          IDel    = ITab( J ) - ITab( ILeft )
          XT( I ) = XTab( J )
          FT( I ) = FTab( J ) * Ten ** IDel
          GT( I ) = GTab( J ) * Ten ** IDel
       END DO

       ! Construct the interpolate

       F      = PZ( X, XT, FT, NAct )
       G      = PZ( X, XT, GT, NAct )
       IPower = ITab( ILeft )
    ENDIF

  END SUBROUTINE IRCINT
END MODULE RefComod
