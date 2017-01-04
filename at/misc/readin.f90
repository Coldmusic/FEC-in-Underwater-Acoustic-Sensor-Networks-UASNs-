SUBROUTINE READIN( FileRoot, Title, Freq, MaxMedia, NMedia, &
     TopOpt, HSTop, NG, sigma, Depth, BotOpt, HSBot, ENVFile, PRTFile )

  ! Reads in the info in ENVFile

  IMPLICIT NONE
  INTEGER, PARAMETER  :: MaxSSP = 2001
  INTEGER                ENVFile, PRTFile, NMedia, NG( * ), MaxMedia, NElts, Medium, Nneeded, iostat
  REAL       (KIND=8) :: alphaR, betaR, alphaI, betaI, rhoR, Freq, &
                         sigma( * ), Depth( * ), rho( 1 ), C, deltaz
  COMPLEX    (KIND=8) :: cP( MaxSSP ), cS( MaxSSP )
  CHARACTER              TopOpt*( * ), BotOpt*( * ), Title*( * )
  CHARACTER ( LEN=1 ) :: SSPType
  CHARACTER ( LEN=2 ) :: AttenUnit
  CHARACTER ( LEN=8 ) :: Task
  CHARACTER ( LEN=* ) :: FileRoot

   ! Halfspace properties
   TYPE HSInfo
      CHARACTER (LEN=1)          :: BC       ! Boundary condition type
      COMPLEX (KIND=8)           :: cP, cS   ! P-wave, S-wave speeds
      REAL    (KIND=8)           :: rho      ! density
   END TYPE

   TYPE( HSInfo )              :: HSTop, HSBot

  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  ! Open the environmental file
  OPEN( UNIT = ENVFile, FILE = TRIM( FileRoot ) // '.env', STATUS = 'OLD', IOSTAT = iostat )
  IF ( IOSTAT /= 0 ) THEN   ! successful open?
     WRITE( PRTFile, * ) 'ENVFile = ', TRIM( FileRoot ) // '.env'
     CALL ERROUT( PrtFile, 'F', 'READIN', 'Unable to open the environmental file' )
  END IF

  ! Open the print file
  OPEN( UNIT = PRTFile, FILE = TRIM( FileRoot ) // '.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )
  IF ( IOSTAT /= 0 ) THEN   ! successful open?
     WRITE( PRTFile, * ) 'PRTFile = ', TRIM( FileRoot ) // '.prt'
     CALL ERROUT( PrtFile, 'F', 'BELLHOP - READIN', 'Unable to open the print file' )
  END IF

  alphaR = 1500.0
  betaR  = 0.0
  rhoR   = 1.0
  alphaI = 0.0
  betaI  = 0.0
  NElts  = 0         ! this is a dummy variable, passed to profil during read of SSP

  WRITE( PRTFile, * ) '_________________________________________________'
  WRITE( PRTFile, * )
  READ(  ENVFile, *, END = 9999 ) Title( 9 : 80 )
  WRITE( PRTFile, * ) Title

  READ( ENVFile, *, END = 9999  ) Freq
  READ( ENVFile, *, END = 9999  ) NMedia
  WRITE( PRTFile, "( ' Frequency = ', G11.4, 'Hz   NMedia = ', I3, // )" ) Freq, NMedia

  IF ( NMedia > MaxMedia ) THEN
     WRITE( PRTFile, * ) 'MaxMedia = ', MaxMedia
     CALL ERROUT( PRTFile, 'F', 'READIN', 'Too many Media' )
  ENDIF

  ! TOP OPTIONS
  READ( ENVFile, * ) TopOpt( 1 : 8 )
  SSPType          = TopOpt( 1 : 1 )
  HSTop%BC         = TopOpt( 2 : 2 )
  AttenUnit        = TopOpt( 3 : 4 )

  SELECT CASE ( SSPType )        ! SSP approximation options
  CASE ( 'N' )
     WRITE( PRTFile, * ) '    N2-LINEAR approximation to SSP'
  CASE ( 'C' )
     WRITE( PRTFile, * ) '    C-LINEAR approximation to SSP'
  CASE ( 'S' )
     WRITE( PRTFile, * ) '    SPLINE approximation to SSP'
  CASE ( 'A' )
     WRITE( PRTFile, * ) '    ANALYTIC SSP option'
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'READIN', 'Unknown option for SSP approximation' )
  END SELECT

  SELECT CASE ( AttenUnit( 1 : 1 ) ) ! Attenuation options
  CASE ( 'N' )
     WRITE( PRTFile, * ) '    Attenuation units: nepers/m'
  CASE ( 'F' )
     WRITE( PRTFile, * ) '    Attenuation units: dB/mkHz'
  CASE ( 'M' ) 
     WRITE( PRTFile, * ) '    Attenuation units: dB/m'
  CASE ( 'W' )
     WRITE( PRTFile, * ) '    Attenuation units: dB/wavelength'
  CASE ( 'Q' )
     WRITE( PRTFile, * ) '    Attenuation units: Q'
  CASE ( 'L' )
     WRITE( PRTFile, * ) '    Attenuation units: Loss parameter'
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'READIN', 'Unknown attenuation units' )
  END SELECT

  SELECT CASE ( AttenUnit( 2 : 2 ) ) !  Added volume attenuation
  CASE ( 'T' )
     WRITE( PRTFile, * ) '    THORP attenuation added'
  END SELECT

  ! CALL TopBot  to read top BC 
  IF ( HSTop%BC == 'A' ) &
       WRITE( PRTFile, "( //, '   z (m)     alphaR (m/s)     betaR   rho (g/cm^3)  alphaI     betaI', / )" )
  CALL TopBot( ENVFile, PRTFile, Freq, AttenUnit, HSTop )

  !  Internal media 
  IF ( HSTop%BC /= 'A' ) &
       WRITE( PRTFile, "( //, '   z (m)     alphaR (m/s)     betaR   rho (g/cm^3)  alphaI     betaI', / )" )

  DO Medium = 1, NMedia
     READ(  ENVFile, *, END = 9999 ) NG( Medium ), sigma( Medium ), Depth( Medium + 1 )
     WRITE( PRTFile, "( /, '          ( Number of pts = ', I5, '  RMS roughness = ', G10.3, ' )')" ) &
          NG( Medium ),            sigma( Medium )

     !  Call PROFIL to read in SSP 
     Task = 'INIT'
     CALL PROFIL( Depth, cP, cS, rho, Medium, NElts, Freq, SSPType, AttenUnit, Task, ENVFile, PRTFile  )

     ! estimate number of points needed
     C = alphar
     IF ( betar > 0.0 ) C = betar     ! shear?
     deltaz  = C / Freq / 20           ! default sampling: 20 points per wavelength
     Nneeded = INT( ( Depth( Medium + 1 ) - Depth( Medium ) ) / deltaz )
     Nneeded = MAX( Nneeded, 10 )     ! require a minimum of 10 points
  
     IF ( NG( Medium ) == 0 ) THEN    ! automatic calculation of f.d. mesh
        NG( Medium ) = Nneeded
        WRITE( PRTFile, * ) 'Number of pts = ', NG( Medium )
     ELSEIF ( NG( Medium ) < Nneeded/2 ) THEN
       CALL ERROUT( PRTFile, 'F', 'READIN', 'Mesh is too coarse' )
     END IF

  END DO   ! next Medium

  ! Bottom properties 
  READ( ENVFile, *, END = 9999 ) BotOpt( 1 : 8 ), sigma( NMedia + 1 )
  HSBot%BC = BotOpt( 1 : 1 )
  WRITE( PRTFile, * )
  WRITE( PRTFile, "( 33X, '( RMS roughness = ', G10.3, ' )' )" ) sigma( NMedia + 1 )

  ! CALL TopBot  to read bottom BC 
  CALL TopBot( ENVFile, PRTFile, Freq, AttenUnit, HSBot )

  RETURN

9999 WRITE( PRTFile, * ) 'End of environmental file'

  CLOSE( ENVFile )
  STOP

END SUBROUTINE READIN
!**********************************************************************!
SUBROUTINE TopBot( ENVFile, PRTFile, Freq, AttenUnit, HS )

  ! Handles top and bottom boundary conditions

  ! Input:
  !     ENVFile: Environmental file
  !     PRTFile: Print file
  !     Freq:    Frequency
  !     HS%BC:  Boundary condition type
  !
  ! Output:
  !    HS%cP:    P-wave speed in halfspace
  !    HS%cS:    S-wave speed in halfspace
  !    HS%rho:   density in halfspace

  !    BumDen:   Bump density
  !    eta:      Principal radius 1
  !    xi:       Principal radius 2

  IMPLICIT NONE
  INTEGER           :: ENVFile, PRTFile
  REAL     (KIND=8) :: alphaR, betaR, alphaI, betaI, rhoR, Freq, zTemp
  COMPLEX  (KIND=8) :: CRCI
  CHARACTER (LEN=2) :: AttenUnit

   ! Halfspace properties
   TYPE HSInfo
      CHARACTER (LEN=1)          :: BC                          ! Boundary condition type
      COMPLEX (KIND=8)           :: cP, cS                      ! P-wave, S-wave speeds
      REAL    (KIND=8)           :: rho, BumpDensity, eta, xi   ! density, boss parameters

   END TYPE

   TYPE( HSInfo )              :: HS

  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  ! Echo to PRTFile user's choice of boundary condition 
  SELECT CASE ( HS%BC )
  CASE ( 'S' )
     WRITE( PRTFile, * ) '    Twersky SOFT BOSS scatter model'
  CASE ( 'H' )
     WRITE( PRTFile, * ) '    Twersky HARD BOSS scatter model'
  CASE ( 'T' )
     WRITE( PRTFile, * ) '    Twersky (amplitude only) SOFT BOSS scatter model'
  CASE ( 'I' )
     WRITE( PRTFile, * ) '    Twersky (amplitude only) HARD BOSS scatter model'
  CASE ( 'V' )
     WRITE( PRTFile, * ) '    VACUUM'
  CASE ( 'R' )
     WRITE( PRTFile, * ) '    Perfectly RIGID'
  CASE ( 'A' )
     WRITE( PRTFile, * ) '    ACOUSTO-ELASTIC half-space'
  CASE ( 'F' )
     WRITE( PRTFile, * ) '    FILE used for reflection loss'
  CASE ( 'W' )
     WRITE( PRTFile, * ) '    Writing an IRC file'
  CASE ( 'P' )
     WRITE( PRTFile, * ) '    reading PRECALCULATED IRC'
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'TopBot', 'Unknown boundary condition type' )
  END SELECT

  ! Read in BC parameters depending on particular choice 
  HS%cP  = 0.0
  HS%cS  = 0.0
  HS%rho = 0.0

  SELECT CASE ( HS%BC )
  CASE ( 'S', 'H', 'T', 'I' )    ! Twersky ice model parameters 
     READ(  ENVFile, *    ) HS%BumpDensity, HS%eta, HS%xi
     WRITE( PRTFile, 1000 ) HS%BumpDensity, HS%eta, HS%xi
1000 FORMAT( /, ' Twersky ice model parameters:', /, &
          ' Bump Density = ', G15.6, '  Eta = ', G11.3, '  Xi = ', G11.3, /)
  CASE ( 'A' )                   !  Half-space properties 
     READ(  ENVFile, *    ) zTemp, alphaR, betaR, rhoR, alphaI, betaI
     WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
                        zTemp, alphaR, betaR, rhoR, alphaI, betaI

     HS%cP  = CRCI( alphaR, alphaI, Freq, AttenUnit )
     HS%cS  = CRCI( betaR,  betaI,  Freq, AttenUnit )
     HS%rho = rhoR
     IF ( HS%cP == 0.0 .OR. HS%rho == 0.0 ) &
          CALL ERROUT( PRTFile, 'F', 'TopBot', 'Sound speed or density vanishes in halfspace' )
  END SELECT

  RETURN
END SUBROUTINE TopBot
!**********************************************************************!
SUBROUTINE PROFIL( Depth, cP, cS, rhoT, Medium, N1, Freq, SSPType, AttenUnit, Task, ENVFile, PRTFile )

  ! Call the particular SSP routine specified by SSPType
  ! PROFIL is expected to perform two Tasks:
  !    Task = 'TAB'  then tabulate cP, cS, rhoT
  !    Task = 'INIT' then initialize
  ! Note that Freq is only need if Task = 'INIT'

  IMPLICIT NONE
  INTEGER              ENVFile, PRTFile, Medium, N1, I
  REAL     (KIND=8) :: rhoT( * ), Depth( * ), Freq, h, z
  COMPLEX  (KIND=8) :: cPT, cST
  COMPLEX  (KIND=8) :: cP( * ), cS( * )
  CHARACTER (LEN=1) :: SSPType
  CHARACTER (LEN=2) :: AttenUnit
  CHARACTER (LEN=8) :: Task

  SELECT CASE ( SSPType )
  CASE ( 'A' )  !  Analytic profile option 
     IF ( Task( 1 : 4 ) == 'INIT' ) THEN
        N1 = 21
        CALL ANALYT( Depth, cP, cS, rhoT, Medium, N1, Freq, AttenUnit, Task )
        H = ( Depth( Medium+1 ) - Depth( Medium ) ) / ( N1 - 1 )

        DO I = 1, N1
           z   = Depth( Medium ) + ( I - 1 ) * H
           cPT =  cP( I )
           cST =  cS( I )
           WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
              z,  REAL( cPT ),  REAL( cST ), rhoT( I ), AIMAG( cPT ), AIMAG( cST )
        END DO
     ELSE
        CALL ANALYT( Depth, cP, cS, rhoT, Medium, N1, Freq, AttenUnit, Task )
     ENDIF
  CASE ( 'N' )  !  N2-linear profile option 
     CALL N2LIN(  Depth, cP, cS, rhoT, Medium, N1, Freq, AttenUnit, Task, ENVFile, PRTFile )
  CASE ( 'C' )  !  C-linear profile option 
     CALL CLIN(   Depth, cP, cS, rhoT, Medium, N1, Freq, AttenUnit, Task, ENVFile, PRTFile )
  CASE ( 'S' )  !  Cubic spline profile option 
     CALL CCUBIC( Depth, cP, cS, rhoT, Medium, N1, Freq, AttenUnit, Task, ENVFile, PRTFile )
  CASE DEFAULT  !  Non-existent profile option 
     WRITE( PRTFile, * ) 'Profile option: ', SSPType(1:1)
     CALL ERROUT( PRTFile, 'F', 'PROFIL', 'Unknown profile option' )
  END SELECT

  RETURN
END SUBROUTINE PROFIL
!**********************************************************************!
FUNCTION CRCI( C, alpha, Freq, AttenUnit )

  ! Converts real wave speed and attenuation to a single
  !  complex wave speed (with positive imaginary part)

  ! 6 CASES:    N for Nepers/meter
  !             M for dB/meter      (M for Meters)
  !             F for dB/m-kHz      (F for Frequency dependent)
  !             W for dB/wavelength (W for Wavelength)
  !             Q for Q
  !             L for Loss parameter
  !
  ! second letter adds volume attenuation according to standard laws:
  !             T for Thorp

  IMPLICIT NONE
  REAL (KIND=8), PARAMETER :: pi = 3.1415926535897932D0
  REAL     (KIND=8) :: Freq, F2, omega, alpha, alphaT, C, Thorpe
  COMPLEX  (KIND=8) :: CRCI
  CHARACTER (LEN=2) :: AttenUnit

  omega = 2.0 * pi * Freq

  !  Convert to Nepers/m 
  alphaT = 0.0
  SELECT CASE ( AttenUnit(1:1) )
  CASE ( 'N' )
     alphaT = alpha
  CASE ( 'M' )
     alphaT = alpha / 8.6858896D0
  CASE ( 'F' )
     alphaT = alpha * Freq / 8685.8896D0
  CASE ( 'W' )
     IF ( C /= 0.0 ) alphaT = alpha * Freq / ( 8.6858896D0 * C )
     !        The following lines give f^1.25 Frequency dependence
     !        FAC = SQRT( SQRT( Freq / 50.0 ) )
     !        IF ( C /= 0.0 ) alphaT = FAC * alpha * Freq / ( 8.6858896D0 * C )
  CASE ( 'Q' )
     IF( C * alpha /= 0.0 ) alphaT = omega / ( 2.0 * C * alpha )
  CASE ( 'L' )   ! loss parameter
     IF ( C /= 0.0        ) alphaT = alpha * omega / C
  END SELECT

  ! added volume attenuation
  SELECT CASE ( AttenUnit(2:2) )
  CASE ( 'T' )
     F2 = ( Freq / 1000.0 ) **2

     ! Original formula from Thorp 1967
     ! Thorpe = 40.0 * F2 / ( 4100.0 + F2 ) + 0.1 * F2 / ( 1.0 + F2 )   ! dB/kyard
     ! Thorpe = Thorpe / 914.4D0                 ! dB / m
     ! Thorpe = Thorpe / 8.6858896D0             ! Nepers / m

     ! Updated formula from JKPS Eq. 1.34
     Thorpe = 3.3d-3 + 0.11 * f2 / ( 1.0 + f2 ) + 44.0 * f2 / ( 4100.0 + f2 ) + 3d-4* f2   ! dB/km
     Thorpe = alpha / 8685.8896 ! Nepers / m

     alphaT = alphaT + Thorpe
  END SELECT

  ! Convert Nepers/m to equivalent imaginary sound speed 
  alphaT = alphaT * C * C / omega
  CRCI   = CMPLX( C, alphaT, KIND=8 )

  RETURN
END FUNCTION CRCI
!**********************************************************************!
SUBROUTINE N2LIN( Depth, cP, cS, rhoT, Medium, N1, Freq, AttenUnit, Task, ENVFile, PRTFile )

  ! Tabulate cP, cS, rho for specified Medium
  ! Uses N2-linear segments for P and S-wave speeds
  ! Uses rho-linear segments for density

  IMPLICIT NONE
  INTEGER ENVFile, PRTFile, MaxMedia, MaxSSP
  PARAMETER ( MaxMedia = 501, MaxSSP = 2001 )
  INTEGER              Loc( MaxMedia ), NSSPPts( MaxMedia ), Medium, N, N1, I, ILoc, Lay
  REAL     (KIND=8) :: Depth( * ), rhoT( * ), z( MaxSSP ), rho( MaxSSP ), Freq, alphaR, betaR, rhoR, alphaI, betaI, H, R, zT
  COMPLEX  (KIND=8) :: cP( * ), cS( * ), alpha( MaxSSP ), beta( MaxSSP ), N2BOT, N2TOP, CRCI
  CHARACTER (LEN=2) :: AttenUnit
  CHARACTER (LEN=8) :: Task

  SAVE z, alpha, beta, rho, Loc, NSSPPts
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  ! If Task = 'INIT' then this is the first call and SSP is read.
  ! Any other call is a request for SSP subtabulation.

  IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization

     ! The variable Loc( Medium ) points to the starting point for the
     ! data in the arrays z, alpha, beta and rho
     IF ( Medium == 1 ) THEN
        Loc( Medium ) = 0
     ELSE
        Loc( Medium ) = Loc( Medium - 1 ) + NSSPPts( Medium - 1 )
     ENDIF
     ILoc = Loc( Medium )

     !  Read in data and convert attenuation to Nepers/m 
     N1 = 1
     DO I = 1, MaxSSP
        READ(  ENVFile, *    ) z( ILoc + I ), alphaR, betaR, rhoR, alphaI, betaI
        WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
                        z( ILoc + I ), alphaR, betaR, rhoR, alphaI, betaI

        alpha( ILoc + I ) = CRCI( alphaR, alphaI, Freq, AttenUnit )
        beta(  ILoc + I ) = CRCI( betaR,  betaI,  Freq, AttenUnit )
        rho(   ILoc + I ) = rhoR

        ! Did we read the last point?
        IF ( ABS( z( ILoc + I ) - Depth( Medium + 1 ) ) < EPSILON( 1.0e0 ) * Depth( Medium + 1 ) ) THEN
           NSSPPts( Medium ) = N1
           IF ( Medium == 1 ) Depth( 1 ) = z( 1 )
           RETURN
        ENDIF

        N1 = N1 + 1
     END DO

     ! Fall through means too many points in the profile
     WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
     CALL ERROUT( PRTFile, 'F', 'N2LIN', 'Number of SSP points exceeds limit' )

  ELSE   ! Task = 'TABULATE'
     ILoc = Loc( Medium )
     N    = N1 - 1
     H    = ( z( ILoc + NSSPPts( Medium ) ) - z( ILoc + 1 ) ) / N
     Lay  = 1

     DO I = 1, N1
        zT = z( ILoc + 1 ) + ( I - 1 ) * H
        IF ( I == N1 ) zT = z( ILoc + NSSPPts( Medium ) )   ! Make sure no overshoot

        DO WHILE ( zT > z( ILoc + Lay + 1 ) )
           Lay = Lay + 1
        END DO

        R = ( zT - z( ILoc + Lay ) ) / ( z( ILoc + Lay+1 ) - z( ILoc + Lay ) )

        ! P-wave
        N2TOP   = 1.0 / alpha( ILoc + Lay     )**2
        N2BOT   = 1.0 / alpha( ILoc + Lay + 1 )**2
        cP( I ) = 1.0 / SQRT( ( 1.0 - R ) * N2TOP + R * N2BOT )

        ! S-wave
        IF ( beta(ILoc + Lay) /= 0.0 ) THEN
           N2TOP   = 1.0 / beta( ILoc + Lay     )**2
           N2BOT   = 1.0 / beta( ILoc + Lay + 1 )**2
           cS( I ) = 1.0 / SQRT( ( 1.0 - R ) * N2TOP + R * N2BOT )
        ELSE
           cS( I ) = 0.0
        ENDIF

        rhoT( I ) = ( 1.0 - R ) * rho( ILoc + Lay ) + R * rho( ILoc + Lay + 1 )
     END DO

  ENDIF

  RETURN
END SUBROUTINE N2LIN
!**********************************************************************!
SUBROUTINE CLIN( Depth, cP, cS, rhoT, Medium, N1, Freq, AttenUnit, Task, ENVFile, PRTFile  )

  ! Tabulate cP, cS, rho for specified Medium

  ! Uses c-linear segments for P and S-wave speeds
  ! Uses rho-linear segments for density
  IMPLICIT NONE
  INTEGER ENVFile, PRTFile, MaxMedia, MaxSSP
  PARAMETER ( MaxMedia = 501, MaxSSP = 2001 )
  INTEGER              Loc( MaxMedia ), NSSPPts( MaxMedia ), Medium, N, N1, I, ILoc, Lay
  REAL    (KIND=8)  :: Depth( * ), rhoT( * ), z( MaxSSP ), rho( MaxSSP ), Freq, alphaR, betaR, rhoR, alphaI, betaI, H, R, zT
  COMPLEX (KIND=8)  :: cP( * ), cS( * ), alpha( MaxSSP ), beta( MaxSSP ), CRCI
  CHARACTER (LEN=2) :: AttenUnit
  CHARACTER (LEN=8) :: Task

  SAVE z, alpha, beta, rho, Loc, NSSPPts
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  ! If Task = 'INIT' then this is the first call and SSP is read.
  ! Any other call is a request for SSP subtabulation.

  IF ( Task(1:4) == 'INIT' ) THEN   ! Task 'INIT' FOR INITIALIZATION
     NSSPPts( Medium ) = N1

     ! The variable Loc(Medium) points to the starting point for the
     ! data in the arrays z, alpha, beta and rho

     IF ( Medium == 1 ) THEN
        Loc( Medium ) = 0
     ELSE
        Loc( Medium ) = Loc( Medium - 1 ) + NSSPPts( Medium - 1 )
     ENDIF
     ILoc = Loc( Medium )

     !  Read in data and convert attenuation to Nepers/m 
     N1 = 1
     DO I = 1, MaxSSP
        READ(  ENVFile, *    ) z( ILoc + I ), alphaR, betaR, rhoR, alphaI, betaI
        WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
                        z( ILoc + I ), alphaR, betaR, rhoR, alphaI, betaI

        alpha( ILoc + I ) = CRCI( alphaR, alphaI, Freq, AttenUnit )
        beta(  ILoc + I ) = CRCI( betaR,  betaI,  Freq, AttenUnit )
        rho(   ILoc + I ) = rhoR

        ! Did we read the last point?
        IF ( ABS( z( ILoc + I ) - Depth( Medium + 1 ) ) <  EPSILON( 1.0e0 ) * Depth( Medium + 1 ) ) THEN
           NSSPPts( Medium ) = N1
           IF ( Medium == 1 ) Depth( 1 ) = z( 1 )
           RETURN
        ENDIF

        N1 = N1 + 1
     END DO

     ! Fall through means too many points in the profile

     WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
     CALL ERROUT( PRTFile, 'F', 'CLIN', 'Number of SSP points exceeds limit' )

  ELSE   ! Task = 'TABULATE'
     ILoc = Loc( Medium )
     N    = N1 - 1
     H    = ( z( ILoc + NSSPPts( Medium ) ) - z( ILoc + 1 ) ) / N
     Lay  = 1

     DO I = 1, N1
        zT = z( ILoc + 1 ) + ( I - 1 ) * H
        IF ( I == N1 ) zT = z( ILoc + NSSPPts( Medium ) )   ! Make sure no overshoot

        DO WHILE ( zT > z( ILoc + Lay + 1 ) )
           Lay = Lay + 1
        END DO

        R = ( zT - z( ILoc + Lay ) ) / ( z( ILoc + Lay + 1 ) - z( ILoc + Lay ) )
        cP(   I ) = ( 1.0 - R ) * alpha( ILoc + Lay ) + R * alpha( ILoc + Lay+1 )
        cS(   I ) = ( 1.0 - R ) *  beta( ILoc + Lay ) + R *  beta( ILoc + Lay+1 )
        rhoT( I ) = ( 1.0 - R ) *   rho( ILoc + Lay ) + R *   rho( ILoc + Lay+1 )
     END DO
  ENDIF

  RETURN
END SUBROUTINE CLIN
!**********************************************************************!
SUBROUTINE CCUBIC( Depth, cP, cS, rhoT, Medium, N1, Freq, AttenUnit, Task, ENVFile, PRTFile  )

  ! Tabulate cP, cS, rho for specified Medium
  ! using cubic spline interpolation

  IMPLICIT NONE
  INTEGER ENVFile, PRTFile, MaxMedia, MaxSSP
  PARAMETER ( MaxMedia = 501, MaxSSP = 2001 )
  INTEGER              Loc( MaxMedia ), NSSPPts( MaxMedia ), Medium, N, N1, I, ILoc, Lay, IBCBeg, IBCEnd
  REAL     (KIND=8) :: Depth( * ), rhoT( * ), z( MaxSSP ), Freq, alphaR, betaR, rhoR, alphaI, betaI, h, HSPLNE, zT
  COMPLEX  (KIND=8) :: cP( * ), cS( * ), ESPLINE, CRCI, alpha( 4, MaxSSP ), beta( 4, MaxSSP ), rho( 4, MaxSSP )
  CHARACTER (LEN=2) :: AttenUnit
  CHARACTER (LEN=8) :: Task


  SAVE z, alpha, beta, rho, Loc, NSSPPts
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  ! If Task = 'INIT' then this is the first call and SSP is read.
  ! Any other call is a request for SSP subtabulation.

  IF ( Task(1:4) == 'INIT' ) THEN   ! --- Task 'INIT' for initialization
     NSSPPts( Medium ) = N1

     ! The variable Loc(Medium) points to the starting point for the
     ! data in the arrays z, alpha, beta and rho

     IF ( Medium == 1 ) THEN
        Loc( Medium ) = 0
     ELSE
        Loc( Medium ) = Loc( Medium - 1 ) + NSSPPts( Medium - 1 )
     ENDIF
     ILoc = Loc( Medium )

     !  Read in data and convert attenuation to Nepers/m 
     N1 = 1

     DO I = 1, MaxSSP
        READ(  ENVFile, *    ) z( ILoc + I ), alphaR, betaR, rhoR, alphaI, betaI
        WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
                        z( ILoc + I ), alphaR, betaR, rhoR, alphaI, betaI

        alpha(1, ILoc + I) = CRCI( alphaR, alphaI, Freq, AttenUnit )
        beta( 1, ILoc + I) = CRCI( betaR,  betaI,  Freq, AttenUnit )
        rho(  1, ILoc + I) = rhoR

        ! Did we read the last point?
        IF ( ABS( z( ILoc + I ) - Depth( Medium+ 1 ) ) <  EPSILON( 1.0e0 ) * Depth( Medium + 1 ) ) THEN
           NSSPPts( Medium ) = N1
           IF ( Medium == 1 ) Depth( 1 ) = z( 1 )

           !   Compute spline coefs 
           IBCBEG = 0
           IBCEND = 0
           CALL CSPLINE( z( ILoc + 1 ), alpha( 1, ILoc + 1 ), NSSPPts( Medium ), IBCBEG, IBCEND, NSSPPts( Medium ) )
           CALL CSPLINE( z( ILoc + 1 ),  beta( 1, ILoc + 1 ), NSSPPts( Medium ), IBCBEG, IBCEND, NSSPPts( Medium ) )
           CALL CSPLINE( z( ILoc + 1 ),   rho( 1, ILoc + 1 ), NSSPPts( Medium ), IBCBEG, IBCEND, NSSPPts( Medium ) )

           RETURN
        ENDIF

        N1 = N1 + 1
     END DO

     ! Fall through means too many points in the profile

     WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
     CALL ERROUT( PRTFile, 'F', 'CCUBIC', 'Number of SSP points exceeds limit' )

  ELSE   ! Task = 'TABULATE'
     ILoc = Loc( Medium )
     N    = N1 - 1
     H    = ( z( ILoc + NSSPPts( Medium ) ) - z( ILoc + 1 ) ) / N
     Lay  = 1

     DO I = 1, N1
        zT = z( ILoc + 1 ) + ( I - 1 ) * H
        IF ( I == N1 ) zT = z( ILoc + NSSPPts( Medium ) )   ! Make sure no overshoot
        DO WHILE ( zT > z( ILoc + Lay + 1 ) )
           Lay = Lay + 1
        END DO

        HSPLNE = zT - z( ILoc + Lay )

        cP(   I ) =       ESPLINE( alpha( 1, ILoc + Lay ), HSPLNE )
        cS(   I ) =       ESPLINE(  beta( 1, ILoc + Lay ), HSPLNE )
        rhoT( I ) = DBLE( ESPLINE(   rho( 1, ILoc + Lay ), HSPLNE ) )

     END DO
  ENDIF

  RETURN
END SUBROUTINE CCUBIC
