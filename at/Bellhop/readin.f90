SUBROUTINE READIN( FileRoot, Title, freq, iSingle, &
     Nimage, Ibwin, deltas, zBox, rBox, epmult, rLoop, &
     TopOpt, BotOpt, HSTop, HSBot, RunType, BeamType, Component )

  ! Routine to read in and echo all the input data

  ! Note that default values of SSP, DENSITY, Attenuation will not work

  USE anglemod
  USE sspmod
  USE SdRdRMod

  IMPLICIT NONE
  INTEGER, PARAMETER :: ENVFile = 5, PRTFile = 6, RAYFile = 21, ARRFile = 36, SSPFile = 40
  REAL,    PARAMETER :: C0 = 1500.0
  INTEGER               NPts, NMedia, NImage, Isingle, IBWin, AllocateStatus, iostat
  REAL                  atten, ZMin, ZMax
  REAL      (KIND=8) :: x( 2 ), c, gradC( 2 ), crr, crz, czz, freq, deltas, sigma,  &
                        epmult, rLoop, rBox, zBox
  CHARACTER (LEN=80) :: Title, FileRoot
  CHARACTER (LEN=6 ) :: TopOpt
  CHARACTER (LEN=3 ) :: BotOpt, BeamType
  CHARACTER (LEN=5 ) :: RunType
  CHARACTER (LEN=1 ) :: SSPType, AttenUnit, Component
  CHARACTER (LEN=10) :: PlotType

  ! Halfspace properties
   TYPE HSInfo
      CHARACTER (LEN=1) :: BC       ! Boundary condition type
      COMPLEX (KIND=8)  :: cP, cS   ! P-wave, S-wave speeds
      REAL    (KIND=8)  :: rho, Depth  ! density, Depth
      REAL    (KIND=8)  :: BumpDensity, eta, xi   ! Twersky boss parameters
   END TYPE

   TYPE( HSInfo )       :: HSTop, HSBot

  DOUBLE PRECISION :: alphaR, betaR, alphaI, betaI, rhoR
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  ! Open the environmental file
  OPEN( UNIT = ENVFile, FILE = TRIM( FileRoot ) // '.env', STATUS = 'OLD',     IOSTAT = iostat )
  IF ( IOSTAT /= 0 ) THEN   ! successful open?
     WRITE( PRTFile, * ) 'ENVFile = ', TRIM( FileRoot ) // '.env'
     CALL ERROUT( PrtFile, 'F', 'BELLHOP - READIN', 'Unable to open the environmental file' )
  END IF

  ! Open the print file
  OPEN( UNIT = PRTFile, FILE = TRIM( FileRoot ) // '.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )
  IF ( IOSTAT /= 0 ) THEN   ! successful open?
     WRITE( PRTFile, * ) 'PRTFile = ', TRIM( FileRoot ) // '.prt'
     CALL ERROUT( PrtFile, 'F', 'BELLHOP - READIN', 'Unable to open the print file' )
  END IF

  alphaR  = 1500.0;  betaR  = 0.0;  rhoR   = 1.0
  alphaI  = 0.0;     betaI  = 0.0
  iSingle = 0

  Title( 1 : 9 ) = 'BELLHOP- '
  READ(  ENVFile, * ) Title( 10 : 80 )
  WRITE( PRTFile, * ) Title

  READ(  ENVFile, *    ) freq
  WRITE( PRTFile, '('' frequency = '', G11.4, '' Hz'', / )' ) freq

  READ(  ENVFile, * ) NMedia
  WRITE( PRTFile, * ) 'Dummy parameter NMedia = ', NMedia
  IF ( NMedia /= 1 ) CALL ERROUT( PRTFile, 'F', 'READIN', &
       'Only one medium or layer is allowed in BELLHOP; sediment layers must be handled using a reflection coefficient' )

  ! ****** Read in SSP data ******

  TopOpt = '      '   ! initialize to blanks
  READ(  ENVFile, * ) TopOpt
  WRITE( PRTFile, * )

  SSPType   = TopOpt( 1 : 1 )
  HSTop%BC  = TopOpt( 2 : 2 )
  AttenUnit = TopOpt( 3 : 3 )

  ! SSP approximation options

  SELECT CASE ( SSPType )
  CASE ( 'N' )
     WRITE( PRTFile, * ) '    N2-LINEAR approximation to SSP'
  CASE ( 'C' )
     WRITE( PRTFile, * ) '    C-LINEAR approximation to SSP'
  CASE ( 'S' )
     WRITE( PRTFile, * ) '    SPLINE approximation to SSP'
  CASE ( 'Q' )
     WRITE( PRTFile, * ) '    QUAD approximation to SSP'
     OPEN ( FILE = TRIM( FileRoot ) // '.ssp', UNIT = SSPFile, FORM = 'FORMATTED' )
  CASE ( 'A' )
     WRITE( PRTFile, * ) '    ANALYTIC SSP option'
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'READIN', 'Unknown option for SSP approximation' )
  END SELECT

  ! Attenuation options

  SELECT CASE ( AttenUnit )
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

  ! optional addition of volume attenuation using standard formulas

  SELECT CASE ( TopOpt( 4 : 4 ) )
  CASE ( 'T' )
     WRITE( PRTFile, * ) '    THORP attenuation added'
  CASE ( ' ' )
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'READIN', 'Unknown top option letter in fourth position' )
  END SELECT

  SELECT CASE ( TopOpt( 5 : 5 ) )
  CASE ( '*' )
     WRITE( PRTFile, * ) '    Altimetry file selected'
  CASE ( ' ' )
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'READIN', 'Unknown top option letter in fifth position' )
  END SELECT

  SELECT CASE ( TopOpt( 6 : 6 ) )
  CASE ( 'I' )
     WRITE( PRTFile, * ) '    Development options enabled'
  CASE ( ' ' )
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'READIN', 'Unknown top option letter in fifth position' )
  END SELECT

  ! Top BC

  IF ( HSTop%BC == 'A' ) WRITE( PRTFile, FMT = "(//, 'Z    alphaR     betaR      rho     alphaI     betaI', /)" )

  CALL TOPBOT( ENVFile, PRTFile, freq, AttenUnit, HSTop )

  ! Ocean SSP

  READ(  ENVFile, * ) NPts, Sigma, HSBot%Depth
  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Depth = ', HSBot%Depth, 'm'

  IF ( TopOpt( 1 : 1 ) == 'A' ) THEN
     WRITE( PRTFile, * ) 'Analytic SSP option'
  ELSE
     x = [ 0.0, REAL( HSBot%Depth ) ]   ! tells SSP Depth to read to
     CALL SSP( x, C, gradC, crr, crz, czz, SSPType, 'INI' )
  ENDIF

  HSTop%Depth = MINVAL( zSSPV )   ! Depth of top boundary is taken from first SSP point

  ! Bottom BC

  BotOpt = '  '   ! initialize to blanks
  READ(  ENVFile, * ) BotOpt, Sigma
  WRITE( PRTFile, * )
  WRITE( PRTFile, FMT = "(33X, '( RMS roughness = ', G10.3, ' )' )" ) Sigma

  SELECT CASE ( BotOpt( 2 : 2 ) )
  CASE ( '*' )
     WRITE( PRTFile, * ) '    Bathymetry file selected'
  CASE( ' ' )
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'READIN', 'Unknown bottom option letter in second position' )
  END SELECT

  HSBot%BC = BotOpt( 1 : 1 )
  CALL TOPBOT( ENVFile, PRTFile, freq, AttenUnit, HSBot )

  ! *** Read source/receiver Depths (shifting just inside boundaries if outside) ***

  ZMin = SNGL( HSTop%Depth )
  ZMax = SNGL( HSBot%Depth )
  CALL SDRD(   ENVFile, PRTFile, ZMin + 100 * SPACING( ZMin ), &
                                 ZMax - 100 * SPACING( ZMax ) )
  CALL RANGES( ENVFile, PRTFile )  ! *** Read receiver ranges ***

  ! Run type

  READ(  ENVFile, * ) RunType( 1 : 5 )
  WRITE( PRTFile, * )

  SELECT CASE ( RunType( 1 : 1 ) )
  CASE ( 'R' )
     WRITE( PRTFile, * ) 'Ray trace run'
  CASE ( 'E' )
     WRITE( PRTFile, * ) 'Eigenray trace run'
  CASE ( 'I' )
     WRITE( PRTFile, * ) 'Incoherent TL calculation'
  CASE ( 'S' )
     WRITE( PRTFile, * ) 'Semi-coherent TL calculation'
  CASE ( 'C' )
     WRITE( PRTFile, * ) 'Coherent TL calculation'
  CASE ( 'A' )
     WRITE( PRTFile, * ) 'Arrivals calculation, ASCII  file output'
  CASE ( 'a' )
     WRITE( PRTFile, * ) 'Arrivals calculation, binary file output'
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'READIN', 'Unknown RunType selected' )
  END SELECT

  SELECT CASE ( RunType( 2 : 2 ) )
  CASE ( 'C' )
     WRITE( PRTFile, * ) 'Cartesian beams'
  CASE ( 'R' )
     WRITE( PRTFile, * ) 'Ray centered beams'
  CASE ( 'S' )
     WRITE( PRTFile, * ) 'Simple gaussian beams'
  CASE ( 'B' )
     WRITE( PRTFile, * ) 'Geometric gaussian beams'
  CASE DEFAULT
     RunType( 2 : 2 ) = 'G'
     WRITE( PRTFile, * ) 'Geometric hat beams'
  END SELECT

  SELECT CASE ( RunType( 4 : 4 ) )
  CASE ( 'R' )
     WRITE( PRTFile, * ) 'Point source (cylindrical coordinates)'
  CASE ( 'X' )
     WRITE( PRTFile, * ) 'Line source (Cartesian coordinates)'
  CASE DEFAULT
     RunType( 4 : 4 ) = 'R'
     WRITE( PRTFile, * ) 'Point source (cylindrical coordinates)'
  END SELECT

  SELECT CASE ( RunType( 5 : 5 ) )
  CASE ( 'R' )
     WRITE( PRTFile, * ) 'Rectilinear receiver grid: Receivers at ( rr( ir ), rd( ir ) ) )'
  CASE ( 'I' )
     WRITE( PRTFile, * ) 'Irregular grid: Receivers at rr( : ) x rd( : )'
  CASE DEFAULT
     WRITE( PRTFile, * ) 'Rectilinear receiver grid: Receivers at rr( : ) x rd( : )'
     RunType( 5 : 5 ) = 'R'
  END SELECT

  ! Beam angles

  NBeams = 0

  IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
     READ( ENVFile, * ) NBeams, iSingle ! option to trace a single beam
  ELSE
     READ( ENVFile, * ) NBeams
  END IF

  IF ( NBeams == 0 ) THEN   ! automatically estimate NBeams to use
     IF ( RunType( 1 : 1 ) == 'R' ) THEN
        NBeams = 50   ! For a ray trace plot, we don't want too many rays ...
     ELSE
        NBeams = MAX( INT( 0.3 * R( NR ) * Freq / c0 ), 300 )
     END IF
  END IF

  ALLOCATE( alpha( NBeams ), STAT = AllocateStatus )
  IF ( AllocateStatus /= 0 ) THEN
     CALL ERROUT( PRTFile, 'F', 'READIN', 'Insufficient memory to store beam angles'  )
  ENDIF

  IF ( NBeams > 2 ) alpha( 3 ) = -999.9
  READ( ENVFile, * ) alpha
  CALL SUBTABd( alpha, NBeams )
  CALL SORTd(   alpha, NBeams )
  
  ! full 360-degree sweep? remove duplicate beam at -180 degrees
  IF ( alpha( NBeams ) .EQ. alpha( 1 ) + 360.0 ) NBeams = NBeams - 1

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Number of beams   = ', NBeams
  IF ( iSingle > 0 ) WRITE( PRTFile, * ) 'Trace only beam number ', iSingle
  WRITE( PRTFile, * ) 'Beam take-off angles (degrees)'

  IF ( NBeams >= 1 ) WRITE( PRTFile, "( 5G14.6 )" ) alpha( 1 : MIN( NBeams, Number_to_Echo ) )
  IF ( NBeams > 1 .AND. alpha( NBeams ) == alpha( 1 ) ) &
          CALL ERROUT( PRTFile, 'F', 'BELLHOP: READIN', &
          'First and last beam take-off angle are identical' )

  IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
     IF ( iSingle < 1 .OR. iSingle > NBEAMS ) &
          CALL ERROUT( PRTFile, 'F', 'BELLHOP: READIN', &
          'Selected beam, iSingl not in [ 1, NBeams ]' )
  END IF

  ! Limits for tracing beams
  READ(  ENVFile, * ) deltas, zBox, rBox
  rBox = 1000.0 * rBox   ! convert km to m

  ! Automatic step size selection

  IF ( deltas == 0.0 ) THEN
     deltas = ( HSBot%Depth - HSTop%Depth ) / 10.0
  ENDIF

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Step length,     deltas = ', deltas, 'm'
  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Maximum ray Depth, zBox = ', zBox, 'm'
  WRITE( PRTFile, * ) 'Maximum ray range, rBox = ', rBox, 'm'

  ! *** Beam characteristics ***

  ! Curvature change can cause overflow in grazing case
  ! Suppress by setting BeamType(2 : 2) = 'Z'

  IF ( SCAN( RunType( 2 : 2 ), 'GBS' ) /= 0 ) THEN
     BeamType( 1 : 2 ) = 'MS'
     rLoop         = 1.0
     epmult        = 1.0
     BeamType( 3 : 3 ) = RunType( 3 : 3 )
     SELECT CASE ( BeamType( 3 : 3 ) )
     CASE ( 'S' )
        WRITE( PRTFile, * ) 'Beam shift in effect'
     CASE DEFAULT
        WRITE( PRTFile, * ) 'No beam shift in effect'
     END SELECT
  ELSE
     READ(  ENVFile, * ) BeamType( 1 : 2 ), epmult, rLoop
     WRITE( PRTFile, * )
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Type of beam = ', BeamType( 1 : 1 )

     SELECT CASE ( BeamType( 2 : 2 ) )
     CASE ( 'D' )
        WRITE( PRTFile, * ) 'Curvature doubling invoked'
     CASE ( 'Z' )
        WRITE( PRTFile, * ) 'Curvature zeroing invoked'
     CASE ( 'S' )
        WRITE( PRTFile, * ) 'Standard curvature condition'
     CASE DEFAULT
        CALL ERROUT( PRTFile, 'F', 'READIN', 'Unknown curvature condition' )
     END SELECT

     WRITE( PRTFile, * ) 'Epsilon multiplier', epmult
     WRITE( PRTFile, * ) 'Range for choosing beam width', rLoop

     ! Images, windows
     READ(  ENVFile, * ) Nimage, Ibwin, Component
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Number of images, Nimage  = ', Nimage
     WRITE( PRTFile, * ) 'Beam windowing parameter  = ', Ibwin
     WRITE( PRTFile, * ) 'Component                 = ', Component

  ENDIF

  WRITE( PRTFile, * )
  CLOSE( ENVFile )

  ! Write appropriate header information

  SELECT CASE ( RunType( 1 : 1 ) )
  CASE ( 'R', 'E' )   ! Ray trace or Eigenrays
     OPEN ( FILE = TRIM( FileRoot ) // '.ray', UNIT = RAYFile, FORM = 'FORMATTED' )
     WRITE( RAYFile, * ) '''', Title( 1 : 50 ), ''''
     WRITE( RAYFile, * ) freq
     WRITE( RAYFile, * ) NBeams
     WRITE( RAYFile, * ) HSTop%Depth
     WRITE( RAYFile, * ) HSBot%Depth
  CASE ( 'A' )        ! arrival file in ascii format
     OPEN ( FILE = TRIM( FileRoot ) // '.arr', UNIT = ARRFile, FORM = 'FORMATTED' )
     WRITE( ARRFile, * ) freq, NSD, NRD, NR
     WRITE( ARRFile, * ) SD( 1 : NSD )
     WRITE( ARRFile, * ) RD( 1 : NRD )
     WRITE( ARRFile, * ) R(  1 : NR  )
  CASE ( 'a' )        ! arrival file in binary format
     OPEN ( FILE = TRIM( FileRoot ) // '.arr', UNIT = ARRFile, FORM = 'UNFORMATTED' )
     WRITE( ARRFile ) SNGL( freq ), NSD, NRD, NR
     WRITE( ARRFile ) SD( 1 : NSD )
     WRITE( ARRFile ) RD( 1 : NRD )
     WRITE( ARRFile ) R(  1 : NR  )
  CASE DEFAULT
     SELECT CASE ( RunType( 5: 5 ) )
     CASE ( 'R' )
        PlotType = 'rectilin  '
     CASE ( 'I' )
        PlotType = 'irregular '
        IF ( Nrd .NE. Nr ) CALL ERROUT( PRTFile, 'F', 'READIN', 'Irregular grid option selected with Nrd not equal to Nr' )
     END SELECT

     ALLOCATE( theta( 1 ) )
     theta( 1 ) = 0   ! dummy bearing angle
     ntheta     = 1
     atten      = 0.0
     CALL WriteHeader( TRIM( FileRoot ) // '.shd', Title, theta, Ntheta, SD, NSD, RD, NRD, R, NR, &
                       SNGL( freq ), atten, PlotType, 0.0, 0.0 )
  END SELECT

END SUBROUTINE READIN

!**********************************************************************!

SUBROUTINE SSP( x, C, gradC, crr, crz, czz, OPT, Task )

  ! Call the particular profil routine based specified by Task

  ! SSP is expected to perform two tasks:
  !   Task = 'TAB'  then tabulate cp, cs, rhoT 
  !   Task = 'INI' then initialize
  !   X( 2 ) is the (r,z) coordinate

  IMPLICIT NONE
  INTEGER, PARAMETER :: PRTFile = 6
  REAL      (KIND=8) :: x( 2 ), c, gradC( 2 ), crr, crz, czz
  CHARACTER (LEN=1)  :: OPT
  CHARACTER (LEN=3)  :: Task

  SELECT CASE ( OPT )
  CASE ( 'N' )
     CALL N2LIN(  x, C, gradC, crr, crz, czz, Task )
  CASE ( 'C' )
     CALL CLIN(   x, C, gradC, crr, crz, czz, Task )
  CASE ( 'S' )
     CALL CCUBIC( x, C, gradC, crr, crz, czz, Task )
  CASE ( 'Q' )
     CALL QUAD(   x, C, gradC, crr, crz, czz, Task )
  CASE ( 'A' )
     CALL ANALYT( x, C, gradC, crr, crz, czz )
  CASE DEFAULT
     WRITE( PRTFile, * ) 'Non-existent profile option: ', OPT
     STOP
  END SELECT

END SUBROUTINE SSP

!**********************************************************************!

SUBROUTINE N2LIN( x, C, gradC, crr, crz, czz, Task )

  ! N2-linear interpolation of SSP data

  USE sspmod

  IMPLICIT NONE
  INTEGER, PARAMETER :: ENVFile = 5, PRTFile = 6
  CHARACTER (LEN=3)  :: Task
  SAVE ztmp

  INTEGER IZ
  REAL (KIND=8) ::   alphaR, betaR, rhoR, alphaI, betaI, x( 2 ), c, gradC( 2 ), crr, crz, czz, ztmp, Depth, W
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  IF ( Task == 'INI' ) THEN

     ! *** Section to read in SSP data ***

     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'N2-Linear SSP option'
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Sound speed profile:'

     Depth = x( 2 )
     NSSP  = 1

     DO IZ = 1, MaxSSP

        READ(  ENVFile, * ) ztmp, alphaR
        zSSPV( IZ ) = ztmp
        cSSPV( IZ ) = alphaR
        WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2 )" ) zSSPV( IZ ), cSSPV( IZ )

        N2V( IZ ) = 1.0 / cSSPV( IZ ) ** 2

        IF ( IZ > 1 ) N2ZV( IZ - 1 ) = ( N2V(   IZ ) - N2V(   IZ - 1 ) ) / &
             ( zSSPV( IZ ) - zSSPV( IZ - 1 ) )

        ! Have we read the last point?
        IF ( ABS( zSSPV( IZ ) - Depth ) < EPSILON( 1.0e0 ) * Depth ) THEN
           Layer = 1
           RETURN
        ENDIF

        NSSP = NSSP + 1

     END DO

     ! Fall through means too many points in the profile

     WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
     CALL ERROUT( PRTFile, 'F', 'N2LIN', 'Number of SSP points exceeds limit' )

  ELSE

     !  *** Section to return SSP info ***

     IF ( x( 2 ) < zSSPV( Layer ) .OR. x( 2 ) > zSSPV( Layer + 1 ) ) THEN
        DO IZ = 2, NSSP   ! Search for bracketting Depths
           IF ( x( 2 ) < zSSPV( IZ ) ) THEN
              Layer = IZ - 1
              EXIT
           ENDIF
        END DO
     ENDIF

     W = ( x( 2 ) - zSSPV( Layer ) ) / ( zSSPV( Layer + 1 ) - zSSPV( Layer ) )

     C      = 1.0 / SQRT( ( 1.0 - W ) * N2V( Layer ) + W * N2V( Layer + 1 ) )
     gradC  = [ 0.0D0, -0.5 * C * C * C * N2ZV( Layer ) ]

     crr = 0.0d0
     crz = 0.0d0
     czz = 3.0d0 * gradc( 2 ) * gradc( 2 ) / C
  ENDIF

END SUBROUTINE N2LIN

!**********************************************************************!

SUBROUTINE CLIN( x, C, gradC, crr, crz, czz, Task )

  ! C-linear interpolation of SSP data

  USE sspmod

  IMPLICIT NONE
  INTEGER, PARAMETER :: ENVFile = 5, PRTFile = 6
  CHARACTER (LEN=3)  :: Task
  SAVE ztmp

  INTEGER iz
  REAL (KIND=8) ::   alphaR, betaR, rhoR, alphaI, betaI, x( 2 ), c, gradC( 2 ), crr, crz, czz, ztmp, Depth
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  IF ( Task == 'INI' ) THEN

     !  *** Section to read in SSP data ***

     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'C-Linear SSP option'
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Sound speed profile:'

     Depth = x( 2 )
     NSSP = 1

     DO IZ = 1, MaxSSP

        READ(  ENVFile, * ) ztmp, alphaR
        zSSPV( IZ ) = ztmp
        cSSPV( IZ ) = alphaR
        WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2 )" ) zSSPV( IZ ), cSSPV( IZ )

        IF ( IZ > 1 ) czV( IZ - 1 )  = ( cSSPV( IZ ) - cSSPV( IZ - 1 ) ) / &
             ( zSSPV( IZ ) - zSSPV( IZ - 1 ) )

        ! Have we read the last point?
        IF ( ABS( zSSPV( IZ ) - Depth ) < EPSILON( 1.0e0 ) * Depth ) THEN
           Layer = 1
           RETURN
        ENDIF

        NSSP = NSSP + 1

     END DO

     ! Fall through means too many points in the profile

     WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
     CALL ERROUT( PRTFile, 'F', 'N2LIN', 'Number of SSP points exceeds limit' )

  ELSE

     ! *** Section to return SSP info ***

     IF ( x( 2 ) < zSSPV( Layer ) .OR. x( 2 ) > zSSPV( Layer + 1 ) ) THEN
        DO IZ = 2, NSSP   ! Search for bracketting Depths
           IF ( x( 2 ) < zSSPV( IZ ) ) THEN
              Layer = IZ - 1
              EXIT
           ENDIF
        END DO
     ENDIF

     C   = cSSPV( Layer ) + ( x( 2 ) - zSSPV( Layer ) ) * czV( Layer )
     gradC  = [ 0.0D0, czV( Layer ) ]

     crr = 0.0d0
     crz = 0.0d0
     czz = 0.0d0

  ENDIF

END SUBROUTINE CLIN

!**********************************************************************!

SUBROUTINE CCUBIC( x, C, gradC, crr, crz, czz, Task )

  ! Cubic spline interpolation

  USE sspmod

  IMPLICIT NONE
  INTEGER, PARAMETER :: ENVFile = 5, PRTFile = 6

  SAVE ztmp

  INTEGER iz, iBCBeg, iBCEnd
  REAL     (KIND=8) :: alphaR, betaR, rhoR, alphaI, betaI, x( 2 ), &
                       c, gradC( 2 ), cz, crr, crz, czz, ztmp, Depth, HSPLNE
  CHARACTER (LEN=3) :: Task
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  IF ( Task == 'INI' ) THEN

     ! *** Task 'INIT' for initialization ***

     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Spline SSP option'
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Sound speed profile:'

     Depth = x( 2 )
     NSSP  = 1

     DO IZ = 1, MaxSSP

        READ(  ENVFile, * ) ztmp, alphaR
        zSSPV( IZ ) = ztmp
        cSSPV( IZ ) = alphaR
        WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2 )" ) zSSPV( IZ ), cSSPV( IZ )

        CVS( 1, IZ ) = cSSPV( IZ )
        IF ( ABS( zSSPV( IZ ) - Depth ) < EPSILON( 1.0e0 ) * Depth ) THEN
           Layer = 1

           ! Compute spline coefs
           IBCBEG = 0;      IBCEND = 0
           CALL CSPLINE( zSSPV, CVS( 1, 1 ), NSSP, IBCBEG, IBCEND, NSSP )
           RETURN
        ENDIF
        NSSP = NSSP + 1
     END DO

     ! Fall through means too many points in the profile
     WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
     CALL ERROUT( PRTFile, 'F', 'N2LIN', 'Number of SSP points exceeds limit' )
  ELSE

     ! *** Section to return SSP info ***

     IF ( x( 2 ) < zSSPV( Layer ) .OR. x( 2 ) > zSSPV( Layer+1 ) ) THEN
        DO IZ = 2, NSSP   ! Search for bracketting Depths
           IF ( x( 2 ) < zSSPV( IZ ) ) THEN
              Layer = IZ - 1
              EXIT
           ENDIF
        END DO

     ENDIF

     HSPLNE = x( 2 ) - zSSPV( Layer )

     ! C   = SPLINE(   CVS( 1, Layer ), HSPLNE )
     ! cz  = SPLINEX(  CVS( 1, Layer ), HSPLNE )
     ! czz = SPLINEXX( CVS( 1, Layer ), HSPLNE )

     CALL SPLINEALL( CVS( 1, Layer ), HSPLNE, C, cz, czz )

     gradC  = [ 0.0D0, cz ]
     crr = 0.0d0
     crz = 0.0d0
  ENDIF

END SUBROUTINE CCUBIC

!**********************************************************************!

SUBROUTINE QUAD( x, C, gradC, crr, crz, czz, Task )

  ! Bilinear quadrilatteral interpolation of SSP data in 2D

  USE sspmod

  IMPLICIT NONE
  INTEGER, PARAMETER :: ENVFile = 5, PRTFile = 6, SSPFile = 40
  CHARACTER (LEN=3)  :: Task
  SAVE ztmp, cMat, czMat, rSeg, NSeg, ISeg
  INTEGER            :: AllocateStatus, iSeg, iSegT, NSeg, iz, iz2
  REAL (KIND=8), ALLOCATABLE :: cMat( :, : ), czMat( :, : ), rSeg( : )
  REAL      (KIND=8) :: alphaR, betaR, rhoR, alphaI, betaI, x( 2 ), c, gradC( 2 ), &
                        crr, crz, czz, c1, c2, cr, cz, Depth, ztmp, s
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  IF ( Task == 'INI' ) THEN

     !  *** Section to read in SSP data ***

     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Quad SSP option'
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Sound speed profile:'

     Depth = x( 2 )
     NSSP  = 1

     DO IZ = 1, MaxSSP

        READ(  ENVFile, * ) ztmp, alphaR
        zSSPV( IZ ) = ztmp
        cSSPV( IZ ) = alphaR
        WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2 )" ) zSSPV( IZ ), cSSPV( IZ )

        IF ( IZ > 1 ) czV( IZ - 1 )  = ( cSSPV( IZ ) - cSSPV( IZ - 1 ) ) / &
                                       ( zSSPV( IZ ) - zSSPV( IZ - 1 ) )

        ! Have we read the last point?
        IF ( ABS( zSSPV( IZ ) - Depth ) < EPSILON( 1.0e0 ) * Depth ) THEN
           Layer = 1
           iSeg  = 1
           ! Read the 2D SSP matrix

           WRITE( PRTFile, * )
           WRITE( PRTFile, * ) 'Reading sound speed profile from file'

           READ( SSPFile,  * ) NSeg
           WRITE( PRTFile, * ) 'Number of segments = ', NSeg

           IF ( NSeg < 2 ) THEN
              CALL ERROUT( PRTFile, 'F', 'READIN: QUAD', 'You must have a least two profiles in your 2D SSP field'  )
           END IF

           ALLOCATE( cMat( NSSP, NSeg ), czMat( NSSP - 1, NSeg ), rSeg( NSeg ), STAT = AllocateStatus )
           IF ( AllocateStatus /= 0 ) THEN
              CALL ERROUT( PRTFile, 'F', 'READIN: QUAD', 'Insufficient memory to store SSP'  )
           ENDIF

           READ( SSPFile,  * ) rSeg( 1 : NSeg )
           WRITE( PRTFile, * )
           WRITE( PRTFile, * ) 'Profile ranges:'
           WRITE( PRTFile, * ) rSeg( 1 : NSeg )
    
           rSeg = 1000.0 * rSeg   ! convert km to m

           WRITE( PRTFile, * )
           WRITE( PRTFile, * ) 'Sound speed matrix:'
           DO iz2 = 1, NSSP
              READ(  SSPFile, * ) cMat( iz2, : )
              WRITE( PRTFile, * ) 'Layer depth = ', zSSPV( iz2 )
              WRITE( PRTFile, * ) cMat( iz2, : )
           END DO

           CLOSE( SSPFile )

           ! calculate cz
           DO iSegt = 1, NSeg
              DO iz2 = 2, NSSP
                 czMat( iz2 - 1, iSegt ) = ( cMat( iz2, iSegt ) - cMat( iz2 - 1, iSegt ) ) / &
                                           ( zSSPV( iz2 )       - zSSPV( iz2 - 1 ) )
              END DO
           END DO

           RETURN
        ENDIF

        NSSP = NSSP + 1

     END DO

     ! Fall through means too many points in the profile

     WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
     CALL ERROUT( PRTFile, 'F', 'QUAD', 'Number of SSP points exceeds limit' )

  ELSE

     ! *** Section to return SSP info ***

     ! check depth-layer contains x( 2 ) in [ zSSPV( Layer ), zSSPV( Layer + 1 ) ]
     IF ( x( 2 ) < zSSPV( Layer ) .OR. x( 2 ) > zSSPV( Layer + 1 ) ) THEN
        DO IZ = 2, NSSP   ! Search for bracketting Depths
           IF ( x( 2 ) < zSSPV( IZ ) ) THEN
              Layer = IZ - 1
              EXIT
           ENDIF
        END DO
     ENDIF

     ! The following tries to be more efficient than the code above by searching away from the current layer
     ! rather than searching through all the layers
     ! However, seems to be no faster
     ! Also, this code caused a problem on at/tests/Gulf for the range-dep. test cases
!!$     IF ( x( 2 ) < zSSPV( Layer ) .AND. Layer > 1 ) THEN
!!$        DO IZ = Layer - 1, 1, -1   ! Search for bracketting Depths
!!$           IF ( x( 2 ) > zSSPV( IZ ) ) THEN
!!$              Layer = IZ
!!$              EXIT
!!$           ENDIF
!!$        END DO
!!$     ENDIF
!!$
!!$     IF ( x( 2 ) > zSSPV( Layer + 1 ) .AND. Layer < NSSP - 2 ) THEN
!!$        DO IZ = Layer + 2, NSSP   ! Search for bracketting Depths
!!$           IF ( x( 2 ) < zSSPV( IZ ) ) THEN
!!$              Layer = IZ - 1
!!$              EXIT
!!$           ENDIF
!!$        END DO
!!$     ENDIF

     ! check range-segment contains x( 1 ) in [ rSeg( iSeg ), rSeg( iSeg + 1 ) )
     IF ( x( 1 ) < rSeg( iSeg ) .OR. x( 1 ) >= rSeg( iSeg + 1 ) ) THEN
       DO iSegT = 2, NSeg   ! Search for bracketting segment ranges
           IF ( x( 1 ) < rSeg( iSegT ) ) THEN
              iSeg = iSegT - 1
              EXIT
           ENDIF
        END DO
     ENDIF

     ! for this depth, x( 2 ) get the sound speed at both ends of the segment
     c1 = cMat( Layer, iSeg     ) + ( x( 2 ) - zSSPV( Layer ) ) * czMat( Layer, iSeg )
     c2 = cMat( Layer, iSeg + 1 ) + ( x( 2 ) - zSSPV( Layer ) ) * czMat( Layer, iSeg + 1 )

     ! s= proportional distance of x( 1 ) in range
     s = ( x( 1 ) - rSeg( iSeg ) ) / ( rSeg( iSeg + 1 ) - rSeg( iSeg ) )
     s = MIN( s, 1.0D0 )   ! force piecewise constant extrapolation for points outside the box
     s = MAX( s, 0.0D0 )   ! "

     c  = ( 1.0 - s ) * c1 + s * c2
     cr = ( c2 - c1 ) / ( rSeg( iSeg + 1 ) - rSeg( iSeg ) )
     cz = ( 1.0 - s ) * czMat( Layer, iSeg ) + s * czMat( Layer, iSeg + 1 )

     gradC = [ cr, cz ]

     crr = 0.0
     crz = 0.0
     czz = 0.0

  ENDIF

END SUBROUTINE QUAD

!**********************************************************************!

SUBROUTINE TOPBOT( ENVFile, PRTFile, freq, AttenUnit, HS )

  ! Handles top and bottom boundary conditions

  ! Input:
  !   ENVFile: Environmental file
  !   PRTFile: Print file
  !   freq:   frequency
  !   BCType: Boundary condition type

  ! Output:
  !  cpHS:    P-wave speed in halfspace
  !  csHS:    S-wave speed in halfspace
  !  rhoHS:   density in halfspace

  !  BumDenN:  Bump density
  !  eta:     Principal radius 1
  !  xi:      Principal radius 2

  IMPLICIT NONE

  INTEGER           :: ENVFile, PRTFile
  REAL     (KIND=8) :: alphaR, betaR, alphaI, betaI, rhoR, freq, ztmp
  COMPLEX  (KIND=8) :: CRCI
  CHARACTER (LEN=1) :: AttenUnit

   ! Halfspace properties
   TYPE HSInfo
      CHARACTER (LEN=1) :: BC       ! Boundary condition type
      COMPLEX (KIND=8)  :: cP, cS   ! P-wave, S-wave speeds
      REAL    (KIND=8)  :: rho, Depth  ! density, depth
      REAL    (KIND=8)  :: BumpDensity, eta, xi   ! Twersky boss parameters
   END TYPE

   TYPE( HSInfo )       :: HS

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
     CALL ERROUT( PRTFile, 'F', 'TOPBOT', 'Unknown boundary condition type' )
  END SELECT

  ! ****** Read in BC parameters depending on particular choice ******

  HS%cp  = 0.0
  HS%cs  = 0.0
  HS%rho = 0.0

  SELECT CASE ( HS%BC )
  CASE ( 'S', 'H', 'T', 'I' )   ! *** Twersky ice model parameters ***
     READ(  ENVFile, *    ) HS%BumpDensity, HS%eta, HS%xi
     WRITE( PRTFile, FMT = "( /, ' Twersky ice model parameters:' )" )
     WRITE( PRTFile, FMT = "(' Bumden = ', G15.6, '  Eta = ', G11.3, '  Xi = ', G11.3, /)" ) &
          HS%BumpDensity, HS%eta, HS%xi
  CASE ( 'A' )                  ! *** Half-space properties ***
     READ(  ENVFile, *    ) ztmp, alphaR, betaR, rhoR, alphaI, betaI
     WRITE( PRTFile, FMT = "( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
          ztmp, alphaR, betaR, rhoR, alphaI, betaI

     HS%cp  = CRCI( alphaR, alphaI, freq, AttenUnit )
     HS%cs  = CRCI( betaR,  betaI,  freq, AttenUnit )
     HS%rho = rhoR
  END SELECT

END SUBROUTINE TOPBOT

!**********************************************************************!

FUNCTION CRCI( c, alpha, freq, AttenUnit )

  ! Converts real wave speed and attenuation to a single complex wave speed
  ! 6 CASES:
  !     N for Nepers/meter
  !     M for dB/meter      (M for Meters)
  !     F for dB/m-kHZ      (F for frequency dependent)
  !     W for dB/wavelength (W for Wavelength)
  !     Q for Q
  !     L for Loss parameter
  !     T for Thorpe

  IMPLICIT NONE
  REAL     (KIND=8) :: omega, freq, alpha, alphaT, c, f2
  COMPLEX  (KIND=8) :: CRCI
  CHARACTER (LEN=1) :: AttenUnit

  omega = 2.0 * 3.1415926535897932 * freq

  ! Convert to Nepers/m

  alphaT = 0.0

  SELECT CASE ( AttenUnit )
  CASE ( 'N' )   ! Nepers/m
     alphaT = alpha
  CASE ( 'M' )   ! dB/meter
     alphaT = alpha / 8.6858896
  CASE ( 'F' )   ! dB/m-kHZ
     alphaT = alpha * freq / 8685.8896
  CASE ( 'W' )   ! dB/wavelength
     IF ( c /= 0.0        ) alphaT = alpha * freq / ( 8.6858896 * c )
  CASE ( 'Q' )
     IF( c * alpha /= 0.0 ) alphaT = omega / ( 2.0 * c * alpha )
  CASE ( 'L' )   ! loss parameter
     IF ( c /= 0.0        ) alphaT = alpha * omega / c
  CASE ( 'T' )   ! Thorpe
     f2     = ( freq / 1000.0 ) **2

     ! Original formula from Thorp 1967
     ! alphaT = 40.0 * f2 / ( 4100.0 + f2 ) + 0.1 * f2 / ( 1.0 + f2 )   ! dB/kyard
     ! alphaT = alphaT / 914.4D0     ! dB / m
     ! alphaT = alphaT / 8.6858896D0 ! Nepers / m

     ! Updated formula from JKPS Eq. 1.34
     alphaT = 3.3d-3 + 0.11 * f2 / ( 1.0 + f2 ) + 44.0 * f2 / ( 4100.0 + f2 ) + 3d-4* f2   ! dB/km
     alphaT = alpha / 8685.8896D0 ! Nepers / m

  END SELECT

  ! Convert Nepers/m to equivalent imaginary sound speed

  alphaT = alphaT * c**2 / omega
  CRCI   = CMPLX( c, alphaT )

END FUNCTION CRCI
