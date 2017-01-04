PROGRAM KRAKENC

  ! Ocean acoustic normal modes.

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

  ! Originally developed as part of the author's dissertation under the supervision
  ! of Prof. Edward L. Reiss, Northwestern University

  USE krakcmod
  USE SdRdRMod
  USE RefCoMod
  IMPLICIT NONE
  INTEGER            :: Min_Loc( 1 ), IFirst, ILast, IRec
  REAL               :: zMin, zMax
  REAL      (KIND=8) :: omega, Error
  CHARACTER (LEN=80) :: FileRoot

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  DO iProf = 1, 9999   ! Loop over a sequence of profiles

     NV( 1 : 5 ) = [ 1, 2, 4, 8, 16 ]

     ! Read in environmental info
     Title = 'KRAKENC-'
     CALL READIN( FileRoot, Title, Freq, MaxMedium, NMedia, TopOpt, HSTop, NG, sigma, Depth, BotOpt, HSBot, ENVFile, PRTFile )

     READ(  ENVFile, * ) cLow, cHigh   ! Spectral limits
     WRITE( PRTFile, "( /, ' cLow = ', G10.5, 'm/s      cHigh = ', G10.5, 'm/s' )" ) cLow, cHigh

     READ(  ENVFile, * ) RMax          ! Maximum range for calculations
     WRITE( PRTFile, * ) 'RMax = ', RMax

     ! Read source/receiver depths
     zMin = SNGL( Depth( 1 ) )
     zMax = SNGL( Depth( NMedia + 1 ) )
     CALL SDRD( ENVFile, PRTFile, zMin, zMax )

     omega2 = ( 2.0D0 * pi * Freq ) ** 2
     M      = MaxM

     ! Optionally read in Bot, Top reflection coefficients
     CALL READRC( FileRoot, HSBot%BC, HSTop%BC, PRTFile )

     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Mesh multiplier   CPU seconds'

     DO ISet = 1, NSets   ! Main loop: solve the problem for a sequence of meshes
        N( 1 : NMedia ) = NG( 1 : NMedia ) * NV( ISet )
        h( 1 : NMedia ) = ( Depth( 2 : NMedia + 1 ) - Depth( 1 : NMedia ) ) / N( 1 : NMedia )
        hV( ISet )      = h( 1 )
        CALL SOLVE( FileRoot, Error )
        IF ( Error * 1000.0 * RMax < 1.0 ) GOTO 3000
     END DO

     ! Fall through indicates failure to converge
     CALL ERROUT( PRTFile, 'W', 'KRAKENC', 'Too many meshes needed: check convergence' )

3000 omega   = SQRT( omega2 )   ! Solution complete: discard modes with phase velocity above cHigh
     Min_Loc = MinLOC( DBLE( Extrap( 1, 1 : M ) ), DBLE( Extrap( 1, 1 : M ) ) > omega2 / cHigh ** 2 )
     M       = Min_Loc( 1 )

     ! Write eigenvalues to PRTFile and MODFile
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) '   I          k             alpha          Phase Speed       Group Speed'

     ! k() contains scatter losses; Extrap contains extrapolated values of wavenumbers, so combine ...
     k( 1 : M ) = SQRT( Extrap( 1, 1 : M ) + k( 1 : M ) )

     DO mode = 1, M
        WRITE( PRTFile, '( I5, 4G18.10 )' ) mode, k( mode ), omega / DBLE( k( mode ) ), VG( mode )
        ! Zero out positive imaginary part which would cause growth in range. Should be small.
        IF ( AIMAG( k( mode ) ) > 0.0D0 ) k( mode ) = REAL( k( mode ) )
     END DO

     WRITE( MODFile, REC = IRecProfile + 4 ) M, LRecordLength

     IFirst = 1
     DO IREC = 1, 1 + ( 2 * M - 1 ) / LRecordLength
        ILast  = MIN( M, IFirst + LRecordLength / 2 - 1 )
        WRITE( MODFile, REC = IRecProfile + 5 + M + IREC ) CMPLX( k( IFirst : ILast ) )
        IFirst = ILast + 1
     END DO

     ! set record pointer to beginning of next mode set
     IRecProfile = IRecProfile + 6 + M + 1 + ( 2 * M - 1 ) / LRecordLength

  END DO   ! next profile
  CLOSE( ENVFile )
  CLOSE( MODFile )

END PROGRAM KRAKENC

!**********************************************************************

SUBROUTINE INIT

  ! Initializes arrays defining difference equations

  USE krakcmod
  IMPLICIT NONE
  LOGICAL           :: ElasticFlag = .FALSE.
  INTEGER           :: IAllocStat, ii, J, Medium, NPoints, N1
  REAL     (KIND=8) :: Two_h
  COMPLEX  (KIND=8) :: cP2, cS2
  COMPLEX (KIND=8), ALLOCATABLE :: cP( : ), cS( : )
  CHARACTER (LEN=8) :: Task

  cMin          = HUGE( cMin )
  FirstAcoustic = 0
  Loc( 1 )      = 0

  ! Allocate storage for finite-difference coefficients

  NPoints = SUM( N( 1 : NMedia ) ) + NMedia

  IF ( ALLOCATED( B1 ) ) DEALLOCATE( B1, B2, B3, B4, rho )
  ALLOCATE ( B1( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints ), rho( NPoints ), &
             cP( NPoints ), cS( NPoints ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) &
       CALL ERROUT( PRTFile, 'F', 'KRAKENC - INIT', 'Insufficient memory to allocate B1, B2, B3, B4 vectors. Reduce mesh.' )

  DO Medium = 1, NMedia   ! Loop over media

     IF ( Medium /= 1 ) Loc( Medium ) = Loc( Medium-1 ) + N( Medium-1 ) + 1
     N1  = N(   Medium ) + 1
     ii  = Loc( Medium ) + 1

     ! PROFIL reads in the data for a medium
     Task = 'TAB'
     CALL PROFIL( Depth, cP( ii ), cS( ii ), rho( ii ), Medium, N1, Freq, TopOpt( 1 : 1 ), TopOpt( 3 : 4 ), Task, ENVFile, PRTFile )

     ! Load diagonals of the finite-difference equations

     IF ( cS( ii ) == ( 0.0, 0.0 ) ) THEN  ! Case of an acoustic medium

        Material( Medium ) = 'ACOUSTIC'
        IF ( FirstAcoustic == 0 ) FirstAcoustic = Medium
        LastAcoustic = Medium

        cMin = MIN( cMin, MINVAL( DBLE( cP( ii : ii + N( Medium ) ) ) ) )

        B1( ii : ii + N( Medium ) ) = -2.0D0 + h( Medium ) ** 2 * omega2 / cP( ii : ii + N( Medium ) ) ** 2

     ELSE                                  ! Case of an elastic medium
        IF ( sigma( Medium ) /= 0.0D0 ) &
           CALL ERROUT( PRTFile, 'F', 'KRAKENC', 'Rough elastic interfaces are not allowed' )

        Material( Medium ) = 'ELASTIC'
        ElasticFlag        = .TRUE.
        Two_h              = 2.0D0 * h( Medium )

        DO J = ii, ii + N( Medium )
           cMin = MIN( DBLE( cS( J ) ), cMin )

           cP2 = cP( J ) ** 2
           cS2 = cS( J ) ** 2

           B1(  J ) = Two_h / ( rho( J ) * cS2 )
           B2(  J ) = Two_h / ( rho( J ) * cP2 )
           B3(  J ) = 4.0D0 * Two_h * rho( J ) * cS2 * ( cP2 - cS2 ) / cP2
           B4(  J ) = Two_h * ( cP2 - 2.0D0 * cS2 ) / cP2
           rho( J ) = Two_h * omega2 * rho( J )
        END DO

     ENDIF

  END DO   ! next Medium

  ! Bottom properties
  IF ( HSBot%BC == 'A' ) THEN
     IF ( HSBot%cS /= ( 0.0, 0.0 ) ) THEN    ! Elastic  bottom half-space
        ElasticFlag = .TRUE.
        cMin   = MIN( cMin, DBLE( HSBot%cS ) )
     ELSE                                    ! Acoustic bottom half-space
        cMin   = MIN( cMin, DBLE( HSBot%cP ) )
     ENDIF
  ENDIF

  ! Top properties
  IF ( HSTop%BC == 'A' ) THEN
     IF ( HSTop%cS /= ( 0.0, 0.0 ) ) THEN    ! Elastic  top half-space
        ElasticFlag = .TRUE.
        cMin   = MIN( cMin, DBLE( HSTop%cS ) )
     ELSE                                    ! Acoustic top half-space
        cMin   = MIN( cMin, DBLE( HSTop%cP ) )
     ENDIF
  ENDIF

  ! If elastic medium then reduce cMin for Scholte wave
  IF ( ElasticFlag ) cMin = 0.85 * cMin
  cLow = MAX( cLow, 0.99 * cMin )

END SUBROUTINE INIT

!**********************************************************************

SUBROUTINE SOLVE( FileRoot, Error )

  ! Solves the eigenvalue problem at the current mesh and produces a new extrapolation

  USE krakcmod
  IMPLICIT NONE
  INTEGER                           :: Min_Loc( 1 ), J, Key
  REAL    (KIND=8)                  :: x1, x2, TStart, TEnd
  COMPLEX (KIND=8)                  :: T1, T2, F1, F2, xTemp( MaxM )
  REAL    (KIND=8),   INTENT( OUT ) :: Error
  CHARACTER (LEN=80), INTENT(  IN ) :: FileRoot


  CALL CPU_TIME( Tstart )
  CALL INIT                 ! set up the finite-difference mesh
  CALL SOLVE2( FileRoot )   ! solve for the eigenvalues

  Extrap( ISet, 1 : M ) = EVMat( ISet, 1 : M )

  ! Remove any eigenvalues outside the spectral limits
  ! Typically the last 'eigenvalue' results from forcing a zero in funct( x ) for x outside the limits
  ! Inverse iteration would usually fail for that mode

  Min_Loc = MINLOC( REAL( Extrap( 1, 1 : M ) ), REAL( Extrap( 1, 1 : M ) ) > omega2 / cHigh ** 2 )
  M       = Min_Loc( 1 )

  IF ( ISet == 1 ) CALL VECTOR( FileRoot )   ! If this is the first mesh, compute the eigenvectors

  xTemp( 1 : M ) = Extrap( Iset, 1 : M )
  CALL ORDER( xTemp, M )    ! order eigenvalues by real part
  EVMat(  ISet, 1 : M ) = xTemp( 1 : M )
  Extrap( ISet, 1 : M ) = xTemp( 1 : M )

  ! Richardson extrapolation to improve the accuracy

  Error = 1.0D10          ! initialize error to a large number
  KEY   = 2 * M / 3 + 1   ! index of element used to check convergence

  IF ( ISet > 1 ) THEN
     T1 = Extrap( 1, KEY )

     DO J = ISet - 1, 1, -1
        DO mode = 1, M
           x1 = NV( J    ) ** 2
           x2 = NV( ISet ) ** 2
           F1 = Extrap( J,     mode )
           F2 = Extrap( J + 1, mode )
           Extrap( J, mode ) = F2 - ( F1 - F2 ) / ( x2 / x1 - 1.0D0 )
        END DO
     END DO

     T2    = Extrap( 1, KEY )
     Error = ABS( T2 - T1 )
  ENDIF

  CALL CPU_TIME( Tend )   ! check elapsed time
  ET( ISet ) = Tend - Tstart
  WRITE( PRTFile, '( 1X, I8, 6X, G15.3, ''s'' )' ) NV( ISet ), ET( ISet )

END SUBROUTINE SOLVE

!**********************************************************************

SUBROUTINE SOLVE2( FileRoot )

  ! Provides initial guess to root finder for each EVMat( I )

  ! A sequence of restarts is done to compute a complete set of modes for the first two sets (ISet = 1, 2).
  ! 'Complete' in this case means that all the modes in the user-specified interval are found.
  ! However, if, for ISet=2, the same number of modes is found as for ISet=1 then it is assumed
  ! that the set is complete and no restarts are performed

  USE krakcmod
  IMPLICIT NONE
  INTEGER            :: ITry, Iteration, MaxIteration, J, MaxTries, NzTab, IFirst, ILast, IRec, ii
  REAL      (KIND=4) :: RVAR1, RVAR2, kLow, kHigh
  REAL      (KIND=8) :: Tolerance, x1, x2
  COMPLEX            :: kIG( MaxM )
  COMPLEX   (KIND=8) :: x, P( 10 ), xTemp( MaxM ), cTry, kTry
  CHARACTER (LEN=80) :: ErrorMessage, FileRoot

  x            = ( 1.0D0, 0.0D0 ) * omega2 / cLow ** 2
  MaxIteration = 1000   ! maximum # of iterations in root-finder

  IF ( TopOpt( 5 : 5 ) == '.' .AND. ISet <= 2 ) THEN
     MaxTries = 5000   ! MaxTries = # of restarts for the root finder
  ELSE
     MaxTries = 1
  ENDIF

  WRITE( PrtFile, * ) 'Max. number of restarts in root finder, MaxTries = ', MaxTries

  ITry = 1   ! counter that keeps track of # of tries

  ! optionally read an existing mode file to get initial guesses
  IF ( ISet == 1 .AND. TopOpt( 4 : 4 ) == 'I' ) THEN
     OPEN ( FILE = 'MODFile', UNIT = MODFile, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED' )

     READ( MODFile, REC = 5 ) M, LRecordLength

     IFirst = 1
     DO IREC = 1, 1 + ( 2 * M - 1 ) / LRecordLength
        ILast  = MIN( M, IFirst + LRecordLength / 2 - 1 )
        READ( MODFile, REC = 6 + M + IREC ) ( kIG( mode ), mode = IFirst, ILast )
        IFirst = ILast + 1
     END DO
  END IF

  ! start looking for each root in succession
  DO mode = 1, M
     ! For first or second meshes, use a high guess
     ! Otherwise use extrapolation to produce an initial guess

5500 CONTINUE
     x = 1.00001 * x

     IF ( ISet == 1 .AND. TopOpt( 4 : 4 ) == 'I' ) THEN
        x = kIG( mode ) ** 2   ! initial guess from MODFile
     ELSE IF ( ISet >= 2 ) THEN
        P( 1 : ISet - 1 ) = EVMat( 1 : ISet - 1, mode )

        IF ( ISet >= 3 ) THEN
           DO ii = 1, ISet - 2
              DO J = 1, ISet - ii - 1
                 x1 = hV( J      ) ** 2
                 x2 = hV( J + ii ) ** 2

                 P( J ) = ( ( hV( ISet ) ** 2 - x2 ) * P( J     ) - &
                            ( hV( ISet ) ** 2 - x1 ) * P( J + 1 ) ) &
                          / ( x1 - x2 )
              END DO
           END DO
           x = P( 1 )
        ENDIF

     ENDIF

     ! Use the secant method to refine the eigenvalue
     Tolerance = ABS( x ) * 10.0D0 ** ( 4.0 - PRECISION( x ) )
     CALL ZSECCX( x, Tolerance, Iteration, MaxIteration, ErrorMessage )

     IF ( ErrorMessage /= ' ' ) THEN   ! any problems in root finder?
        ! CALL ERROUT( PRTFile, 'W', 'KRAKENC-ZSECCX', ErrorMessage )
        x = TINY( REAL( x ) )   ! make sure value discarded later
     ENDIF

     ! Is the mode inside the user specified spectral limits?

     IF ( SQRT( omega2 ) / cHigh > REAL( SQRT( x ) ) ) THEN
        ! Mode outside, restart at a random point
        CALL RANDOM_NUMBER( RVAR1 )
        CALL RANDOM_NUMBER( RVAR2 )

        !cTry  = cLow + RVAR1 * ( cHigh - cLow ) + ( 0.0, 0.01 ) * RVAR2 * cLow
        !x     = omega2 / cTry ** 2

        kLow  = sqrt( omega2 ) / cHigh
        kHigh = sqrt( omega2 ) / cLow
        kTry  =  kLow + RVAR1 * ( kHigh - kLow ) - ( 0, 0.01D0 ) * RVAR2 * kLow
        x     = kTry ** 2
        cTry  = sqrt( omega2 ) / kTry

        MaxIteration = 30 ! but don't let root-finder search long

        IF ( ITry < MaxTries ) THEN
           ! MaxTries = MAX( MaxTries, 5 * mode )
           ITry     = ITry + 1
           ! WRITE( *, * ) 'Restart at phase speed = ', ITry, cTry
           GOTO 5500
        ELSE
           M = mode - 1
           EXIT   ! done searching for modes
        ENDIF
     ENDIF

     EVMat( ISet, mode ) = x
     xTemp( mode )       = x

  END DO   ! next mode

  ! If no modes, open a dummy MODFile for use by FIELD3D
  IF ( M == 0 ) THEN
     LRecordLength = 32
     NzTab         = 0

     ! Open MODFile and write header
     IF ( iProf == 1 ) THEN
        OPEN ( FILE = TRIM( FileRoot) //'.mod', UNIT = MODFile, ACCESS = 'DIRECT', RECL = 4 * LRecordLength, FORM = 'UNFORMATTED' )
     END IF

     WRITE( MODFile, REC = 1 ) LRecordLength, Title, REAL( Freq ), 1, NzTab, NzTab
     WRITE( MODFile, REC = 5 ) M, LRecordLength
     CALL ERROUT( PRTFile, 'F', 'KRAKENC', 'No modes for given phase speed interval' )
  ENDIF

  IF ( M == MaxM ) THEN   ! Have we hit the ceiling on max # modes
     WRITE( PRTFile, * ) 'Number of modes = ', M
     WRITE( PRTFile, * ) 'Number of modes allowable = ', MaxM
     CALL ERROUT( PRTFile, 'W', 'KRAKENC', 'Too many modes: Program will compute as many as it can' )
  ENDIF

  CALL ORDER( xTemp, M )    ! order eigenvalues by real part
  EVMat( ISet, 1 : M ) = xTemp( 1 : M )

END SUBROUTINE SOLVE2

!**********************************************************************

SUBROUTINE FUNCT( x, Delta, IPower )

  ! FUNCT( x ) = 0 is the dispersion relation

  USE krakcmod
  IMPLICIT NONE
  INTEGER,           PARAMETER :: IPowerR = 50, IPowerF = -50
  REAL (KIND=8),     PARAMETER :: Roof = 1.0D+50, Floor = 1.0D-50
  INTEGER                      :: IPower, IPower1, J
  COMPLEX (KIND=8), INTENT(IN) :: x
  COMPLEX (KIND=8)             :: Delta, F, G, F1, G1

!!$  IF ( REAL( x ) <= omega2 / cHigh ** 2 ) THEN    ! For a k below the spectrum limit, force a zero
!!$     Delta  = 0.0D0
!!$     IPower = 0
!!$     RETURN
!!$  ENDIF

  CALL BCImpedance(    x, 'BOT', HSBot, F,  G,  IPower  ) ! Bottom impedance
  CALL AcousticLayers( x,               F,  G,  IPower  ) ! Shoot through acoustic layers
  CALL BCImpedance(    x, 'TOP', HSTop, F1, G1, IPower1 ) ! Top impedance

  Delta  = F * G1 - G * F1
  IPower = IPower + IPower1

  ! Delta = f /g - f1 / g1
  ! Delta = g / f
  ! Delta = 1 / ( log( f * g1 - g * f1 ) + log( 10.0d0 ) * ( IPower + IPower1 ) )
  ! IPower = 0

  ! Deflate previous roots

  IF ( mode > 1 ) THEN
     DO J = 1, mode - 1
        Delta = Delta / ( x - EVMat( ISet, J ) )

        ! Scale if necessary
        DO WHILE ( ABS( DBLE( Delta ) ) < Floor .AND. ABS( Delta ) > 0.0D0 )
           Delta  = Roof * Delta
           IPower = IPower - IPowerR
        END DO

        DO WHILE ( ABS( DBLE( Delta ) ) > Roof )
           Delta  = Floor * Delta
           IPower = IPower - IPowerF
        END DO

     END DO
  ENDIF

END SUBROUTINE FUNCT

!**********************************************************************

SUBROUTINE AcousticLayers( x, F, G, IPower )

  ! Shoot through acoustic layers

  USE krakcmod
  IMPLICIT NONE
  INTEGER,           PARAMETER :: IPowerF = -50
  REAL (KIND=8),     PARAMETER :: Roof = 1.0D+50, Floor = 1.0D-50
  INTEGER                      :: IPower, ii, Medium
  REAL (KIND=8)                :: rhoMedium
  COMPLEX (KIND=8), INTENT(IN) :: x
  COMPLEX (KIND=8)             :: F, G, p0, p1, p2, h2k2
  ! for 2-layer code at end: COMPLEX (KIND=8) :: gamma1, gamma2

  IF ( FirstAcoustic == 0 ) RETURN

  DO Medium = LastAcoustic, FirstAcoustic, -1   ! Loop over successive acoustic media

     h2k2      = h( Medium ) ** 2 * x
     ii        = Loc( Medium ) + N( Medium ) + 1
     rhoMedium = rho(  Loc( Medium ) + 1  )   ! density is made homogeneous using value at top of each medium

     p1   = -2.0D0 * G
     p2   = ( B1( ii ) - h2k2 ) * G - 2.0D0 * h( Medium ) * F * rhoMedium

     ! Shoot (towards surface) through a single medium
     DO ii = Loc( Medium ) + N( Medium ), Loc( Medium ) + 1, -1

        p0 = p1
        p1 = p2
        p2 = ( h2k2 - B1( ii ) ) * p1 - p0

        DO WHILE ( ABS(  DBLE( p2 ) ) > Roof )  ! Scale if necessary
           p0     = Floor * p0
           p1     = Floor * p1
           p2     = Floor * p2
           IPower = IPower - IPowerF
        END DO

     END DO

     ! F = P' / rho and G = -P since F P + G P' / rho = 0
     rhoMedium = rho( Loc( Medium ) + 1 )   ! density at top of layer
     F         = -( p2 - p0 ) / ( 2.0D0 * h( Medium ) ) / rhoMedium
     G         = -p1
  END DO

  ! here's a little two-layer analytic formula as a test
  !d1 =  5000
  !d2 = 10000
  !gamma1 = SQRT( ( 2 + b1(    1 ) ) / h( 1 ) ** 2 - x )
  !gamma2 = SQRT( ( 2 + b1( 6000 ) ) / h( 2 ) ** 2 - x )
  !rho2 = 1.8
  !f = 0
  !g = SIN( gamma1 * d1 ) * COS( gamma2 * (d1-d2) ) / gamma1 - COS( gamma1 * d1 ) * SIN( gamma2 * ( d1-d2 ) ) / gamma2
  !IPower = 0

END SUBROUTINE AcousticLayers

!**********************************************************************

SUBROUTINE VECTOR( FileRoot )

  ! Do inverse iteration to compute each of the eigenvectors and write these to the disk file

  USE krakcmod
  USE SdRdRMod
  IMPLICIT NONE
  INTEGER              :: IAllocStat, IErr, IPower, ii, ITP, J, jj, L, Medium, NTotal, NTotal1, NzTab
  INTEGER, ALLOCATABLE :: IzTab( : )
  REAL                 :: zTab( NSD + NRD )
  REAL    (KIND=8)     :: h_rho
  COMPLEX, ALLOCATABLE :: PhiTab( : )
  COMPLEX (KIND=8)     :: xh2, x, F, G
  REAL, ALLOCATABLE    :: z( : ), WTS( : )
  COMPLEX (KIND=8), ALLOCATABLE :: RV1( : ), RV2( : ), RV3( : ), RV4( : ), Phi( : ), d( : ), e( : )
  CHARACTER (LEN=80)   :: FileRoot

  ! Tabulate z-coordinates and off-diagonals of matrix
  NTotal  = SUM( N( FirstAcoustic : LastAcoustic ) )
  NTotal1 = NTotal + 1

  ALLOCATE( z( NTotal1 ), e( NTotal1 + 1 ), d( NTotal1 ), Phi( NTotal1 ), &
          RV1( NTotal1 ), RV2( NTotal1 ), RV3( NTotal1 ), RV4( NTotal1 ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'KRAKENC - VECTOR', 'Insufficient memory: Reduce mesh.' )

  J      = 1
  z( 1 ) = SNGL( Depth( FirstAcoustic ) )

  DO Medium = FirstAcoustic, LastAcoustic

     h_rho = h( Medium ) * rho( Loc( Medium ) + 1 )        ! density at the top of each layer

     e( J + 1 : J + N( Medium ) ) = 1.0D0 / h_rho
     z( J + 1 : J + N( Medium ) ) = z( J ) + SNGL( h( Medium ) * [ (jj, jj = 1, N( Medium ) ) ] )

     J = J + N( Medium )
  END DO

  e( NTotal1 + 1 ) = 1.0D0 / h_rho       ! Dummy value; never used

  ! Calculate the indices, weights, ... for mode interpolation
  CALL MERGEV( SD, NSD, RD, NRD, zTab, NzTab )
  ALLOCATE( WTS( NzTab ), IzTab( NzTab ), PhiTab( NzTab ) )
  CALL WEIGHT( z, NTotal1, zTab, NzTab, WTS, IzTab )

  ! Open MODFile and write header

  IF ( iProf == 1 ) THEN
     ! LRecordLength must not increase between profiles !!!
     LRecordLength = MAX( 2 * NzTab, 32, 3 * ( LastAcoustic - FirstAcoustic + 1 ) )   ! Logical record length in `longwords' (4 bytes)
     OPEN ( FILE = TRIM( FileRoot) //'.mod', UNIT = MODFile, ACCESS = 'DIRECT', RECL = 4 * LRecordLength, FORM = 'UNFORMATTED' )
  END IF

  WRITE( MODFile, REC = IRecProfile     ) LRecordLength, Title, REAL( Freq ), LastAcoustic - FirstAcoustic + 1, NzTab, NzTab
  WRITE( MODFile, REC = IRecProfile + 1 ) ( N( Medium ), Material( Medium ), Medium = FirstAcoustic, LastAcoustic )
  WRITE( MODFile, REC = IRecProfile + 2 ) HSTop%BC, CMPLX( HSTop%cP ), CMPLX( HSTop%cS ), REAL( HSTop%rho ), &
                                          REAL( Depth( 1          ) ), &
                                          HSBot%BC, CMPLX( HSBot%cP ), CMPLX( HSBot%cS ), REAL( HSBot%rho ), &
                                          REAL( Depth( NMedia + 1 ) )
  WRITE( MODFile, REC = IRecProfile + 3 ) ( REAL( Depth( Medium ) ), REAL( rho( Loc( Medium ) + 1 ) ), &
                                          Medium = FirstAcoustic, LastAcoustic )
  WRITE( MODFile, REC = IRecProfile + 5 ) zTab( 1 : NzTab )

  ! Main loop: for each eigenvalue call SINVIT to get eigenvector

  DO mode = 1, M
     x = EVMat( 1, mode )

     ! Corner elt requires top impedance
     CALL BCImpedance( x, 'TOP', HSTop, F, G, IPower )

     IF ( G == 0.0D0 ) THEN
        d( 1 ) = 1.0D0
        e( 2 ) = 0.0D0
     ELSE
        L      = Loc( FirstAcoustic ) + 1
        xh2    = x * h( FirstAcoustic ) * h( FirstAcoustic )
        h_rho  = h( FirstAcoustic ) * rho( L )
        d( 1 ) = ( B1( L ) - xh2 ) / h_rho / 2.0D0 + F / G
     ENDIF

     ! Set up the diagonal
     ITP = NTotal
     J   = 1
     L   = Loc( FirstAcoustic ) + 1

     DO Medium = FirstAcoustic, LastAcoustic
        xh2   = x * h( Medium ) ** 2
        h_rho = h( Medium ) * rho( Loc( Medium ) + 1 )

        IF ( Medium >= FirstAcoustic + 1 ) THEN
           L      = L + 1
           d( J ) = ( d( J ) + ( B1( L ) - xh2 ) / h_rho ) / 2.0D0
        ENDIF

        DO ii = 1, N( Medium )
           J      = J + 1
           L      = L + 1
           d( J ) = ( B1( L ) - xh2 ) / h_rho

           IF ( REAL( B1( L ) - xh2 ) + 2.0D0 > 0.0D0 ) THEN   ! Find index of turning point nearest top
              ITP = MIN( J, ITP )
           ENDIF
        END DO

     END DO

     ! Corner elt requires bottom impedance
     CALL BCImpedance( x, 'BOT', HSBot, F, G, IPower )

     IF ( G == 0.0D0 ) THEN
        d( NTotal1 ) = 1.0D0
        e( NTotal1 ) = 0.0D0
     ELSE
        d( NTotal1 ) = d( NTotal1 ) / 2.0D0 - F / G
     ENDIF

     CALL SINVIT( NTotal1, d, e, IERR, RV1, RV2, RV3, RV4, Phi )   ! Inverse iteration to compute eigenvector

     IF ( IERR /= 0 ) THEN
        WRITE( PRTFile, * ) 'mode = ', mode
        CALL ERROUT( PRTFile, 'W', 'KRAKENC-SINVIT', 'Inverse iteration failed to converge' )
        ! EVMat(  1, mode ) = 0.0D0   ! remove that eigenvalue
        ! Extrap( 1, mode ) = 0.0D0
     ELSE
        CALL Normalize( Phi, ITP, NTotal1, x )   ! Normalize the eigenvector

        ! Tabulate the modes at the source/rcvr depths and write to disk
        PhiTab = CMPLX( Phi( IzTab ) ) + WTS * CMPLX( Phi( IzTab + 1 ) - Phi( IzTab ) )
        WRITE( MODFile, REC = IRecProfile + 5 + mode ) PhiTab
     END IF
  END DO

  DEALLOCATE( z, e, d, Phi, RV1, RV2, RV3, RV4 )
  DEALLOCATE( WTS, IzTab, PhiTab )

END SUBROUTINE VECTOR

!**********************************************************************

SUBROUTINE Normalize( Phi, ITP, NTotal1, x )

  ! Normalize the eigenvector:
  ! SqNorm = Integral(Phi ** 2) by the trapezoidal rule:
  ! Integral( F ) = H * ( F(1) + ... + F(N-1) + 0.5 * ( F(0) + F(N) ) )

  ! Compute perturbation due to material absorption
  ! Call ScatterLoss to figure interfacial scatter loss

  USE krakcmod
  IMPLICIT NONE
  INTEGER,           INTENT( IN    ) :: NTotal1, ITP
  INTEGER                            :: IPower, J, J1, L, L1, Medium
  REAL     (KIND=8)                  :: rhoMedium, rho_omega_h2
  COMPLEX  (KIND=8), INTENT( IN    ) :: x
  COMPLEX  (KIND=8), INTENT( INOUT ) :: Phi( NTotal1 )
  COMPLEX  (KIND=8)                  :: x1, F1, G1, x2, F2, G2, DrhoDx, DetaDx, SqNorm, RN, ScaleFactor, Slow

  SqNorm = 0.0D0
  Slow   = 0.0D0

  ! Compute contribution from the top half-space
  IF ( HSTop%BC == 'A' ) THEN
     Slow  = Slow + Phi( 1 ) ** 2 / ( 2 * SQRT( x - omega2 / HSTop%cP ** 2 ) ) / ( HSTop%rho * HSTop%cP ** 2 )
  ENDIF

  ! Compute contribution from the volume
  L = Loc( FirstAcoustic )
  J = 1

  DO Medium = FirstAcoustic, LastAcoustic
     L            = L + 1
     rhoMedium    = rho( L )
     rho_omega_h2 = rhoMedium * omega2 * h( Medium ) ** 2

     ! top interface
     SqNorm = SqNorm + 0.5D0 * h( Medium ) *                    Phi( J ) ** 2 / rhoMedium
     Slow   = Slow   + 0.5D0 * h( Medium ) * ( B1( L ) + 2. ) * Phi( J ) ** 2 / rho_omega_h2

     ! medium
     L1 = L + 1
     L  = L + N( Medium ) - 1
     J1 = J + 1
     J  = J + N( Medium ) - 1

     SqNorm = SqNorm + h( Medium ) * SUM(                         Phi( J1 : J ) ** 2 ) / rhoMedium
     Slow   = Slow   + h( Medium ) * SUM( ( B1( L1 : L ) + 2. ) * Phi( J1 : J ) ** 2 ) / rho_omega_h2

     ! bottom interface
     L = L + 1
     J = J + 1

     SqNorm = SqNorm + 0.5D0 * h( Medium ) *                    Phi( J ) ** 2 / rhoMedium
     Slow   = Slow   + 0.5D0 * h( Medium ) * ( B1( L ) + 2. ) * Phi( J ) ** 2 / rho_omega_h2

  END DO

  ! Compute contribution from the bottom half-space

  IF ( HSBot%BC == 'A' ) THEN
     Slow  = Slow + Phi( J ) ** 2 / ( 2 * SQRT( x - omega2 / HSBot%cP ** 2 ) ) / ( HSBot%rho * HSBot%cP ** 2 )
  ENDIF

  ! Compute derivative of top admitance
  x1 = 0.9999999D0 * x
  x2 = 1.0000001D0 * x

  CALL BCImpedance( x1, 'TOP', HSTop, F1, G1, IPower )
  CALL BCImpedance( x2, 'TOP', HSTop, F2, G2, IPower )

  DrhoDx = 0.0D0
  IF ( G1 /= 0.0D0 ) DrhoDx = ( F2 / G2 - F1 / G1 ) / ( x2 - x1 )

  ! Compute derivative of bottom admitance

  CALL BCImpedance( x1, 'BOT', HSBot, F1, G1, IPower )
  CALL BCImpedance( x2, 'BOT', HSBot, F2, G2, IPower )

  DetaDx = 0.0D0
  IF ( G1 /= 0.0D0 ) DetaDx = ( F2 / G2 - F1 / G1 ) / ( x2 - x1 )

  ! Scale the mode
  RN          = SqNorm - DrhoDx * Phi( 1 ) ** 2 + DetaDx * Phi( NTotal1 ) ** 2
  ScaleFactor = 1.0D0 / SQRT( RN )
  IF ( REAL( Phi( ITP ) ) < 0.0D0 ) ScaleFactor = -ScaleFactor  ! make sign consistent at mode turning point

  Phi        = ScaleFactor * Phi
  Slow       = ScaleFactor ** 2 * Slow * SQRT( omega2 / x )
  VG( mode ) = DBLE( 1 / Slow )

  CALL ScatterLoss( Phi, x )   ! Compute interfacial scatter loss

END SUBROUTINE Normalize

!**********************************************************************

SUBROUTINE ScatterLoss( Phi, x )

  ! Figure scatter loss

  USE krakcmod
  IMPLICIT NONE
  INTEGER           :: J, L, Medium
  REAL     (KIND=8) :: rho1, rho2, rhoInside, h2
  COMPLEX  (KIND=8) :: Phi( * ), x, PERK, eta1Sq, eta2Sq, KupIng, U, PekerisRoot

  PERK = 0.0D0
  J = 1
  L = Loc( FirstAcoustic )

  DO Medium = FirstAcoustic - 1, LastAcoustic   ! Loop over media

     ! Calculate rho1, eta1Sq, Phi, U

     IF ( Medium == FirstAcoustic - 1 ) THEN   ! Top properties
        SELECT CASE ( HSTop%BC )
        CASE ( 'A' )          ! Acousto-elastic
           rho1      = HSTop%rho
           eta1Sq    = x - omega2 / HSTop%cP ** 2
           U         = PekerisRoot( eta1Sq ) * Phi( 1 ) / HSTop%rho
        CASE ( 'V' )          ! Vacuum
           rho1      = 1.0D-9
           eta1Sq    = 1.0D0
           rhoInside = rho( Loc( FirstAcoustic ) + 1 )
           U         = Phi( 2 ) / h( FirstAcoustic ) / rhoInside
        CASE ( 'R' )          ! Rigid
           rho1      = 1.0D+9
           eta1Sq    = 1.0D0
           U         = 0.0D0
        CASE DEFAULT          ! Tabulated
           rho1      = 0.0D0
           eta1Sq    = 0.0D0
           U         = 0.0D0
        END SELECT
     ELSE
        h2 = h( Medium ) ** 2
        J  = J + N( Medium )
        L  = Loc( Medium ) + N( Medium ) + 1

        rho1   = rho( L )
        eta1Sq = ( 2.0D0 + B1( L ) ) / h2 - x
        U      = ( -Phi( J - 1 ) - 0.5D0 * ( B1( L ) - h2 * x ) * Phi( J ) ) / ( h( Medium ) * rho1 )
     ENDIF

     ! Calculate rho2, eta2

     IF ( Medium == LastAcoustic ) THEN   ! Bottom properties
        SELECT CASE ( HSBot%BC )
        CASE ( 'A' )          ! Acousto-elastic
           rho2   = HSBot%rho
           eta2Sq = omega2 / HSBot%cP ** 2 - x
        CASE ( 'V' )          ! Vacuum
           rho2   = 1.0D-9
           eta2Sq = 1.0D0
        CASE ( 'R' )          ! Rigid
           rho2   = 1.0D+9
           eta2Sq = 1.0D0
        CASE DEFAULT          ! Tabulated
           rho2   = 0.0D0
           eta2Sq = 0.0D0
        END SELECT
     ELSE
        rho2   = rho( L + 1 )
        eta2Sq = ( 2.0D0 + B1( L + 1 ) ) / h( Medium + 1 ) ** 2 - x
     ENDIF

     PERK = PERK + KupIng( sigma( Medium + 1 ), eta1Sq, rho1, eta2Sq, rho2, Phi( J ), U )

  END DO

  k( mode ) = PERK

END SUBROUTINE ScatterLoss

!**********************************************************************

SUBROUTINE ORDER( x, N )

  ! Does an insertion sort of the complex vector x in order of decreasing real part

  ! At the Ith step, the first I-1 positions contain a sorted vector.
  ! We shall insert the Ith value into its place in that
  ! vector, shifting up to produce a new vector of length I.

  IMPLICIT NONE
  INTEGER, INTENT( IN ) :: N
  INTEGER               :: ILeft, IMiddle, IRight, I
  COMPLEX (KIND=8)      :: x( N ), xTemp

  IF ( N == 1 ) RETURN

  DO I = 2, N

     xTemp = x( I )

     IF ( REAL( xTemp ) > REAL( x( 1 ) ) ) THEN
        x( 2 : I ) = x( 1 : I - 1 )
        x( 1 )     = xTemp  ! goes in the first position
     ELSE IF ( REAL( xTemp ) > REAL( x( I - 1 ) ) ) THEN

        ! Binary search for its place

        IRight = I - 1
        ILeft  = 1

        DO WHILE ( IRight > ILeft + 1 )
           IMiddle = ( ILeft + IRight ) / 2

           IF ( REAL( xTemp ) > REAL( x( IMiddle ) ) ) THEN
              IRight = IMiddle
           ELSE
              ILeft  = IMiddle
           ENDIF
        END DO

        x( IRight + 1 : I ) = x( IRight : I - 1 )
        x( IRight ) = xTemp

     ENDIF

  END DO

END SUBROUTINE ORDER
