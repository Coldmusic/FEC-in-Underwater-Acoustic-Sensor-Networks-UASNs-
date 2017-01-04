PROGRAM BELLHOP

  ! BELLHOP Beam tracing for ocean acoustics

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

  ! First version (1983) originally developed with Homer Bucker, Naval Ocean Systems Center

  USE bellMod
  USE RefCoMod
  USE bdryMod
  USE angleMod
  USE SdRdRMod
  USE ArrMod
  USE BeamPatternMod

  IMPLICIT NONE
  INTEGER,       PARAMETER :: SHDFile = 25, RAYFile = 21, ArrivalsStorage = 2000000
  REAL (KIND=8), PARAMETER :: DegRad = pi / 180.0D0

  INTEGER              :: IBPvec( 1 ), is, ibeam, iSingle, iBeamWindow, iBeamWindow2, ibp, Ird1, Irec, NBeamsOpt, NImage
  REAL                 :: Tstart, Tstop, DeltaR
  REAL        (KIND=8) :: freq, Amp0, Dalpha, alpha0, DalphaOpt, xs( 2 ), RadMax, epsMultiplier, deltas, Rloop, s, &
                          c, gradc( 2 ), crr, crz, czz, rBox, zBox
  COMPLEX, ALLOCATABLE :: U( :, : )
  COMPLEX     (KIND=8) :: eps, PickEpsilon
  CHARACTER ( LEN=80 ) :: Title, FileRoot
  CHARACTER ( LEN=3  ) :: BotOpt, BeamType
  CHARACTER ( LEN=5  ) :: RunType
  CHARACTER ( LEN=1  ) :: Component

  CALL CPU_TIME( Tstart )

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  ! Read in control data

  CALL READIN(  FileRoot, Title, freq, iSingle, NImage, iBeamWindow, deltas, zBox, rBox, epsMultiplier, rLoop,  &
                TopOpt, BotOpt, HSTop, HSBot, RunType, BeamType, Component )
  CALL READATI( FileRoot, TopOpt( 5 : 5 ), HSTop%Depth, rBox, PRTFile )    ! READ AlTImetry
  CALL READBTY( FileRoot, BotOpt( 2 : 2 ), HSBot%Depth, rBox, PRTFile )    ! READ BaThYmetry
  CALL READRC(  FileRoot, BotOpt( 1 : 1 ), TopOpt(2 : 2),     PRTFile )    ! READ Reflection Coefficients (top and bottom)
  CALL READPAT( FileRoot, RunType( 3 : 3 ),                   PRTFile )    ! Read Source Beam Pattern

  IF ( RunType( 5 : 5 ) == 'I' ) THEN
     Nrd_per_range = 1     ! irregular grid
  ELSE
     Nrd_per_range = Nrd   ! rectilinear grid
  ENDIF

  ! for a TL calculation, allocate space for the pressure matrix
  IF ( SCAN( 'CSI', RunType( 1 : 1 ) ) /= 0 ) THEN
     ALLOCATE ( U( Nrd_per_range, Nr ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) &
          CALL ERROUT( PRTFile, 'F', 'BELLHOP', 'Insufficient memory for TL matrix: reduce Nr * Nrd'  )
  ELSE
     ALLOCATE ( U( 1, 1 ), Stat = IAllocStat )
  ENDIF

  ! for an arrivals run, allocate space for arrivals matrices
  IF ( SCAN( 'Aa', RunType( 1 : 1 ) ) /= 0 ) THEN
     MaxNArr = MAX( ArrivalsStorage / ( Nrd_per_range * Nr ), 10 )   ! allow space for at least 10 arrivals
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) '( Maximum # of arrivals = ', MaxNArr, ')'

     ALLOCATE ( Arr( Nrd_per_range, Nr, MaxNArr ), NArr( Nrd_per_range, Nr ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP', &
          'Insufficient memory to allocate arrivals matrix; reduce parameter ArrivalsStorage' )
  ELSE
     MaxNArr = 1
     ALLOCATE ( Arr( Nrd_per_range, Nr, 1 ), NArr( Nrd_per_range, Nr ), Stat = IAllocStat )
  END IF

  IF ( Nr > 1 ) THEN
     DeltaR = r( Nr ) - r( Nr - 1 )
  ELSE
     DeltaR = 0.0
  ENDIF

  omega  = 2.0 * pi * freq
  alpha  = DegRad * alpha   ! convert to radians
  Dalpha = 0.0
  IF ( NBeams /= 1 ) Dalpha = ( alpha( NBeams ) - alpha( 1 ) ) / ( NBeams - 1 )  ! angular spacing between beams

  WRITE( *, * )
  DO IS = 1, Nsd   ! Loop over source depth
     xs = [ 0.0, sd( IS ) ]   ! source coordinate

     IF ( SCAN( 'CSI', RunType( 1 : 1 ) ) /= 0 ) U    = 0.0 ! For a TL run, zero out pressure matrix
     IF ( SCAN( 'Aa',  RunType( 1 : 1 ) ) /= 0 ) NArr = 0   ! For an arrivals run, zero out arrival matrix

     CALL SSP( xs, c, gradc, crr, crz, czz, TopOpt, 'TAB' )
     RadMax = 5 * c / freq  ! 5 wavelength max radius

     ! Are there enough beams?
     DalphaOpt = SQRT( c / ( 6.0 * freq * r( Nr ) ) )
     NBeamsOpt = 2 + INT( ( alpha( NBeams ) - alpha( 1 ) ) / DalphaOpt )

     IF ( RunType( 1 : 1 ) == 'C' .AND. NBeams < NBeamsOpt ) THEN
        CALL ERROUT( PRTFile, 'W', 'BELLHOP', 'Too few beams' )
        WRITE( PRTFile, * ) 'NBeams should be at least = ', NBeamsOpt
     ENDIF

     ! Trace successive beams

     DO ibeam = 1, NBeams

        IF ( iSingle == 0 .OR. ibeam == iSingle ) THEN    ! Single beam run?

           alpha0 = alpha( ibeam ) * 180.0 / pi   ! take-off angle in degrees

           IBPvec = maxloc( SrcBmPat( :, 1 ), mask = SrcBmPat( :, 1 ) < alpha0 )       ! index of ray angle in beam pattern
           IBP    = IBPvec( 1 )
           IBP    = MAX( IBP, 1 )               ! don't go before beginning of table
           IBP    = MIN( IBP, NSBPPts - 1 )     ! don't go past end of table

           ! linear interpolation to get amplitude
           s    = ( alpha0  - SrcBmPat( IBP, 1 ) ) / ( SrcBmPat( IBP + 1, 1 ) - SrcBmPat( IBP, 1 ) )
           Amp0 = ( 1 - s ) * SrcBmPat( IBP, 2 ) + s * SrcBmPat( IBP + 1, 2 )

           ! show progress ...
           IF ( MOD( ibeam, max( NBeams / 50, 1 ) ) == 1 ) THEN
              WRITE( *, * ) 'Tracing beam ', ibeam, alpha0
           END IF

           CALL TraceRay( deltas, xs, alpha( ibeam ), Amp0, BeamType, zBox, rBox, RunType )   ! Trace a ray

           IF ( RunType( 1 : 1 ) == 'R' ) THEN     ! Write the ray trajectory to RAYFile
              CALL WriteRay( alpha0, Nsteps )
           ELSE                                ! Compute the contribution to the field

              Eps = PickEpsilon( BeamType( 1 : 1 ), omega, c, gradc, alpha( ibeam ), Dalpha, rLoop, epsMultiplier ) ! 'optimal' beam constant

              SELECT CASE ( RunType( 2 : 2 ) )
              CASE ( 'R' )
                 iBeamWindow2 = iBeamWindow **2
                 RadMax       = 50 * c / freq  ! 50 wavelength max radius
                 CALL InfluenceR(   U, DeltaR, Eps, alpha( IBeam ), NImage, IBeamWindow2, RunType, RadMax, BeamType, Component )
              CASE ( 'C' )
                 iBeamWindow2 = iBeamWindow **2
                 RadMax       = 50 * c / freq  ! 50 wavelength max radius
                 CALL InfluenceC(   U, DeltaR, Eps, alpha( IBeam ), NImage, IBeamWindow2, RunType, RadMax, BeamType )
              CASE ( 'S' )
                 CALL InfluenceSGB(         U,            alpha( IBeam ), RunType, Dalpha, deltas )
              CASE ( 'B' )
                 CALL InfluenceGeoGaussian( U,  sd( IS ), alpha( IBeam ), RunType, Dalpha, deltas )
              CASE DEFAULT
                 CALL InfluenceGeoHat(      U,  sd( IS ), alpha( IBeam ), RunType, Dalpha         )
              END SELECT

           END IF
        END IF
     END DO ! Next beam

     ! write results to disk

     SELECT CASE ( RunType( 1 : 1 ) )
     CASE ( 'C', 'S', 'I' )   ! TL calculation
        CALL ScalePressure( Dalpha, ray( 1 )%c, R, U, Nrd_per_range, Nr, RunType, TopOpt, freq )
        IRec  = 7 + Nrd_per_range * ( IS - 1 )
        DO Ird1 = 1, Nrd_per_range
           IRec = IRec + 1
           WRITE( SHDFile, REC = IRec ) U( Ird1, 1 : Nr )
        END DO

     CASE ( 'A' )             ! arrivals calculation, ascii
        CALL WriteArrivalsASCII(  R, Nrd_per_range, Nr, TopOpt, freq, RunType( 4 : 4 ) )
     CASE ( 'a' )             ! arrivals calculation, binary
        CALL WriteArrivalsBinary( R, Nrd_per_range, Nr, TopOpt, freq, RunType( 4 : 4 ) )
     END SELECT

  END DO    ! Next source depth

  ! close all files
  SELECT CASE ( RunType( 1 : 1 ) )
  CASE ( 'C', 'S', 'I' )      ! TL calculation
     CLOSE( SHDFile )
  CASE ( 'A' )                ! arrivals calculation, ascii
     CLOSE( ARRFile )
  CASE ( 'a' )                ! arrivals calculation, binary
     CLOSE( ARRFile )
  CASE ( 'R' )                ! ray trace
     CLOSE( RAYFile )
  END SELECT

  ! Display run time
  CALL CPU_TIME( Tstop )
  WRITE( PRTFile, "( /, ' CPU Time = ', G15.3, 's' )" ) Tstop - Tstart

END PROGRAM BELLHOP

! **********************************************************************!

FUNCTION PickEpsilon( BeamType, omega, c, gradc, alpha, Dalpha, rLoop, EpsMultiplier )

  ! Picks the optimum value for epsilon

  IMPLICIT NONE
  INTEGER, PARAMETER :: PRTFile = 6
  COMPLEX, PARAMETER :: i = ( 0.0, 1.0 )
  LOGICAL, SAVE      :: INIFlag = .TRUE.
  REAL      (KIND=8) :: omega, alpha, Dalpha, epsMultiplier, Rloop, c, gradc( 2 ), cz, halfwidth
  COMPLEX   (KIND=8) :: PickEpsilon, EpsOpt
  CHARACTER (LEN=1 ) :: BeamType
  CHARACTER (LEN=40) :: TAG

  SELECT CASE ( BeamType )
  CASE ( 'F' )
     TAG       = 'Space filling beams'
     halfwidth = 2.0 / ( ( omega / c ) * Dalpha )
     EpsOpt    = i * 0.5 * omega * halfwidth ** 2
  CASE ( 'M' )
     TAG       = 'Minimum width beams'
     halfwidth = SQRT( 2.0 * c * 1000.0 * rLoop / omega )
     EpsOpt    = i * 0.5 * omega * halfwidth ** 2
  CASE ( 'W' )
     TAG       = 'WKB beams'
     halfwidth = HUGE( halfwidth )
     cz        = gradc( 2 )
     IF ( cz == 0.0 ) THEN
        EpsOpt = 1.0E10
     ELSE
        EpsOpt = ( -SIN( alpha ) / COS( alpha ** 2 ) ) * c * c / cz
     ENDIF
  CASE ( 'C' )
     TAG    = 'Cerveny style beam'
  END SELECT

  PickEpsilon = EpsMultiplier * EpsOpt

  ! On first call write info to prt file
  IF ( INIFlag ) THEN
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) TAG
     WRITE( PRTFile, * ) 'halfwidth  = ', halfwidth
     WRITE( PRTFile, * ) 'EpsOpt     = ', EpsOpt
     WRITE( PRTFile, * ) 'EpsMult    = ', EpsMultiplier
     WRITE( PRTFile, * )
     INIFlag = .FALSE.
  END IF

END FUNCTION PickEpsilon

! **********************************************************************!

SUBROUTINE TraceRay( deltas, xs, alpha, Amp0, BeamType, zBox, rBox, RunType )

  ! Traces the beam corresponding to a particular take-off angle

  USE bellMod
  USE bdryMod
  USE RefCoMod

  IMPLICIT NONE
  INTEGER           :: IsegTopT( 1 ), IsegBotT( 1 ), IsegTop, IsegBot, is, iStep
  REAL     (KIND=8) :: Amp0, alpha, c, gradc( 2 ), crr, crz, czz,  &
                       xs( 2 ), dEndTop( 2 ), dEndBot( 2 ), TopnInt( 2 ), BotnInt( 2 ), &
                       ToptInt( 2 ), BottInt( 2 ), DistBegTop, DistEndTop, DistBegBot, DistEndBot, &
                       rTopseg( 2 ), rBotseg( 2 ), deltas, sss, rBox, zBox
  CHARACTER (LEN=3) :: BeamType
  CHARACTER (LEN=5) :: RunType

  ! Initial conditions

  CALL SSP( xs, c, gradc, crr, crz, czz, TopOpt, 'TAB' )
  NumTopBnc( 1 ) = 0
  NumBotBnc( 1 ) = 0
  ray(  1    )%c = c
  ray(  1    )%x = xs
  ray(  1    )%t = [ COS( alpha ), SIN( alpha ) ] / c
  ray( 1     )%p = [ 1.0, 0.0 ]
  ray( 1     )%q = [ 0.0, 1.0 ]
  ray( 1   )%tau = 0.0
  ray( 1   )%Amp = Amp0
  ray( 1 )%Phase = 0.0

  ! second component of qv is not used in geometric beam tracing
  ! set I.C. to 0 in hopes of saving run time
  IF ( RunType( 2 : 2 ) == 'G' ) ray( 1 )%q = [ 0.0, 0.0 ]

  ! identify the top segment above the source

  IsegTopT = MAXLOC( Top( 1 : NatiPts )%x( 1 ), Top( 1 : NatiPts )%x( 1 ) <= xs( 1 ) )

  IF ( IsegTopT( 1 ) > 0 .AND. IsegTopT( 1 ) < NatiPts ) THEN
     IsegTop  = IsegTopT( 1 )   ! IsegTop MUST LIE IN [ 1, NatiPts-1 ]
     rTopSeg = [ Top( IsegTop )%x( 1 ), Top( IsegTop + 1 )%x( 1 ) ]   ! segment limits in range
  ELSE
     CALL ERROUT( PRTFile, 'F', 'TraceRay', 'Top altimetry undefined above the source' )
  ENDIF

  ! identify the bottom segment below the source

  IsegBotT = MAXLOC( Bot( 1 : NbtyPts )%x( 1 ), Bot( 1 : NbtyPts )%x( 1 ) <= xs( 1 ) )

  IF ( IsegBotT( 1 ) > 0 .AND. IsegBotT( 1 ) < NbtyPts ) THEN
     IsegBot  = IsegBotT( 1 )   ! IsegBot MUST LIE IN [ 1, NbtyPts-1 ]
     rBotSeg = [ Bot( IsegBot )%x( 1 ), Bot( IsegBot + 1 )%x( 1 ) ]   ! segment limits in range
  ELSE
     CALL ERROUT( PRTFile, 'F', 'TraceRay', 'Bottom bathymetry undefined below the source' )
  ENDIF

  ! Trace the beam
  ! (note that Reflect alters the step index is)
  is         = 0
  CALL Distances( ray( 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
                              Top( IsegTop )%n, Bot( IsegBot )%n, DistBegTop, DistBegBot )
  ! !!! note above distance is really a negative distance (throughout the code)

  IF ( DistBegTop >= 0 .OR. DistBegBot >= 0 ) THEN
     Nsteps = 1
     RETURN       ! source must be within the medium
  END IF

  Stepping: DO IStep = 1, MaxN - 1
     is = is + 1
     CALL Step( ray( is ), ray( is + 1 ),  &
          Top( IsegTop )%x, Top( IsegTop )%n, &
          Bot( IsegBot )%x, Bot( IsegBot )%n, rTopSeg, rBotSeg, deltas, TopOpt )

     NumTopBnc( is + 1 ) = NumTopBnc( is )
     NumBotBnc( is + 1 ) = NumBotBnc( is )

     ! New altimetry segment?
     IF ( ray( is + 1 )%x( 1 ) < rTopSeg( 1 ) .OR. ray( is + 1 )%x( 1 ) > rTopSeg( 2 ) ) THEN
        IsegTopT = MAXLOC( Top( : )%x( 1 ), Top( : )%x( 1 ) < ray( is + 1 )%x( 1 ) )
        IF ( IsegTopT( 1 ) > 0 .AND. IsegTopT( 1 ) < NatiPts ) IsegTop  = IsegTopT( 1 )  ! IsegTop MUST LIE IN [ 1, NatiPts-1 ]
        rTopSeg = [ Top( IsegTop )%x( 1 ), Top( IsegTop + 1 )%x( 1 ) ]   ! segment limits in range
     END IF

     ! New bathymetry segment?
     IF ( ray( is + 1 )%x( 1 ) < rBotSeg( 1 ) .OR. ray( is + 1 )%x( 1 ) > rBotSeg( 2 ) ) THEN
        IsegBotT = MAXLOC( Bot( : )%x( 1 ), Bot( : )%x( 1 ) < ray( is + 1 )%x( 1 ) )
        IF ( IsegBotT( 1 ) > 0 .AND. IsegBotT( 1 ) < NbtyPts ) IsegBot  = IsegBotT( 1 )  ! IsegBot MUST LIE IN [ 1, NbtyPts-1 ]
        rBotSeg = [ Bot( IsegBot )%x( 1 ), Bot( IsegBot + 1 )%x( 1 ) ]   ! segment limits in range
     END IF

     ! Reflections?
     ! Tests that ray at step i is inside, and ray at step i+1 is outside
     ! to detect only a crossing from inside to outside

     CALL Distances( ray( IS + 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
                                      Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     IF      ( DistBegTop < 0.0d0 .AND. DistEndTop >= 0.0d0 ) THEN  ! test top reflection
        IF ( atiType == 'C' ) THEN
           sss = DOT_PRODUCT( dEndTop, Top( IsegTop )%t ) / Top( IsegTop )%Len   ! proportional distance along segment
           TopnInt = ( 1 - sss ) * Top( IsegTop )%Noden + sss * Top( 1 + IsegTop )%Noden
           ToptInt = ( 1 - sss ) * Top( IsegTop )%Nodet + sss * Top( 1 + IsegTop )%Nodet
        ELSE
           TopnInt = Top( IsegTop )%n   ! normal is constant in a segment
           ToptInt = Top( IsegTop )%t
        END IF

        CALL Reflect( is, BeamType, HSTop, 'TOP', ToptInt, TopnInt, &
             Top( IsegTop )%kappa, RTop, NTopPTS )

        NumTopBnc( is + 1 ) = NumTopBnc( is ) + 1

        CALL Distances( ray( IS + 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
                                         Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     ELSE IF ( DistBegBot < 0.0d0 .AND. DistEndBot >= 0.0d0 ) THEN  ! test bottom reflection
        IF ( btyType == 'C' ) THEN
           sss = DOT_PRODUCT( dEndBot, Bot( IsegBot )%t ) / Bot( IsegBot )%Len   ! proportional distance along segment
           BotnInt = ( 1 - sss ) * Bot( IsegBot )%Noden + sss * Bot( 1 + IsegBot )%Noden
           BottInt = ( 1 - sss ) * Bot( IsegBot )%Nodet + sss * Bot( 1 + IsegBot )%Nodet
        ELSE
           BotnInt = Bot( IsegBot )%n   ! normal is constant in a segment
           BottInt = Bot( IsegBot )%t
        END IF

        CALL Reflect( is, BeamType, HSBot, 'BOT', BottInt, BotnInt, &
             Bot( IsegBot )%kappa, RBot, NBotPTS )

        NumBotBnc( is + 1 ) = NumBotBnc( is ) + 1

        CALL Distances( ray( IS + 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
                                         Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     END IF

     ! Has the ray left the box, lost its energy, escaped the boundaries, or exceeded storage limit?
     IF ( ABS( ray( is + 1 )%x( 1 ) ) > rBox .OR. ABS( ray( is + 1 )%x( 2 ) ) > zBox .OR. ray( is + 1 )%Amp < 0.005 .OR. &
         ( DistBegTop > 0.0 .AND. DistEndTop > 0.0 ) .OR. &
         ( DistBegBot > 0.0 .AND. DistEndBot > 0.0 ) ) THEN
        Nsteps = is + 1
        EXIT Stepping
     ELSE IF ( is >= MaxN - 3 ) THEN
        CALL ERROUT( PRTFile, 'W', 'TraceRay', 'Insufficient storage for ray trajectory' )
        Nsteps = is
        EXIT Stepping
     END IF

     DistBegTop = DistEndTop
     DistBegBot = DistEndBot

  END DO Stepping   ! Next step

END SUBROUTINE TraceRay

! **********************************************************************!

SUBROUTINE Distances( rayx, Topx, Botx, dTop, dBot, Topn, Botn, DistTop, DistBot )

! Input:
!    rayx   ray coordinate
!    Topx   top    boundary coordinate
!    Botx   bottom boundary coordinate
!    Topn   top    boundary normal
!    Botn   bottom boundary normal
!
! Output:
!    dTop    vector from top    boundary to the ray coordinate
!    dBot    vector from bottom boundary to the ray coordinate
!    DistTop normal 'distance' to top    boundary
!    DistBot normal 'distance' to bottom boundary

  IMPLICIT NONE
  REAL (KIND=8), INTENT( IN  ) :: rayx( 2 ), Topx( 2 ), Botx( 2 ), Topn( 2 ), Botn( 2 )
  REAL (KIND=8), INTENT( OUT ) :: dTop( 2 ), dBot( 2 ), DistTop, DistBot 

  dTop    = rayx - Topx  ! vector pointing from top    to ray
  dBot    = rayx - Botx  ! vector pointing from bottom to ray
  DistTop = DOT_PRODUCT( Topn, dTop )
  DistBot = DOT_PRODUCT( Botn, dBot )

END SUBROUTINE Distances

! **********************************************************************!

SUBROUTINE Step( ray0, ray2, xTop, nTop, xBot, nBot, rTopSeg, rBotSeg, deltas, TopOpt )

  ! Does a single step along the ray
  ! x denotes the ray coordinate, (r,z)
  ! t denotes the scaled tangent to the ray (previously (rho, zeta))
  ! c * t would be the unit tangent

  USE sspMod
  IMPLICIT NONE
  TYPE rayPt
      REAL (KIND=8 ) :: x( 2 ), t( 2), p( 2 ), q( 2 ), tau, c, Amp, Phase
  END TYPE
  TYPE( rayPt )      :: ray0, ray1, ray2
  INTEGER            :: LAYER0
  REAL     (KIND=8 ) :: gradc0( 2 ), gradc1( 2 ), gradc2( 2 ), &
                        c0, crr0, crz0, czz0, csq0, cnn0_csq0, &
                        c1, crr1, crz1, czz1, csq1, cnn1_csq1, &
                        c2, crr2, crz2, czz2, &
                        xTop( 2 ), nTop( 2 ), xBot( 2 ), nBot( 2 ), rTopSeg( 2 ), rBotSeg( 2 ), &
                        deltas, h, halfh, hw0, hw1, ray2n( 2 ), RM, RN, gradcjump( 2 ), cnjump, csjump, w0, w1 
  CHARACTER (LEN=6) :: TopOpt

  ! The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
  ! to the Heun (second order Runge-Kutta method).
  ! However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

  ! *** Phase 1 (an Euler step)

  CALL SSP( ray0%x, c0, gradc0, crr0, crz0, czz0, TopOpt, 'TAB' )

  csq0      = c0 * c0
  cnn0_csq0 = crr0 * ray0%t( 2 )**2 - 2.0 * crz0 * ray0%t( 1 ) * ray0%t( 2 ) + czz0 * ray0%t( 1 )**2
  Layer0    = Layer     ! make note of current layer

  h = deltas            ! initially set the step h, to the basic one, deltas
  CALL ReduceStep( ray0%x, ray0%t, Layer0, c0, xTop, nTop, xBot, nBot, rTopSeg, rBotSeg, deltas, h ) ! reduce h to land on boundary
  halfh = 0.5 * h   ! first step of the modified polygon method is a half step

  !$OMP PARALLEL SHARED( halfh, c0 )
  !$OMP WORKSHARE
  ray1%x = ray0%x + halfh * c0 * ray0%t
  ray1%t = ray0%t - halfh * gradc0 / csq0
  ray1%p = ray0%p - halfh * cnn0_csq0 * ray0%q
  ray1%q = ray0%q + halfh * c0        * ray0%p
  !$OMP END WORKSHARE
  !$OMP END PARALLEL

  ! *** Phase 2

  CALL SSP( ray1%x, c1, gradc1, crr1, crz1, czz1, TopOpt, 'TAB' )
  csq1      = c1 * c1
  cnn1_csq1 = crr1 * ray1%t( 2 )**2 - 2.0 * crz1 * ray1%t( 1 ) * ray1%t( 2 ) + czz1 * ray1%t( 1 )**2

  ! The Munk test case with a horizontally launched ray caused problems.
  ! The ray vertexes on an interface and can ping-pong around that interface.
  ! Have to be careful in that case about big changes to the stepsize (that invalidate the leap-frog scheme) in phase II.
  ! A modified Heun or Box method could also work.

  !h = deltas            ! initially set the step h, to the basic one, deltas
  CALL ReduceStep( ray0%x, ray1%t, Layer0, c1, xTop, nTop, xBot, nBot, rTopSeg, rBotSeg, deltas, h ) ! reduce h to land on boundary

  ! use blend of f' based on proportion of a full step used.
  w1  = h / ( 2.0d0 * halfh )
  w0  = 1.0d0 - w1
  hw0 = h * w0
  hw1 = h * w1

  !$OMP PARALLEL SHARED( c0, c1, hw0, hw1 )
  !$OMP WORKSHARE
  ray2%x   = ray0%x   + hw0 * c0 * ray0%t        + hw1 * c1 * ray1%t
  ray2%t   = ray0%t   - hw0 * gradc0 / csq0      - hw1 * gradc1 / csq1
  ray2%p   = ray0%p   - hw0 * cnn0_csq0 * ray0%q - hw1 * cnn1_csq1 * ray1%q
  ray2%q   = ray0%q   + hw0 * c0        * ray0%p + hw1 * c1        * ray1%p
  ray2%tau = ray0%tau + hw0 / c0                 + hw1 / c1
  !$OMP END WORKSHARE
  !$OMP END PARALLEL
  ray2%Amp   = ray0%Amp
  ray2%Phase = ray0%Phase

  ! If we crossed an interface, apply jump condition

  CALL SSP( ray2%x, c2, gradc2, crr2, crz2, czz2, TopOpt, 'TAB' )
  ray2%c = c2

  IF ( Layer /= Layer0 ) THEN
     gradcjump =  gradc2 - gradc0  ! this is precise only for c-linear layers
     ray2n     = [ -ray2%t( 2 ), ray2%t( 1 ) ]

     cnjump    = DOT_PRODUCT( gradcjump, ray2n  )
     csjump    = DOT_PRODUCT( gradcjump, ray2%t )

     RM        = ray2%t( 1 ) / ray2%t( 2 )
     RN        = RM * ( 2 * cnjump - RM * csjump ) / c2
     RN        = -RN
     ray2%p    = ray2%p + ray2%q * RN

  END IF

END SUBROUTINE Step

! **********************************************************************!

SUBROUTINE ReduceStep( x0, rayt, Layer0, c, xTop, nTop, xBot, nBot, rTopSeg, rBotSeg, deltas, h )

  USE sspMod
  IMPLICIT NONE
  INTEGER          Layer0
  REAL (KIND=8) :: c, x( 2 ), x0( 2 ), rayt( 2 ), &
                   xTop( 2 ), nTop( 2 ), xBot( 2 ), nBot( 2 ), rTopSeg( 2 ), rBotSeg( 2 ), &
                   d( 2 ), e( 2 ), h, h1, h2, h3, h4, h5, deltas

  ! Detect interface or boundary crossing and reduce step, if necessary, to land on that crossing.
  ! Need to keep in mind possibility that user put source right on an interface
  ! and that multiple events can occur (crossing interface, top, and bottom in a single step).
  ! reminder: rayt is a scaled tangent; c * rayt is the unit tangent

  x = x0 + h * c * rayt ! make a trial step
  !$OMP PARALLEL

  !$OMP SECTIONS
  !$OMP SECTION
  ! interface crossing in depth
  h1 = huge( h1 )
  IF ( rayt( 2 ) /= 0.0d0 ) THEN
     IF      ( zSSPV( Layer0     ) > x(  2 ) ) THEN
        h1 = ( zSSPV( Layer0     ) - x0( 2 ) ) / ( rayt( 2 ) * c )
     ELSE IF ( zSSPV( Layer0 + 1 ) < x(  2 ) ) THEN
        h1 = ( zSSPV( Layer0 + 1 ) - x0( 2 ) ) / ( rayt( 2 ) * c )
     END IF
  END IF

  !$OMP SECTION
  ! top crossing
  h2 = huge( h2 )
  d  = x - xTop              ! vector from top to ray
  IF ( DOT_PRODUCT( nTop, d ) > 0.0d0 ) THEN
     e  = x0 - xTop   ! vector from top    node to ray origin
     h2 = -DOT_PRODUCT( e, nTop ) / ( c * DOT_PRODUCT( rayt, nTop ) )
  END IF

  !$OMP SECTION
  ! bottom crossing
  h3 = huge( h3 )
  d  = x - xBot              ! vector from bottom to ray
  IF ( DOT_PRODUCT( nBot, d ) > 0.0d0 ) THEN
     e  = x0 - xBot   ! vector from bottom node to ray origin
     h3 = -DOT_PRODUCT( e, nBot ) / ( c * DOT_PRODUCT( rayt, nBot ) )
  END IF

  !$OMP SECTION
  ! top segment crossing in range
  h4 = huge( h4 )
  IF ( rayt( 1 ) /= 0.0d0 ) THEN
     IF       ( x(  1 ) < rTopSeg( 1 ) ) THEN
        h4 = -( x0( 1 ) - rTopSeg( 1 ) ) / ( rayt( 1 ) * c )
     ELSE IF  ( x(  1 ) > rTopSeg( 2 ) ) THEN
        h4 = -( x0( 1 ) - rTopSeg( 2 ) ) / ( rayt( 1 ) * c )
     END IF
  END IF

  !$OMP SECTION
  ! bottom segment crossing in range
  h5 = huge( h5 )
  IF ( rayt( 1 ) /= 0.0d0 ) THEN
     IF       ( x(  1 ) < rBotSeg( 1 ) ) THEN
        h5 = -( x0( 1 ) - rBotSeg( 1 ) ) / ( rayt( 1 ) * c )
     ELSE IF  ( x(  1 ) > rBotSeg( 2 ) ) THEN
        h5 = -( x0( 1 ) - rBotSeg( 2 ) ) / ( rayt( 1 ) * c )
     END IF
  END IF
  !$OMP END SECTIONS

  !$OMP END PARALLEL

  h = MIN( h, h1, h2, h3, h4, h5 )  ! take limit set by shortest distance to a crossing
  h = MAX( h, 1.0d-8 * deltas )     ! make sure we make some motion

END SUBROUTINE ReduceStep

! **********************************************************************!

SUBROUTINE Reflect( is, BeamType, HS, BotTop, tbdry, nbdry, kappa, RefC, Npts )

  USE bellMod
  USE sspMod
  USE RefCoMod
  IMPLICIT NONE

  INTEGER           :: is, is1, Npts
  REAL     (KIND=8) :: ck, co, si, cco, ssi, pdelta, rddelta, sddelta, &
                       c, gradc( 2 ), crr, crz, czz, &
                       tbdry( 2 ), nbdry( 2 ), kappa, cN, cS, RM, RN, Tg, Th, &
                       rayn( 2 )
  COMPLEX  (KIND=8) :: gamma1, gamma2, gamma1Sq, gamma2Sq, GK, Refl, ch, a, b, d, sb, delta, ddelta
  CHARACTER (LEN=3) :: BeamType
  CHARACTER (LEN=3) :: BotTop
  TYPE( HSInfo )    :: HS
  TYPE(ReflectionCoef), INTENT( IN ) :: RefC( NPts )
  TYPE(ReflectionCoef) :: RInt

  is  = is + 1
  is1 = is + 1

  NumTopBnc( is1 ) = NumTopBnc( is )
  NumBotBnc( is1 ) = NumBotBnc( is )

  ! here's the geometric part, changing the ray direction
  ray( is1 )%x = ray( is )%x

  Tg = DOT_PRODUCT( ray( is )%t, TBdry )  ! component of ray tangent, along boundary
  Th = DOT_PRODUCT( ray( is )%t, NBdry )  ! component of ray tangent, normal to boundary

  ray( is1 )%t =  ray( is )%t - 2.0 * Th * NBdry
  rayn         = [ -ray( is )%t( 2 ), ray( is )%t( 1 ) ]

  ! Calculate the change in curvature
  ! Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).

  CALL SSP( ray( is1 )%x, c, gradc, crr, crz, czz, TopOpt, 'TAB' )

  ! rayn and ray%t should be multiplied by c to produce unit normals and tangents
  ! Therefore cn, cs are off by that same factor
  cn       = DOT_PRODUCT( gradc, rayn )
  cs       = DOT_PRODUCT( gradc, ray( is )%t )

  RN = 2 * kappa / c**2 / Th    ! boundary curvature correction

  IF ( BotTop == 'TOP' ) THEN
     cn = -cn    ! flip sign for top reflection
     RN = -RN
  END IF

  RM = Tg / Th
  RN = RN + RM * ( 4 * cn - 2 * RM * cs ) / c   ! dividing by c instead of c^2 compensates for the missing factor in cn, cs

  SELECT CASE ( BeamType(2 : 2) )
  CASE ( 'D' )
     RN = 2.0 * RN
  CASE ( 'Z' )
     RN = 0.0
  END SELECT

  ray( is1 )%c = c
  ray( is1 )%p = ray( is )%p + ray( is )%q * RN
  ray( is1 )%q = ray( is )%q

  ! account for phase change

  SELECT CASE ( HS%BC )
  CASE ( 'R' )                 ! rigid
     ray( is1 )%tau   = ray( is )%tau
     ray( is1 )%Amp   = ray( is )%Amp
     ray( is1 )%Phase = ray( is )%Phase
  CASE ( 'V' )                 ! vacuum
     ray( is1 )%tau   = ray( is )%tau
     ray( is1 )%Amp   = ray( is )%Amp
     ray( is1 )%Phase = ray( is )%Phase + pi
  CASE ( 'F' )                 ! file
     RInt%theta = RadDeg * ABS( ATAN2( Th, Tg ) )   ! angle of incidence (relative to normal to bathymetry)
     IF ( RInt%theta > 90 ) RInt%theta = 180. - RInt%theta  ! reflection coefficient is symmetric about 90 degrees
     CALL RefCO( RInt, RefC, Npts, PRTFile )
     ray( is1 )%tau   = ray( is )%tau
     ray( is1 )%Amp   = ray( is )%Amp * RInt%R
     ray( is1 )%Phase = ray( is )%Phase + RInt%phi
  CASE ( 'A' )                 ! half-space
     GK       = omega * Tg   ! wavenumber in direction parallel to bathymetry
     gamma1Sq = ( omega / c     ) ** 2 - GK ** 2 - i * tiny( omega )   ! tiny prevents g95 giving -zero, and wrong branch cut
     gamma2Sq = ( omega / HS%cP ) ** 2 - GK ** 2 - i * tiny( omega )
     gamma1   = SQRT( -gamma1Sq )
     gamma2   = SQRT( -gamma2Sq )

     Refl = ( HS%rho * gamma1 - gamma2 ) / ( HS%rho * gamma1 + gamma2 )

     IF ( ABS( Refl ) < 1.0E-5 ) THEN   ! kill a ray that has lost its energy in reflection
        ray( is1 )%tau   = ray( is )%tau
        ray( is1 )%Amp   = 0.0
        ray( is1 )%Phase = ray( is )%Phase
     ELSE
        ray( is1 )%tau   = ray( is )%tau
        ray( is1 )%Amp   = ABS( Refl ) * ray(  is )%Amp
        ray( is1 )%Phase = ray( is )%Phase + ATAN2( AIMAG( Refl ), REAL( Refl ) )

        ! compute beam-displacement Tindle, Eq. (14)
        ! needs a correction to beam-width as well ...
        !  IF ( REAL( gamma2Sq ) < 0.0 ) THEN
        !     rhoW   = 1.0   ! density of water
        !     rhoWSq  = rhoW  * rhoW
        !     rhoHSSq = rhoHS * rhoHS
        !     DELTA = 2 * GK * rhoW * rhoHS * ( gamma1Sq - gamma2Sq ) /
        ! &( gamma1 * i * gamma2 *
        ! &( -rhoWSq * gamma2Sq + rhoHSSq * gamma1Sq ) )
        !     RV( is + 1 ) = RV( is + 1 ) + DELTA
        !  END IF

        if ( BeamType( 3 : 3 ) == 'S' ) then   ! beam displacement & width change (Seongil's version)

           ch = ray( is )%c / conjg( HS%cP )
           co = ray( is )%t( 1 ) * ray( is )%c
           si = ray( is )%t( 2 ) * ray( is )%c
           ck = omega / ray( is )%c

           a   = 2 * HS%rho * ( 1 - ch * ch )
           b   = co * co - ch * ch
           d   = HS%rho * HS%rho * si * si + b
           sb  = sqrt( b )
           cco = co * co
           ssi = si * si

           delta   = a * co / si / ( ck * sb * d )    
           pdelta  = real( delta ) / ( ray( is )%c / co)

           ddelta  = -a / ( ck*sb*d ) - a*cco / ssi / (ck*sb*d) + a*cco / (ck*b*sb*d) &
                     -a*co / si / (ck*sb*d*d) * (2* HS%rho * HS%rho *si*co-2*co*si)
           rddelta = -real( ddelta )
           sddelta = rddelta / abs( rddelta )        

           ray( is1 )%x( 1 ) = ray( is1 )%x( 1 ) + real( delta )   ! displacement
           ray( is1 )%tau    = ray( is1 )%tau + pdelta             ! phase change
           ray( is1 )%q      = ray( is1 )%q + sddelta * rddelta * si * c * ray( is )%p   ! beam-width change
        endif

     ENDIF
  END SELECT

END SUBROUTINE Reflect

!*******************************************************************!

SUBROUTINE InfluenceR( U, deltar, eps, alpha, NImage, iBeamWindow2, RunType, RadMax, BeamType, Component )

  ! Computes the beam influence, i.e. the contribution of a single beam to the complex pressure
  ! This routine is for the Cerveny-style beams
  ! The 'R' in InfluenceR is for beams in Ray-centered coordinates

  USE bellMod
  USE SdRdRMod

  IMPLICIT NONE
  INTEGER              is, id, ir, ir1, ir2, KMAHV( MaxN ), KMAH, Nimage, image, iBeamWindow2
  REAL              :: deltar
  REAL     (KIND=8) :: nA, nB, nSq, rA, rB, c, W, DS, zr, znv, rnv, tau, ratio1, alpha, n, RadMax, HERMITE
  COMPLEX  (KIND=8) :: pVB( MaxN ), qVB( MaxN ), q, epsV( MaxN ), contri, gammaV( MaxN ), gamma, P_n, P_s, eps
  COMPLEX           :: U( Nrd_per_range, Nr )
  CHARACTER (LEN=5) :: RunType
  CHARACTER (LEN=3) :: BeamType
  CHARACTER (LEN=1) :: Component

  ! need to add logic related to Nrd_per_range

  ! Note that during reflection imag(q) is constant and
  ! adjacent normals cannot bracket a segment of the TL
  ! line, so no special treatment is necessary

  DS = 2 * ( SIN( omega * ray( 1 )%x( 2 ) * ray( 1 )%t( 2 ) ) ) ** 2   ! Lloyd mirror pattern
  
  IF ( BeamType( 1 : 1 ) == 'C' ) THEN
     epsV( 1 : NSteps ) = i * ABS( ray( 1 : NSteps )%q( 1 ) / ray( 1 : NSteps )%q( 2 ) )
  ELSE
     epsV( 1 : Nsteps ) = eps
  ENDIF

  pVB(    1 : NSteps ) = ray( 1 : NSteps )%p( 1 ) + epsV( 1 : NSteps ) * ray( 1 : NSteps )%p( 2 )
  qVB(    1 : NSteps ) = ray( 1 : NSteps )%q( 1 ) + epsV( 1 : NSteps ) * ray( 1 : NSteps )%q( 2 )
  gammaV( 1 : NSteps ) = pVB( 1 : NSteps ) / qVB( 1 : NSteps )

  IF ( RunType( 4 : 4 ) == 'R' ) THEN
     Ratio1 = SQRT( ABS( COS( alpha ) ) )  ! point source
  ELSE
     Ratio1 = 1                            ! line  source
  END IF

  ! compute KMAH index
  ! Following is incorrect for 'Cerveny'-style beamwidth (narrow as possible)
  KMAHV(  1 ) = 1
  DO is = 2, Nsteps
     KMAHV(  is ) = KMAHV( is - 1 )
     CALL BranchCut( qVB( is - 1 ), qVB( is ), BeamType, KMAHV( is ) )
  END DO

  DO id = 1, Nrd_per_range ! Loop over receivers
     zR = RD( id )

     DO Image = 1, NImage ! Loop over images
        ir1 = HUGE( ir1 )

        DO is = 2, Nsteps  ! Loop over steps
           IF ( ABS( ray( is )%x( 1 ) - ray( is - 1 )%x( 1 ) ) < TINY( ray( is )%x( 1 ) ) ) CYCLE   ! don't process duplicate points

           ! Compute ray-centered coordinates, (zNV, RNV)
           zNV = -ray( is )%t( 1 ) * ray( is )%c
           RNV =  ray( is )%t( 2 ) * ray( is )%c
           IF ( ABS( zNV ) < tiny( zNV ) ) THEN   ! Check for normal parallel to TL-line
              CYCLE   ! skip to next step on ray
           ENDIF

           SELECT CASE ( Image )     ! Images of beams
           CASE ( 1 )                ! True beam
              nB  = ( zR -                  ray( is )%x( 2 )   ) / zNV
           CASE ( 2 )                ! Surface-reflected beam
              RNV = -RNV
              nB  = ( zR - ( 2.0 * HSTop%Depth - ray( is )%x( 2 ) ) ) / zNV
           CASE ( 3 )                ! Bottom-reflected beam
              RNV = -RNV
              nB  = ( zR - ( 2.0 * HSBot%Depth - ray( is )%x( 2 ) ) ) / zNV
           END SELECT

           rB  = ray( is )%x( 1 ) + nB * RNV
           ir2 = MAX( MIN( INT( ( rB - r( 1 ) ) / deltar ) + 1, Nr ), 1 ) ! index of receiver

           IF ( ir1 >= ir2 ) THEN
              rA  = rB
              nA  = nB
              ir1 = ir2
              CYCLE
           END IF

           DO ir = ir1 + 1, ir2    ! Compute influence for each rcvr
              W     = ( r( ir ) - rA ) / ( rB - rA )
              q     =    qVB( is-1 ) + W * (    qVB( is ) -    qVB( is-1 ) )
              gamma = gammaV( is-1 ) + W * ( gammaV( is ) - gammaV( is-1 ) )
              n     = nA + W * ( nB - nA )
              nSq   = n * n
              IF ( AIMAG( gamma ) > 0 ) THEN
                 WRITE( PRTFile, * ) 'Unbounded beam'
                 CYCLE   ! next receiver depth
              ENDIF

              IF ( -0.5 * omega * AIMAG( gamma ) * nSq < iBeamWindow2 ) THEN   ! Within beam window?
                 c      = ray( is-1 )%c
                 tau    = ray( is-1 )%tau + W * ( ray( is )%tau - ray( is-1 )%tau )
                 contri = ratio1 * ray( is )%Amp * SQRT( c * ABS( epsV( is ) ) / q ) * &
                      EXP( -i * ( omega * ( tau + 0.5 * gamma * nSq ) - ray( is )%phase ) )

                 SELECT CASE ( Component )
                 CASE ( 'P' )   ! pressure
                 CASE ( 'V' )   ! vertical component
                    P_n    = -i * omega * gamma * n * contri
                    P_s    = -i * omega / c         * contri
                    contri = c * DOT_PRODUCT( [ P_n, P_s ], ray( is )%t ) 
                 CASE ( 'H' )   ! horizontal component
                    P_n    = -i * omega * gamma * n * contri
                    P_s    = -i * omega / c         * contri
                    contri = c * ( -P_n * ray( is )%t( 2 ) + P_s * ray( is )%t( 1 ) ) 
                 END SELECT

                 KMAH = KMAHV( is - 1 )
                 CALL BranchCut( qVB( is - 1 ), q, BeamType, KMAH ) ! Get correct branch of SQRT

                 IF ( KMAH  < 0  ) contri = -contri
                 IF ( Image == 2 ) contri = -contri

                 SELECT CASE ( RunType( 1 : 1 ) )
                 CASE ( 'I' )    ! Incoherent TL
                    contri =      ABS(contri)
                 CASE ( 'S' )    ! Semi-coherent TL
                    contri = DS * ABS(contri)
                 END SELECT
                 U( id, ir ) = U( id, ir ) + CMPLX( HERMITE( n, RadMax, 2 * RadMax ) * contri )

              ENDIF
           END DO   ! next ir
           rA  = rB
           nA  = nB
           ir1 = ir2
        END DO   ! Next step, is
     END DO   ! Next image
  END DO   ! Next receiver depth

END SUBROUTINE InfluenceR
! **********************************************************************!
SUBROUTINE InfluenceC( U, deltar, eps, alpha, NImage, iBeamWindow2, RunType, RadMax, BeamType )

  ! Computes the beam influence, i.e. 
  ! the contribution of a single beam to the complex pressure
  ! This version uses a beam representation in Cartesian coordinates

  USE bellMod
  USE SdRdRMod

  IMPLICIT NONE
  INTEGER           :: KMAHV( MaxN ), KMAH, is, id, ir, ir1, ir2, NImage, IBeamWindow2
  REAL              :: deltar
  REAL     (KIND=8) :: x( 2 ), rayt( 2 ), rayn( 2 ), DS, Tr, Tz, rA, rB, zr, W, tau, ratio1, &
                       alpha, c, cs, cn, csq, gradc( 2 ), crr, crz, czz, RadMax, Hermite, deltaz
  COMPLEX  (KIND=8) :: pVB( MaxN ), qVB( MaxN ), q, epsV( MaxN ), contri, gammaV( MaxN ), gamma, const, eps
  COMPLEX           :: U( Nrd_per_range, Nr )
  CHARACTER (LEN=5) :: RunType
  CHARACTER (LEN=3) :: BeamType

  ! missing logic for Nrd_per_range

  ! need to add logic related to Nrd_per_range

  ! Note that during reflection imag(q) is constant and
  ! adjacent normals cannot bracket a segment of the TL
  ! line, so no special treatment is necessary

  DS = 2 * ( SIN( omega * ray( 1 )%x( 2 ) * ray( 1 )%t( 2 ) ) ) ** 2   ! Lloyd mirror pattern
  
  IF ( BeamType( 1 : 1 ) == 'C' ) THEN
     epsV( 1 : NSteps ) = i * ABS( ray( 1 : NSteps )%q( 1 ) / ray( 1 : NSteps )%q( 2 ) )
  ELSE
     epsV( 1 : Nsteps ) = eps
  ENDIF

  pVB(    1 : NSteps ) = ray( 1 : NSteps )%p( 1 ) + epsV( 1 : NSteps ) * ray( 1 : NSteps )%p( 2 )
  qVB(    1 : NSteps ) = ray( 1 : NSteps )%q( 1 ) + epsV( 1 : NSteps ) * ray( 1 : NSteps )%q( 2 )

  ! Note that during reflection imag(q) is constant and
  ! adjacent normals cannot bracket a segment of the TL
  ! line, so no special treatment is necessary

  IF ( RunType( 4 : 4 ) == 'R' ) THEN
     Ratio1 = SQRT( ABS( COS( alpha ) ) )  ! point source
  ELSE
     Ratio1 = 1                            ! line source
  END IF

  ! Form gamma and KMAH index
  ! Treatment of KMAH index is incorrect for 'Cerveny' style beam width BeamType

  DO is = 1, Nsteps

     rayt = ray( is )%c * ray( is )%t   ! unit tangent
     rayn = [ rayt( 2 ), -rayt( 1 ) ]   ! unit normal

     CALL SSP( ray( is )%x, c, gradc, crr, crz, czz, TopOpt, 'TAB' )

     csq = c * c
     cS  = DOT_PRODUCT( gradc, rayt )
     cN  = DOT_PRODUCT( gradc, rayn )

     Tr  = rayt(  1 )
     Tz  = rayt(  2 )

     gammaV( is ) = 0.0
     IF ( qVB( is ) /= 0.0 ) gammaV( is ) = 0.5 * ( pVB( is ) / qVB( is ) * Tr**2 + &
          2.0 * cN / csq * Tz * Tr - cS / csq * Tz**2 )

     IF ( is == 1 ) THEN
        KMAHV( 1 ) = 1
     ELSE
        KMAHV( is ) = KMAHV( is - 1 )
        CALL BranchCut( qVB( is - 1 ), qVB( is ), BeamType, KMAHV( is ) )
     ENDIF

     !RLTEMP = SQRT( -2.0 / ( omega * AIMAG( pVB( is ) / qVB( is ) ) ) )
     !RKTEMP = -cV( is ) * REAL( pVB( is ) / qVB( is ) )
     !write( *, * ) is, rltemp, rktemp

  END DO

  ! Loop over steps

  DO is = 3, Nsteps
     IF ( ray( is     )%x( 1 ) > r( Nr ) ) RETURN
     rA = ray( is - 1 )%x( 1 )
     rB = ray( is     )%x( 1 )
     IF ( ABS( rB - rA ) < TINY( rB ) ) CYCLE   ! don't process duplicate points

     ! Compute upper index on rcvr line
     ! Note: assumes r is a vector of equally spaced points
     ir1 = MAX( MIN( INT( ( rA - r( 1 ) ) / deltar ) + 1, Nr ), 1 ) ! should be ", 0 )" ?
     ir2 = MAX( MIN( INT( ( rB - r( 1 ) ) / deltar ) + 1, Nr ), 1 )

     IF ( ir1 >= ir2 ) CYCLE

     DO ir = ir1 + 1, ir2

        W     = ( r( ir ) - rA ) / ( rB - rA )

        x     = ray(    is-1 )%x    + W * ( ray(       is )%x    -  ray(    is-1 )%x )
        rayt  = ray(    is-1 )%t    + W * ( ray(       is )%t    -  ray(    is-1 )%t )
        c     = ray(    is-1 )%c    + W * ( ray(       is )%c    -  ray(    is-1 )%c )
        q     = qVB(    is-1 )      + W * ( qVB(       is )      -  qVB(    is-1 ) )
        tau   = ray(    is-1 )%tau  + W * ( ray(       is )%tau  -  ray(    is-1 )%tau )
        gamma = gammaV( is-1 )      + W * ( gammaV(    is )      -  gammaV( is-1 ) )

        IF ( AIMAG( gamma ) > 0 ) THEN
           WRITE( PRTFile, * ) 'Unbounded beam'
           WRITE( PRTFile, * ) gammaV( is-1 ), gammaV( is ), gamma
           CYCLE   ! next receiver
        ENDIF

        const = ratio1 * SQRT( c * ABS( epsV( is - 1 ) ) / q )

        ! Get correct branch of SQRT
        KMAH = KMAHV( is - 1 )
        CALL BranchCut( qVB( is - 1 ), q, BeamType, KMAH )
        IF ( KMAH < 0 ) const = -const

        DO id = 1, Nrd_per_range  ! Loop over receivers
           zR = RD( id )

                                   ! True beam
              deltaz = zR - x( 2 )
              IF ( omega * AIMAG( gamma ) * deltaz**2 < iBeamWindow2 ) &
                   contri =           ray( is )%Amp * HERMITE( deltaz, RadMax, 2.0 * RadMax ) * &
                   EXP( -i * ( omega * ( tau + rayt( 2 ) * deltaz + gamma * deltaz**2 ) - ray( is )%Phase ) )

           IF ( NImage >= 2 ) THEN ! Surface reflected beam
              deltaz = -zR + 2.0 * HSTop%Depth - x( 2 )
              IF ( omega * AIMAG( gamma ) * deltaz**2 < iBeamWindow2 ) &
                   contri =  contri - ray( is )%Amp * HERMITE( deltaz, RadMax, 2.0 * RadMax ) * &
                   EXP( -i * ( omega * ( tau + rayt( 2 ) * deltaz + gamma * deltaz**2 ) - ray( is )%Phase ) )
           ENDIF

           IF ( NImage >= 3 ) THEN ! Bottom reflected beam
              deltaz = -zR + 2.0 * HSBot%Depth - x( 2 )
              IF ( omega * AIMAG( gamma ) * deltaz**2 < iBeamWindow2 ) &
                   contri =  contri + ray( is )%Amp * HERMITE( deltaz, RadMax, 2.0 * RadMax ) * &
                   EXP( -i * ( omega * ( tau + rayt( 2 ) * deltaz + gamma * deltaz**2 ) - ray( is )%Phase ) )
           ENDIF

           ! contribution to field
           SELECT CASE( RunType( 1 : 1 ) )
           CASE ( 'C' )   ! coherent
              contri = const * contri
           CASE ( 'I' )   ! incoherent
              contri =ABS( const * contri )
           CASE ( 'S' )   ! semi-coherent
              contri = DS * ABS( const * contri )
           END SELECT
           U( id, ir ) = U( id, ir ) + CMPLX( contri )
        END DO   ! Next receiver depth
     END DO   ! Next receiver range
  END DO   ! Next step along the ray

END SUBROUTINE InfluenceC

! **********************************************************************!

SUBROUTINE InfluenceGeoHat( U, zs, alpha, RunType, Dalpha )

  ! Computes the beam influence, i.e.
  ! the contribution of a single beam to the complex pressure
  ! This version uses geometric, hat-shaped beams

  USE bellMod
  USE SdRdRMod
  USE ArrMod

  IMPLICIT NONE
  INTEGER           :: is, id, ir, irT( 1 ), irTT
  REAL              :: zs
  REAL     (KIND=8) :: xray( 2 ), rayt( 2 ), rayn( 2 ), xrcvr( 2 ), s, n, phaseInt, &
                       rA, rB, rLen, RadMax, zMin, zMax, q0, qold, A, const, ratio1, W, DS, &
                       Amp, phase, delay, dqds, dtauds, q, alpha, Dalpha, SrcAngle, RcvrAngle
  COMPLEX           :: U( Nrd_per_range, Nr )
  CHARACTER (LEN=5) :: RunType
  CHARACTER (LEN=1) :: RunTypeE

  ! some ugly code to have RunType 'a' treated like 'A'
  RunTypeE = RunType( 1 : 1 )
  IF ( RunTypeE == 'a' ) RunTypeE = 'A'

  DS       = SQRT( 2.0 ) * SIN( omega * zs * ray( 1 )%t( 2 ) )   ! Lloyd mirror pattern
  q0       = ray( 1 )%c / Dalpha   ! Reference for J = q0 / q
  SrcAngle = RadDeg * alpha        ! take-off angle in degrees
  phase    = 0.0
  qOld     = ray( 1 )%q( 1 )
  rA       = ray( 1 )%x( 1 )       ! range at start of ray

  ! what if never satistified?
  irT = MINLOC( r( 1 : nr ), MASK = r( 1 : nr ) .GT. rA )   ! find index of first receiver to the right of rA
  ir  = irT( 1 )
  IF ( ray( 1 )%t( 1 ) .LT. 0.0d0 ) ir = ir - 1  ! if ray is traveling to the left, then we want the first receiver to the left of rA
  !ir = min( 1, nr )

  IF ( RunType( 4 : 4 ) == 'R' ) THEN
     Ratio1 = SQRT( ABS( COS( alpha ) ) )  ! point source
  ELSE
     Ratio1 = 1                            ! line  source
  END IF

  DO is = 2, Nsteps  ! Loop over steps
     rB   = ray( is     )%x( 1 )
     xray = ray( is - 1 )%x
     IF ( ABS( rB - rA ) < TINY( rA ) ) CYCLE   ! jump to next step if duplicate point

     dqds   = ray( is )%q( 1 ) -  ray( is-1 )%q( 1 )
     dtauds = ray( is )%tau     - ray( is-1 )%tau

     ! initialize the index of the receiver range
     IF ( is == 2 ) THEN
        IF ( rB > rA ) THEN   ! ray is moving right
           ir = 1   ! index all the way to the left
        ELSE                  ! ray is moving left
           ir = nr  ! index all the way to the right
        END IF
     END IF

     ! compute normalized tangent (compute it because we need to measure the step length)
     rayt = ray( is )%x - ray( is - 1 )%x
     rlen = sqrt( rayt( 1 ) **2 + rayt( 2 ) **2 )
     !IF ( rlen < 1d-7 * deltas ) CYCLE   ! if duplicate point in ray, skip to next step along the ray
     IF ( rlen < 100 * TINY( ray( is )%x( 1 ) ) ) CYCLE   ! if duplicate point in ray, skip to next step along the ray

     rayt = rayt / rlen                    ! unit tangent to ray
     rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal  to ray

     ! phase shifts at caustics
     q  = ray( is-1 )%q( 1 )
     IF ( q <= 0.0d0 .AND. qOld > 0.0d0 .OR. q >= 0.0d0 .AND. qOld < 0.0d0 ) phase = phase + pi / 2.  ! phase shifts at caustics
     qold = q

     RadMax = MAX( ABS( ray( is -1 )%q( 1 ) ), ABS( ray( is )%q( 1 ) ) ) / q0 / abs( rayt( 1 ) ) ! beam radius projected onto vertical line
     zmin   = min( ray( is - 1 )%x( 2 ), ray( is )%x( 2 ) ) - RadMax  ! min depth of ray segment
     zmax   = max( ray( is - 1 )%x( 2 ), ray( is )%x( 2 ) ) + RadMax  ! max depth of ray segment

     ! is this a steep ray? Then don't try to get depth limits: it's too complicated
     IF ( ABS( rayt( 1 ) ) < 0.5 ) THEN
        zmin = -HUGE( zmin )
        zmax = +HUGE( zmax )
     END IF

     ! computed beam influence for this segment of the ray
     DO
        ! is r( ir ) contained in [ rA, rB }? Then compute beam influence
        !IF ( ABS( r( ir ) - rA ) + ABS( rB - r( ir ) ) <= ABS( rB - rA ) ) THEN
        IF ( r( ir ) >= MIN( rA, rB ) .AND. r( ir ) < MAX( rA, rB ) ) THEN
           DO id = 1, Nrd_per_range  ! Loop over receiver depths
              IF ( RunType( 5 : 5 ) == 'I' ) THEN
                 xrcvr = [ r( ir ), rd( ir ) ]   ! irregular   grid
              ELSE
                 xrcvr = [ r( ir ), rd( id ) ]   ! rectilinear grid
              ENDIF
              IF ( xrcvr( 2 ) < zmin .OR. xrcvr( 2 ) > zmax ) CYCLE

              s      =      dot_product( xrcvr - xray, rayt ) / rlen ! proportional distance along ray
              n      = abs( dot_product( xrcvr - xray, rayn ) )      ! normal distance to ray
              q      = ray( is-1 )%q( 1 ) + s * dqds                 ! interpolated amplitude
              RadMax = ABS( q / q0 )                                 ! beam radius
              IF ( n < RadMax ) THEN
                 A        = 1 / RadMax
                 delay    = ray( is-1 )%tau + s * dtauds             ! interpolated delay
                 const    = Ratio1 * SQRT( ray( is )%c / ABS( q ) ) * A * ray( is )%Amp
                 phaseInt = phase
                 IF ( q <= 0.0d0 .AND. qOld > 0.0d0 .OR. &
                      q >= 0.0d0 .AND. qOld < 0.0d0 ) phaseInt = phase + pi / 2. ! phase shifts at caustics

                 IF ( RunTypeE == 'S' ) const = DS * const   ! semi-coherent TL
                 Amp   = const * ( RadMax - n )

                 SELECT CASE( RunTypeE )
                 CASE ( 'E' )      ! eigenrays
                    CALL WriteRay( SrcAngle, is )
                 CASE ( 'A' )      ! arrivals
                    RcvrAngle  = RadDeg * ATAN2( rayt( 2 ), rayt( 1 ) )
                    CALL AddArr( omega, id, ir, Amp, ray( is - 1 )%Phase + phaseInt, delay, &
                         SrcAngle, RcvrAngle, NumTopBnc( is ), NumBotBnc( is ) )
                 CASE ( 'C'  )     ! coherent TL
                    U( id, ir ) = U( id, ir ) + CMPLX( Amp * EXP( -i * ( omega * delay - ray( is - 1 )%Phase - phaseInt ) ) )
                 CASE DEFAULT      ! incoherent/semi-coherent TL
                    W           = ( RadMax - n ) / RadMax   ! hat function: 1 on center, 0 on edge
                    U( id, ir ) = U( id, ir ) + SNGL( ( Amp / W ) ** 2 * W )
                 END SELECT
              ENDIF
           END DO   ! Next receiver depth
        END IF

        ! bump receiver index, ir, towards rB
        IF ( r( ir ) < rB ) THEN
           IF ( ir >= nr ) EXIT   ! jump out of the search and go to next step on ray
           irTT = ir + 1          ! bump right
           IF ( r( irTT ) >= rB ) EXIT
        ELSE
           IF ( ir <= 1  ) EXIT   ! jump out of the search and go to next step on ray
           irTT = ir - 1          ! bump left
           IF ( r( irTT ) <= rB ) EXIT
        END IF
        ir = irTT
     END DO   ! Next receiver range

     rA = rB
  END DO   ! Next step along the ray

END SUBROUTINE InfluenceGeoHat


! **********************************************************************!

SUBROUTINE InfluenceGeoGaussian( U, zs, alpha, RunType, Dalpha, deltas )

  ! Computes the beam influence, i.e.
  ! the contribution of a single beam to the complex pressure
  ! This version uses geometric, Gaussian beams

  USE bellMod
  USE SdRdRMod
  USE ArrMod

  IMPLICIT NONE
  INTEGER, PARAMETER :: BeamWindow = 4   ! beam window: kills beams outside e**(-0.5 * ibwin**2 )
  INTEGER            :: is, id, ir, irT( 1 ), irTT
  REAL               :: zs
  REAL      (KIND=8) :: xray( 2 ), rayt( 2 ), rayn( 2 ), xrcvr( 2 ), s, n, phaseInt, &
                        rA, rB, rLen, RadMax, zMin, zMax, q, q0, qold, sigma, lambda, A, const, ratio1, W, DS, &
                        alpha, dalpha, deltas, Amp, phase, delay, SrcAngle, RcvrAngle
  COMPLEX            :: U( Nrd_per_range, Nr )
  CHARACTER  (LEN=5) :: RunType
  CHARACTER  (LEN=1) :: RunTypeE

  ! some ugly code to have RunType 'a' treated like 'A'
  RunTypeE = RunType( 1 : 1 )
  IF ( RunTypeE == 'a' ) RunTypeE = 'A'

  DS       = SQRT( 2.0 ) * SIN( omega * zs * ray( 1 )%t( 2 ) )   ! Lloyd mirror pattern
  q0       = ray( 1 )%c / Dalpha     ! Reference for J = q0 / q
  SrcAngle = RadDeg * alpha          ! take-off angle in degrees
  phase    = 0
  qOld     = ray( 1 )%q( 1 )
  rA       = ray( 1 )%x( 1 )         ! range at start of ray

  ! what if never satistified?
  irT       = MINLOC( r( 1 : nr ), MASK = r( 1 : nr ) .GT. rA )      ! find index of first receiver to the right of rA
  ir        = irT( 1 )

  IF ( ray( 1 )%t( 1 ) .LT. 0.0d0 ) ir = ir - 1  ! if ray is left-traveling, get the first receiver to the left of rA

  ! factor 1.2535 represents a sum of Gaussians in free space
  IF ( RunType( 4 : 4) == 'R' ) THEN
     Ratio1 = SQRT( ABS( COS( alpha ) ) ) / 1.2535   ! point source
  ELSE
     Ratio1 = 1 / 1.2535                             ! line  source
  END IF

  DO is = 2, Nsteps   ! Loop over steps

     rB   = ray( is     )%x( 1 )
     xray = ray( is - 1 )%x

     ! compute normalized tangent (compute it because we need to measure the step length)
     rayt = ray( is )%x - ray( is - 1 )%x
     rlen = sqrt( rayt( 1 )**2 + rayt( 2 ) **2 )
     rayt = rayt / rlen
     rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal to ray

     ! phase shifts at caustics
     q  = ray( is-1 )%q( 1 )
     IF ( q <= 0.0 .AND. qOld > 0.0 .OR. q >= 0.0 .AND. qOld < 0.0 ) phase = phase + pi / 2.  ! phase shifts at caustics
     qold = q

     lambda = ray( is - 1 )%c / ( omega / ( 2 * pi ) )

     sigma  = MAX( ABS( ray( is -1 )%q( 1 ) ), ABS( ray( is )%q( 1 ) ) ) / q0 / abs( rayt( 1 ) ) ! beam radius projected onto vertical line
     sigma  = MAX( sigma, MIN( 0.2 * is * deltas / lambda, pi * lambda ) )
     RadMax = BeamWindow * sigma

     zmin   = min( ray( is - 1 )%x( 2 ), ray( is )%x( 2 ) ) - RadMax   ! min depth of ray segment
     zmax   = max( ray( is - 1 )%x( 2 ), ray( is )%x( 2 ) ) + RadMax   ! max depth of ray segment

     ! is this a steep ray?
     ! If so, don't try to get depth limits: it's too complicated
     IF ( ABS( rayt( 1 ) ) < 0.5 ) THEN
        zmin = -HUGE( zmin )
        zmax = +HUGE( zmax )
     END IF

    ! computed beam influence for this segment of the ray
     DO
        ! is r( ir ) contained in [ rA, rB )?
        ! If yes, then compute beam influence

        IF ( r( ir ) >= MIN( rA, rB ) .AND. r( ir ) < MAX( rA, rB ) ) THEN

        DO id = 1, Nrd_per_range  ! Loop over receiver depths
           IF ( RunType( 5 : 5 ) == 'I' ) THEN
              xrcvr = [ r( ir ), rd( ir ) ]   ! irregular   grid
           ELSE
              xrcvr = [ r( ir ), rd( id ) ]   ! rectilinear grid
           ENDIF
           IF ( xrcvr( 2 ) < zmin .OR. xrcvr( 2 ) > zmax ) CYCLE

           s      = dot_product( xrcvr - xray, rayt ) / rlen  ! proportional distance along ray
           n      = dot_product( xrcvr - xray, rayn )         ! normal distance to ray
           n      = ABS( n )
           q      = ray( is-1 )%q( 1 ) + s * ( ray( is )%q( 1 ) -  ray( is-1 )%q( 1 ) )   ! interpolated amplitude
           sigma  = ABS( q / q0 )                             ! beam radius

           ! calculate the beamwidth (must be at least pi * lambda, except in the nearfield)
           lambda = ray( is - 1 )%c / ( omega / ( 2 * pi ) )
           sigma  = MAX( sigma, MIN( 0.2 * is * deltas / lambda, pi * lambda ) )

           IF ( n < BeamWindow * sigma ) THEN   ! Within beam window?
              A        = ABS( q0 / q )
              !c        = ray( is-1 )%c   + s * ( ray( is )%c -   ray( is-1 )%c )       ! sound speed
              delay    = ray( is-1 )%tau + s * ( ray( is )%tau - ray( is-1 )%tau )     ! interpolated delay
              const    = Ratio1 * SQRT( ray( is )%c / abs( q ) ) * ray( is )%Amp
              phaseInt = phase
              IF ( q <= 0.0d0 .AND. qOld > 0.0d0 .OR. q >= 0.0d0 .AND. qOld < 0.0d0 ) phaseInt = phase + pi / 2.  ! phase shifts at caustics
              IF ( RunTypeE == 'S' ) const = DS * const   ! semi-coherent TL
              Amp      = const * EXP( -0.5 * ( n / sigma )**2 ) / ( 2. * sigma * A )

              SELECT CASE( RunTypeE )
              CASE ( 'E' )                ! eigenrays
                 CALL WriteRay( SrcAngle, is )
              CASE ( 'A' )                ! arrivals
                 RcvrAngle  = RadDeg * ATAN2( rayt(2), rayt( 1 ) )
                 CALL AddArr( omega, id, ir, Amp, ray( is )%Phase + phaseInt, delay, &
                      SrcAngle, RcvrAngle, NumTopBnc( is ), NumBotBnc( is ) )
              CASE ( 'C' )                ! coherent TL
                 U( id, ir ) = U( id, ir ) + CMPLX( Amp * EXP( -i * ( omega * delay - ray( is )%Phase - phaseInt ) ) )
              CASE DEFAULT                ! incoherent/semicoherent TL
                 W           =  EXP( -0.5 * ( n / sigma )**2 ) / ( 2. * sigma * A )   ! Gaussian decay
                 U( id, ir ) = U( id, ir ) + SNGL( ( Amp / W ) ** 2 * W )
              END SELECT
           END IF
        END DO   ! Next receiver depth
        END IF

        ! receiver not bracketted; bump receiver index, ir, towards rB
        IF ( rB > r( ir ) ) THEN
           IF ( ir >= nr ) EXIT   ! jump out of the search and go to next step on ray
           irTT = ir + 1          ! bump right
           IF ( r( irTT ) >= rB ) EXIT
        ELSE
           IF ( ir <= 1  ) EXIT   ! jump out of the search and go to next step on ray
           irTT = ir - 1          ! bump left
           IF ( r( irTT ) <= rB ) EXIT
        END IF
        ir = irTT

     END DO   ! Next receiver range

     rA = rB
  END DO   ! Next step along the ray

END SUBROUTINE InfluenceGeoGaussian


! **********************************************************************!

SUBROUTINE InfluenceSGB( U, alpha, RunType, Dalpha, deltas )

  ! Computes the beam influence, i.e. 
  ! the contribution of a single beam to the complex pressure
  ! This version uses a beam representation in Cartesian coordinates
  ! This is an implementation of Bucker's Simple Gaussian Beams

  USE bellMod
  USE SdRdRMod
  USE ArrMod

  IMPLICIT NONE
  INTEGER           :: id, ir, is
  REAL (KIND=8)     :: x( 2 ), rayt( 2 ), A, beta, cn, CPA, deltaz, DS, phase, q, qOld, rA,rB, Ratio1, &
                       sint, SX, thet, W, tau, SrcAngle, alpha, Dalpha, deltas, delay
  COMPLEX           :: U( Nrd_per_range, Nr )
  COMPLEX (KIND=8)  :: contri
  CHARACTER (LEN=5) :: RunType

  Ratio1 = SQRT(  COS( alpha ) )
  phase  = 0
  qOld   = 1.0
  BETA   = 0.98  ! Beam Factor
  A      = -4.0 * LOG( BETA ) / Dalpha**2
  CN     = Dalpha * SQRT( A / pi )
  rA     = ray( 1 )%x( 1 )
  ir     = 1

  DO is = 2, Nsteps    ! Loop over steps

     rB = ray( is )%x( 1 )
     ! phase shifts at caustics
     q  = ray( is - 1 )%q( 1 )
     IF ( q < 0.0d0 .AND. qOld >= 0.0d0 .OR. q > 0.0d0 .AND. qOld <= 0.0d0 ) phase = phase + pi / 2.  ! phase shifts at caustics
     qold = q

     DO WHILE ( ABS( rB - rA ) > TINY( rA ) .AND. rB > r( ir ) )   ! Loop over bracketted receiver ranges

        W     = ( r( ir ) - rA ) / ( rB - rA )
        x     = ray( is-1 )%x      + W * ( ray( is )%x      - ray( is-1 )%x )
        rayt  = ray( is-1 )%t      + W * ( ray( is )%t      - ray( is-1 )%t )
        q     = ray( is-1 )%q( 1 ) + W * ( ray( is )%q( 1 ) - ray( is-1 )%q( 1 ) )
        tau   = ray( is-1 )%tau    + W * ( ray( is )%tau    - ray( is-1 )%tau )

        ! following is incorrect because ray doesn't always use a step of deltas
        SINT  =  ( is - 1 ) * deltas + W * deltas

        IF ( q < 0.0d0 .AND. qOld >= 0.0d0 .OR. q > 0.0d0 .AND. qOld <= 0.0d0 ) phase = phase + pi / 2. ! phase shifts at caustics

        DO id = 1, Nrd_per_range   ! Loop over receiver depths
           deltaz =  RD( id ) - x( 2 )   ! ray to rcvr distance
           ! Adeltaz    = ABS( deltaz )
           ! IF ( Adeltaz < RadMax ) THEN
           SELECT CASE( RunType( 1 : 1 ) )
           CASE ( 'E' )         ! eigenrays
              SrcAngle = RadDeg * alpha   ! take-off angle in degrees
              CALL WriteRay( SrcAngle, is )

           CASE DEFAULT         ! coherent TL
              CPA    = ABS( deltaz * ( rB - rA ) ) / SQRT( ( rB - rA )**2 + ( ray(is)%x( 2 ) - ray(is-1 )%x( 2 ) )**2  )
              DS     = SQRT( deltaz**2 - CPA**2 )
              SX     = SINT + DS
              thet   = ATAN( CPA / SX )
              delay  = tau + rayt( 2 ) * deltaz
              contri = Ratio1 * CN * ray( is )%Amp * EXP(-A * thet ** 2 - &
                       i * ( omega * delay - ray( is )%Phase - phase ) ) / SQRT( SX )
              U( id, ir ) = U( id, ir ) + CMPLX( contri )

           END SELECT
           ! ENDIF
        END DO   ! Next receiver depth

        qOld = q
        ir   = ir + 1
        IF ( ir > Nr ) RETURN
     END DO   ! Next receiver range

     rA = rB
  END DO   ! Next step along the ray

END SUBROUTINE InfluenceSGB

! **********************************************************************!

SUBROUTINE BranchCut( q1C, q2C, BeamType, KMAH )

  ! Checks for a branch cut crossing and updates KMAH accordingly

  IMPLICIT NONE
  INTEGER   KMAH
  REAL     (KIND=8) :: q1, q2
  COMPLEX  (KIND=8) :: q1C, q2C
  CHARACTER (LEN=3) :: BeamType

  SELECT CASE ( BeamType( 1 : 1 ) )
  CASE ( 'W' )   ! WKBeams
     q1 = REAL( q1C )
     q2 = REAL( q2C )
     IF    ( ( q1 < 0.0 .AND. q2 >= 0.0 ) .OR. ( q1 > 0.0 .AND. q2 <= 0.0 ) ) KMAH = -KMAH
  CASE DEFAULT
     IF ( REAL( q2C ) < 0.0 ) THEN
        q1 = AIMAG( q1C )
        q2 = AIMAG( q2C )
        IF ( ( q1 < 0.0 .AND. q2 >= 0.0 ) .OR. ( q1 > 0.0 .AND. q2 <= 0.0 ) ) KMAH = -KMAH
     ENDIF
  END SELECT

END SUBROUTINE BranchCut

! **********************************************************************!

FUNCTION Hermite( x, x1, x2 )

  ! Calculates a smoothing function based on the h0 hermite cubic
  ! x is the point where the function is to be evaluated
  ! returns:
  ! [  0, x1  ] = 1
  ! [ x1, x2  ] = cubic taper from 1 to 0
  ! [ x2, inf ] = 0

  IMPLICIT NONE
  REAL (KIND=8 ) :: Hermite, x, x1, x2, Ax, u
  Ax  = ABS( x  )

  IF ( Ax <= x1 ) THEN
     HERMITE = 1.0d0
  ELSE IF ( Ax >= x2 ) THEN
     HERMITE = 0.0d0
  ELSE
     u       = ( Ax - x1 ) / ( x2 - x1 )
     HERMITE = ( 1.0d0 + 2.0d0 * u ) * ( 1.0d0 - u ) ** 2
  ENDIF

  !hermit = hermit / ( 0.5 * ( x1 + x2 ) )

END FUNCTION Hermite

! **********************************************************************!

SUBROUTINE ScalePressure( Dalpha, c, r, U, Nrd, Nr, RunType, TopOpt, freq )

  ! Scale the pressure field

  IMPLICIT NONE
  INTEGER           :: Nrd, Nr, ir
  REAL, PARAMETER   :: pi = 3.14159265
  REAL              :: r( Nr )
  REAL     (KIND=8) :: freq, f2, alpha, dalpha, c, const, factor
  COMPLEX           :: U( Nrd, Nr )
  CHARACTER (LEN=6) :: TopOpt
  CHARACTER (LEN=5) :: RunType

  ! Compute scale factor for field
  SELECT CASE ( RunType(2 : 2) )
  CASE ( 'C' )
     const = -Dalpha * SQRT( freq ) / c
  CASE ( 'R' )
     const = -Dalpha * SQRT( freq ) / c
  CASE DEFAULT
     const = -1.0
  END SELECT

  IF ( TopOpt( 4 : 4 ) == 'T' ) THEN ! Thorpe attenuation?
     f2    = ( freq / 1000.0 ) ** 2
     ! Updated formula from JKPS Eq. 1.34
     alpha = 3.3d-3 + 0.11 * f2 / ( 1.0 + f2 ) + 44.0 * f2 / ( 4100.0 + f2 ) + 3d-4* f2   ! dB/km
     alpha = alpha / 8685.8896 ! Nepers / m
  ELSE
     alpha = 0.0
  ENDIF

  IF ( RunType( 1 : 1 ) /= 'C' ) U = SQRT( REAL( U ) ) ! For incoherent run, convert intensity to pressure

  ! add in attenuation
  DO ir = 1, Nr
     IF ( RunType( 4 : 4 ) == 'X' ) THEN   ! line source
        factor = -4.0 * SQRT( pi ) * const
     ELSE                              ! point source
        IF ( r ( ir ) == 0 ) THEN
           factor = 0.0D0              ! avoid /0 at origin, return pressure = 0
        ELSE
           factor = const * EXP( -alpha * r( ir ) ) / SQRT( abs( r( ir ) ) )
        END IF
     END IF
     U( :, ir ) = SNGL( factor ) * U( :, ir )
  END DO

END SUBROUTINE ScalePressure

! **********************************************************************!

SUBROUTINE WriteRay( alpha0, Nsteps1 )

  ! Compress the ray data keeping every iSkip point, points near surface or bottom, and last point.
  ! Write to RAYFile.
  ! During an eigenray calculation, subsets of the full ray may be passed
  ! These have lengths Nsteps1 vs. Nsteps for the entire ray

  USE bellMod
  IMPLICIT NONE
  INTEGER, PARAMETER :: RAYFile = 21
  INTEGER            :: is, N2, iSkip, Nsteps1
  REAL      (KIND=8) :: alpha0

  ! compression

  N2    = 1
  iSkip = MAX( Nsteps1 / 5000, 1 )   ! max #pts written is about 5000

  DO is = 2, Nsteps1
     ! ensure that we always write ray pts. near bdry reflections
     IF ( MIN( HSBot%Depth - ray( is )%x( 2 ),  ray( is )%x( 2 ) - HSTop%Depth ) < 0.2 .OR. &
          MOD( is, iSkip ) == 0 .OR. is == Nsteps ) THEN
        N2 = N2 + 1
        ray( N2 )%x = ray( is )%x
     END IF
  END DO

  ! write to ray file

  WRITE( RAYFile, * ) alpha0
  WRITE( RAYFile, * ) N2, NumTopBnc( Nsteps1 ), NumBotBnc( Nsteps1 )
  DO is = 1, N2
     WRITE( RAYFile, * ) SNGL( ray( is )%x )
  END DO

END SUBROUTINE WriteRay

