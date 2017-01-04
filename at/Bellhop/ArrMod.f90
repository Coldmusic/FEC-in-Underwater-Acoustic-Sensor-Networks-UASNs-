MODULE ArrMod

  ! Variables for arrival information

  INTEGER, PARAMETER    :: ARRFile = 36
  INTEGER               :: MaxNArr
  INTEGER,  ALLOCATABLE :: NArr( :, : )

  TYPE Arrival
     INTEGER :: NTopBnc, NBotBnc
     REAL    :: delay, SrcAngle, RcvrAngle, A, Phase
  END TYPE

  TYPE(Arrival), ALLOCATABLE :: Arr( :, :, : )

CONTAINS

  SUBROUTINE AddArr( omega, id, ir, Amp, Phase, delay, SrcAngle, RcvrAngle, NumTopBnc, NumBotBnc )

    ! ADDs the amplitude and delay for an ARRival into a matrix of same.
    ! Extra logic included to keep only the strongest arrivals.

    IMPLICIT NONE
    REAL,      PARAMETER :: PhaseTol = 0.5  ! arrivals with essentially the same phase are grouped into one
    LOGICAL                 NewRay
    INTEGER                 IArr( 1 ), id, ir, Nt
    INTEGER              :: NumTopBnc, NumBotBnc
    REAL                 :: AmpTot
    REAL    ( KIND = 8 ) :: omega, Amp, Phase, delay, SrcAngle, RcvrAngle

    Nt     = NArr( id, ir )    ! # of arrivals
    NewRay = .TRUE.

    ! Is this the second bracketting ray of a pair?
    ! If so, we want to combine the arrivals to conserve space.
    ! (test this by seeing if the arrival time is close to the previous one)
    ! (also need that the phase is about the same to make sure surface and direct paths are not joined)

    IF ( Nt >= 1 ) THEN
       IF( omega * ABS( delay - Arr( id, ir, Nt )%delay ) < PhaseTol .AND. &
           ABS( Arr( id, ir, Nt )%phase - Phase )       < PhaseTol ) NewRay = .FALSE.
    END IF

    IF ( NewRay ) THEN
       IF ( Nt >= MaxNArr ) THEN       ! space available to add an arrival?
          IARR = MINLOC( Arr( id, ir, : )%A )                       ! no: replace weakest arrival
          IF ( Amp > Arr( id, ir, IArr( 1 ) )%A ) THEN
             Arr( id, ir, IArr( 1 ) )%A         = SNGL( Amp )       ! amplitude
             Arr( id, ir, IArr( 1 ) )%Phase     = SNGL( Phase )     ! phase
             Arr( id, ir, IArr( 1 ) )%delay     = SNGL( delay )     ! delay time
             Arr( id, ir, IArr( 1 ) )%SrcAngle  = SNGL( SrcAngle )  ! angle
             Arr( id, ir, IArr( 1 ) )%RcvrAngle = SNGL( RcvrAngle ) ! angle
             Arr( id, ir, IArr( 1 ) )%NTopBnc   = NumTopBnc         ! Number of top     bounces
             Arr( id, ir, IArr( 1 ) )%NBotBnc   = NumBotBnc         !   "       bottom
          ENDIF
       ELSE
          NArr( id, ir         )           = Nt + 1              ! # of arrivals
          Arr(  id, ir, Nt + 1 )%A         = SNGL( Amp )         ! amplitude
          Arr(  id, ir, Nt + 1 )%Phase     = SNGL( Phase )       ! phase
          Arr(  id, ir, Nt + 1 )%delay     = SNGL( delay )       ! delay time
          Arr(  id, ir, Nt + 1 )%SrcAngle  = SNGL( SrcAngle )    ! angle
          Arr(  id, ir, Nt + 1 )%RcvrAngle = SNGL( RcvrAngle )   ! angle
          Arr(  id, ir, Nt + 1 )%NTopBnc   = NumTopBnc           ! Number of top     bounces
          Arr(  id, ir, Nt + 1 )%NBotBnc   = NumBotBnc           !   "       bottom
       ENDIF
    ELSE      ! not a new ray
       !PhaseArr(   id, ir, Nt ) = PhaseArr( id, ir, Nt )
       AmpTot = Arr( id, ir, Nt )%A + SNGL( Amp )
       Arr( id, ir, Nt )%delay     = ( Arr( id, ir, Nt )%A * Arr( id, ir, Nt )%delay + SNGL( Amp * delay ) ) / AmpTot ! weighted sum
       Arr( id, ir, Nt )%A         = AmpTot
       Arr( id, ir, Nt )%SrcAngle  = SNGL( SrcAngle )
       Arr( id, ir, Nt )%RcvrAngle = SNGL( RcvrAngle )
    ENDIF

    RETURN
  END SUBROUTINE AddArr

  ! **********************************************************************!

  SUBROUTINE WriteArrivalsASCII( R, Nrd, Nr, TopOpt, freq, SourceType )

    ! Writes the arrival data (Amplitude, delay for each eigenray)
    ! ASCII output file

    IMPLICIT NONE
    INTEGER           :: Nrd, Nr, ir, id, IArr
    REAL,   PARAMETER :: PI = 3.14159265, RadDeg = 180 / PI
    REAL                 r( Nr ), AArrT
    REAL     (KIND=8) :: freq, f2, alpha, factor
    CHARACTER (LEN=6) :: TopOpt
    CHARACTER (LEN=1) :: SourceType

    ! *** Thorpe attenuation? ***

    IF ( TopOpt( 4 : 4 ) == 'T' ) THEN
       f2 = ( freq / 1000.0d0 ) ** 2
       ! Updated formula from JKPS Eq. 1.34
       alpha = 3.3d-3 + 0.11d0 * f2 / ( 1.0d0 + f2 ) + 44.0d0 * f2 / ( 4100.0d0 + f2 ) + 3d-4* f2   ! dB/km
       alpha = alpha / 8685.8896 ! Nepers / m
    ELSE
       alpha = 0.0
    ENDIF


    WRITE( ARRFile, * ) MAXVAL( NArr( 1 : Nrd, 1 : Nr ) )

    DO id = 1, Nrd
       DO ir = 1, Nr
          IF ( SourceType == 'X' ) THEN   ! line source
             factor =  4.0 * SQRT( PI ) * EXP( -alpha * r( ir ) )
          ELSE                            ! point source
             IF ( r ( ir ) == 0 ) THEN
                factor = 1e5       ! avoid /0 at origin
             ELSE
                factor = EXP( -alpha * r( ir ) ) / SQRT( r( ir ) )  ! vol. atten. and cyl. spreading
             END IF
          END IF

          WRITE( ARRFile, * ) NArr( id, ir )
          DO IARR = 1, NArr( id, ir )
             AArrT = SNGL( factor * Arr( id, ir, IARR )%A )
             ! you can compress the output file a lot by putting in an explicit format statement here ...
             ! However, you'll need to make sure you keep adequate precision
             WRITE( ARRFile, * ) AArrT, RadDeg * Arr( id, ir, IARR )%Phase, Arr( id, ir, IARR )%delay, &
                  Arr(  id, ir, IARR)%SrcAngle,  Arr( id, ir, IARR )%RcvrAngle, &
                  Arr( id, ir, IARR )%NTopBnc, Arr( id, ir, IARR )%NBotBnc
          END DO  ! next arrival
       END DO  ! next receiver depth
    END DO  ! next range

    RETURN
  END SUBROUTINE WriteArrivalsASCII

  ! **********************************************************************!

  SUBROUTINE WriteArrivalsBinary( R, Nrd, Nr, TopOpt, freq, SourceType )

    ! Writes the arrival data (amplitude, delay for each eigenray)
    ! Binary output file

    IMPLICIT NONE
    INTEGER           :: Nrd, Nr, ir, id, IArr
    REAL,   PARAMETER :: PI = 3.14159265, RadDeg = 180 / PI
    REAL                 r( Nr ), AArrT
    REAL     (KIND=8) :: freq, f2, alpha, factor
    CHARACTER (LEN=6) :: TopOpt
    CHARACTER (LEN=1) :: SourceType

    ! *** Thorpe attenuation? ***

    IF ( TopOpt(4:4) == 'T' ) THEN
       f2 = ( freq / 1000.0d0 ) ** 2
       ! Original formula from Thorp 1967
       ! alpha = 40.0 * f2 / ( 4100.0 + f2 ) + 0.1 * f2 / ( 1.0 + f2 )   ! dB/kyard
       ! alpha = alpha / 914.4     ! dB / m
       ! alpha = alpha / 8.6858896 ! Nepers / m

       ! Updated formula from JKPS Eq. 1.34
       alpha = 3.3d-3 + 0.11d0 * f2 / ( 1.0d0 + f2 ) + 44.0d0 * f2 / ( 4100.0d0 + f2 ) + 3d-4* f2   ! dB/km
       alpha = alpha / 8685.8896 ! Nepers / m
    ELSE
       alpha = 0.0
    ENDIF

    WRITE( ARRFile ) MAXVAL( NArr( 1 : Nrd, 1 : Nr ) )

    DO id = 1, Nrd
       DO ir = 1, Nr
          IF ( SourceType == 'X' ) THEN   ! line source
             factor =  4.0 * SQRT( PI ) * EXP( -alpha * r( ir ) )
          ELSE                            ! point source
             IF ( r ( ir ) == 0 ) THEN
                factor = 1e5       ! avoid /0 at origin
             ELSE
                factor = EXP( -alpha * r( ir ) ) / SQRT( r( ir ) )  ! vol. atten. and cyl. spreading
             END IF
          END IF

          WRITE( ARRFile ) NArr( id, ir )

          DO IARR = 1, NArr( id, ir )
             AArrT = SNGL( factor * Arr( id, ir, IARR )%A )
             ! integers written out as reals below for fast reading in Matlab
             WRITE( ARRFile ) AArrT, RadDeg * Arr( id, ir, IARR )%Phase, Arr( id, ir, IARR )%delay, &
                  Arr( id, ir, IARR)%SrcAngle, Arr( id, ir, IARR )%RcvrAngle,  &
                  REAL( Arr( id, ir, IARR )%NTopBnc ), REAL( Arr( id, ir, IARR )%NBotBnc )
          END DO   ! next arrival
       END DO   ! next receiver depth
    END DO   ! next range

    RETURN
  END SUBROUTINE WriteArrivalsBinary

END MODULE ArrMod
