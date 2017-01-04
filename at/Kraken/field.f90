PROGRAM FIELD

  ! Generate file of replica vectors
  ! Useage: Field FileRoot
  ! where
  ! FileRoot.mod contains the modes and
  ! FileRoot.shd contains the output shade file

  USE SdRdRMod

  IMPLICIT NONE
  INTEGER,   PARAMETER :: FLPFile = 5, PRTFile = 6, ShdFile = 25, MaxM = 20000 ! MaxM also in Eval, EvalAD, EvalCM
  INTEGER                 I, IR, iProf, NProf, Nrr, IAllocStat, iRec, IS, M, MLimit, MSrc
  REAL                    zMin, zMax, freq, atten
  COMPLEX                 k( MaxM )
  CHARACTER   (LEN=50) :: Opt
  CHARACTER   (LEN=80) :: SHDTitle, Title, FileRoot
  CHARACTER   (LEN=1 ) :: Comp
  CHARACTER   (LEN=10) :: PlotType
  REAL,    ALLOCATABLE :: rr( : ), rProf( : )
  COMPLEX, ALLOCATABLE :: phiS( :, : ), phiR( :, : ), P( :, : ), C( : ), Ptemp( : )

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) '__________________________________________________'
  WRITE( PRTFile, * ) 'Running FIELD'
  WRITE( PRTFile, * ) 'Sums modes, producing pressure'
  WRITE( PRTFile, * )

  SHDTitle( 1 : 1 ) = '$'

  READ( FLPFile, * ) SHDTitle
  READ( FLPFile, * ) Opt
  READ( FLPFile, * ) MLimit

  IF ( NProf == 1      ) THEN
     WRITE( PRTFile, * ) 'Range-independent calculation'
  ELSE
     WRITE( PRTFile, * ) 'Range-dependent calculation'
     IF ( Opt( 2 : 2 ) == 'C' ) THEN
        WRITE( PRTFile, * ) 'Coupled modes'
     ELSE
        WRITE( PRTFile, * ) 'Adiabatic modes'
     END IF
  END IF

  Comp   = Opt( 3 : 3 )
  MLimit = MIN( MaxM, MLimit )

  ! *** Read profile ranges ***

  READ( FLPFile, * ) NProf
  ALLOCATE( rProf( MAX( 3, NProf + 1 ) ), Stat = IAllocStat )   ! NProf + 1 profiles (one added at `infinite' range)

  IF ( IAllocStat /= 0 ) THEN
     WRITE( PRTFile, * ) 'NProf = ', NProf
     CALL ERROUT( PRTFile, 'F', 'FIELD', 'Dynamic memory allocation failed: Too many profiles' )
  ENDIF

  IF ( NProf > 2 ) rProf( 3 ) = -999.9
  READ( FLPFile, * ) rProf( 1 : NProf )
  CALL SUBTAB( rProf, NProf )

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Number of profiles   = ', NProf
  IF ( NProf >= 1  ) WRITE( PRTFile, "( 5G14.6 )" ) ( rProf( iProf ), iProf = 1, MIN( NProf, Number_to_Echo ) )
  IF ( NProf > Number_to_Echo  ) WRITE( PRTFile, * ) ' ... ', rProf( NProf )

  ! EVALAD/EVALCM need a profile at zero range
  IF ( rProf( 1 ) /= 0.0 ) CALL ERROUT( PRTFile, 'F', 'FIELD', 'The first profile must be at 0 km' )

  CALL RANGES( FLPFile, PRTFile )           ! Read receiver ranges
  zMin = -HUGE( zMin )
  zMax = +HUGE( zMax )
  CALL SDRD( FLPFile, PRTFile, zMin, zMax ) ! Read source/receiver depths

  ALLOCATE( phiS( MaxM, Nsd ), phiR( MaxM, Nrd ), C( MaxM ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) THEN
     WRITE( PRTFile, * ) 'Nrd = ', Nrd
     CALL ERROUT( PRTFile, 'F', 'FIELD', 'Dynamic memory allocation failed: Too many receiver depths' )
  ENDIF

  ! *** Read receiver ranges (offsets from vertical) ***

  READ( FLPFile, * ) Nrr

  IF ( Nrr /= Nrd ) THEN
     WRITE( PRTFile, * ) 'Nrr, Nrd = ', Nrr, Nrd
     CALL ERROUT( PRTFile, 'W', 'FIELD', 'Nrr being set to Nrd' )
     Nrr = Nrd
  ENDIF

  ALLOCATE( rr( Nrr ) )
  IF ( Nrr > 1 ) rr( 2 ) = -999.9
  IF ( Nrr > 2 ) rr( 3 ) = -999.9
  READ( FLPFile, * ) rr( 1 : Nrr )

  CALL SUBTAB( rr, Nrd )
  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Number of receiver range offsets = ', Nrd
  WRITE( PRTFile, * ) 'Receiver range offsets (km)'

  IF ( Nrd >= 1 ) WRITE( PRTFile,  "( 5G14.6 )" ) ( rr( IR ), IR = 1, MIN( Nrd, Number_to_Echo ) )
  IF ( Nrd > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', rr( Nrd )

  WRITE( PRTFile, * )

  ALLOCATE ( P( Nrd, Nr ), Ptemp( Nr ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) STOP "Fatal Error: Insufficient memory to allocate P( Nrd, Nr )"

  !  *** Read in modes ***

  IProf = 1
  CALL GetModes( FileRoot, iProf, MaxM, sd, Nsd, 'N' , k, phiS, MSrc, freq, Title )
  CALL GetModes( FileRoot, iProf, MaxM, rd, Nrd, Comp, k, phiR, MSrc, freq, Title )

  ! Generate header

  IF ( SHDTitle( 1 : 1 ) == '$' ) SHDTitle = Title
  PlotType   = '          '
  atten      = 0.0
  ALLOCATE( theta( 1 ) )
  theta( 1 ) = 0   ! dummy bearing angle
  Ntheta     = 1
  CALL WriteHeader( TRIM( FileRoot ) // '.shd', SHDTitle, theta, Ntheta, sd, Nsd, rd, Nrd, R, Nr, freq, atten, PlotType, 0.0, 0.0 )
  iRec = 7

  ! *** MAIN LOOP: For each source evaluate and write the field ***

  DO IS = 1, Nsd
     M = MIN( MLimit, MSrc )   ! Set number of propagating modes

     C( 1 : MSrc ) = phiS( 1 : MSrc, is )
     IF    ( NProf == 1      ) THEN   ! Range-independent case
        CALL EVAL(                              C, phiR,     Nrd, R, Nr, rr, k, M, Opt, P )
     ELSE
        IF ( Opt(2:2) == 'C' ) THEN   ! Range-  dependent case
           ! Coupled mode theory
           CALL EVALCM( FileRoot, rProf, NProf, C, phiR, rd, Nrd, R, Nr,     k, M, Opt, P )
        ELSE
           ! Adiabatic mode theory
           CALL EVALAD( FileRoot, rProf, NProf, C, phiR, rd, Nrd, R, Nr,        M, Opt, P )
        ENDIF
     ENDIF

     ! write out the field
     DO I = 1, Nrd
        IRec  = IRec + 1
        Ptemp = P( i, 1 : Nr )   ! temporary variable, avoids compiler warning
        WRITE( SHDFile, REC = IRec ) Ptemp
     END DO

  END DO   ! next source depth
  CLOSE( ShdFile )

END PROGRAM FIELD
