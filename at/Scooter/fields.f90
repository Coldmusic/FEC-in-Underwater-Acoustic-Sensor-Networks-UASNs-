PROGRAM FIELDS

  !  Compute pressure from Green's function.

  ! The transform parameters are chosen to satisfy certain sampling requirements.

  ! It is assumed that the user chose the maximum possible 
  ! spacing in computing g(k) and so deltak is not allowed to be any larger.

  ! The other part of the kernel is the exp(ikr) term.
  ! This must be sampled on the k-axis finely enough to have about 6 pts per wavelength.

  ! deltaR implies kmax. kmax is the upper limit of integration and beyond the last point computed.
  ! This assumes the user had chosen kmax as small as possible while still covering the support of G(k).

  IMPLICIT NONE
  INTEGER, PARAMETER   :: FLPFile = 5, PRTFile = 6, MaxNk =  1000000, MaxNsd = 501, MaxNrd = 2001
  REAL,    PARAMETER   :: pi = 3.14159265
  INTEGER              :: IAllocStat, IRatio, IRatiodeltar, J, Nk, Nrd, Nsd, Nrr, NrrLast, Nrrsubsample, Nt, Nt2
  REAL                 :: k( MaxNk ), rd( MaxNrd ), sd( MaxNsd ), Atten, AttInt, Freq, kmax, deltak, deltakInterp, &
                          Rmin, Rmax, RMinKM, RMaxKM, deltar
  COMPLEX              :: G( MaxNk )
  CHARACTER   (LEN=80) :: PlotTitle, FileRoot
  CHARACTER   (LEN=3 ) :: Option
  REAL,    ALLOCATABLE :: R( : ), kInterp( : )
  COMPLEX, ALLOCATABLE :: Ginterp( : ), P( : ), Temp( : )

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) '__________________________________________________'
  WRITE( PRTFile, * ) 'Running FIELDS'
  WRITE( PRTFile, * ) 'Transforms Green''s function FILE, producing pressure'
  WRITE( PRTFile, * )

  ! Begin by reading in the data
  READ( FLPFile, * ) Option
  SELECT CASE ( Option( 1 : 1 ) )
  CASE ( 'R' )
     WRITE( PRTFile, * ) 'Cylindrical coordinates/point source'
  CASE ( 'X' )
     WRITE( PRTFile, * ) 'Cartesian coordinates/line source'
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Option( 1 : 1 ) should select Line vs. point source' )
  END SELECT

  SELECT CASE ( Option( 2 : 2 ) )
  CASE ( 'P' )
     WRITE( PRTFile, * ) 'Using only positive part of wavenumber spectrum'
  CASE ( 'N' )
     WRITE( PRTFile, * ) 'Using only negative part of wavenumber spectrum'
  CASE ( 'B' )
     WRITE( PRTFile, * ) 'Using both positive and negative part of wavenumber spectrum'
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Option(2:2) should select positive, negative, or both parts of the spectrum' )
  END SELECT

  SELECT CASE ( Option( 3 : 3 ) )
  CASE ( 'O' )
     WRITE( PRTFile, * ) 'Polynomial interpolation'
  CASE ( 'P' )
     WRITE( PRTFile, * ) 'Pade interpolation'
  CASE DEFAULT
     Option( 3 : 3 ) = 'O'   ! use Polynomial interpolation as default
     !CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Option(3:3) should select Polynomial (''O'') or Pade (''P'') interpolation' )
  END SELECT

  READ( FLPFile, * ) RminKM, RmaxKM, Nrr
  Rmin   = 1000.0 * RminKM   ! convert km to m
  Rmax   = 1000.0 * RmaxKM
  deltaR = ( Rmax - Rmin ) / ( Nrr - 1 )

  ! Read the header records from GRNFile
  CALL RDHead( TRIM( FileRoot ) // '.grn', PlotTitle, sd, Nsd, MaxNsd, rd, Nrd, MaxNrd, &
       k, Nk, MaxNk, deltak, Atten, Freq )
  ! Set up for transform: need deltak, NT
  ! deltak is what scooter used; deltakInterp is for interpolation
  ! If deltaR is too big, take submultiple

  kmax = k( 1 ) + 2.0 * pi / deltaR
  IRatiodeltar = 1
  IF ( kmax < k( Nk ) ) THEN
     IRatiodeltar = INT( ( k( Nk ) - k( 1 ) ) / ( kmax - k( 1 ) ) ) + 1
     deltaR       = deltaR / IRatiodeltar
     Nrr          = IRatiodeltar * ( Nrr - 1 ) + Nrr
     kmax         = k( 1 ) + 2 * pi / deltaR
     write( PRTfile, * ) 'IRatiodeltar', IRatiodeltar
     WRITE( PRTFile, * ) 'Number or ranges, Nrr, increased so that wavenumber limit exceeds kmax used by SCOOTER', Nrr
  END IF

  ! Compute NT based on deltak (= 1/Rmax) sampling requirement
  NT2 = NINT( Rmax * ( kmax - k( 1 ) ) / ( 2 * pi ) )       ! deltak limit
  IF ( NT2 > Nrr ) THEN
     WRITE( PRTFile, * ) 'NT bumped to', NT2
     WRITE( PRTFile, * ) 'Thus we are zero filling the wavenumber spectrum to effectively interpolate on a finer grid'
  END IF
  NT = MAX( Nrr, NT2 )

  ! bump Nt if necessary to make sure deltakInterp is not coarser than deltak grid
  deltak       = ( k( Nk ) - k( 1 ) ) / Nk;
  deltakInterp = 2 * pi / ( Nt * deltar );
  IF ( deltakInterp > deltak ) THEN
     IRatio = NINT( deltakInterp / deltak );
     Nt     = IRatio * Nt;
     WRITE( PRTFile, * ) 'Transform size, Nt, bumped to ensure deltak sampling is fine enough', Nt
  END IF

  ! Parameters when N must be a power of 2
  NT = 2 ** ( INT( LOG10( REAL( NT ) ) / 0.301 ) + 1 )
  deltakInterp = 2 * pi / ( NT * deltaR )

  WRITE( PRTFile, * ) 'NT  used = ', NT
  WRITE( PRTFile, * ) 'Nrr used = ', Nrr
  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'deltak       = ', deltak
  WRITE( PRTFile, * ) 'deltakInterp = ', deltakInterp
  WRITE( PRTFile, * ) 'deltaR       = ', deltaR
  WRITE( PRTFile, * )

  ALLOCATE( kInterp( NT ), Ginterp( NT ), P( NT ), Temp( 2 * NT ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Insufficient memory to allocate kInterp, Ginterp, ...' )

  ! Compute optimal stabilizing attenuation.
  ! Atten is the value used by SCOOTER
  ! AttInt is the optimal value for the interpolated k-values

  AttInt = ( deltakInterp / deltak ) * Atten

  IF ( PlotTitle( 1 : 5 ) == 'SPARC' ) THEN
     WRITE( PRTFile, * ) 'SPARC RUN: AttInt SET TO ZERO'
     AttInt = 0.0
     IF ( Atten /= 0.0 ) CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Stabilizing attenuation must vanish' )
  ENDIF

  kInterp = k( 1 ) +[ ( J, J=0, Nt - 1 ) ] * deltakInterp   ! Set up vector of kInterp points
  kmax    = kInterp( Nt ) + deltakInterp   ! because the fft goes from 0 to kmax but G(kmax) is not computed
  deltar  = 2 * pi / ( kmax - kInterp( 1 ) )
  Nrr     = Nt
  NrrLast = min( NINT( ( Rmax - Rmin ) / deltar )+1, Nt )

  ! set up vector of range points (subsampled to satisfy user request)
  Nrrsubsample = ( NrrLast - 1 ) / Iratiodeltar + 1
  ALLOCATE( R( Nrrsubsample ) )
  R = Rmin + [ ( J, J=0, Nrrsubsample - 1 ) ] * Iratiodeltar * deltar

  ! Construct shade file
  CALL SHADE( FileRoot, k, Atten, Nk, G, NT, kInterp, deltakInterp, Ginterp, &
       &  R, Nrr, NrrLast, Nrrsubsample, Iratiodeltar, deltaR, &
       &  sd, Nsd, rd, Nrd, PlotTitle, Freq, AttInt, Option, P, Temp )

END PROGRAM FIELDS

!**********************************************************************C

SUBROUTINE RDHead( FileName, PlotTitle, sd, Nsd, MaxNsd, rd, Nrd, MaxNrd, k, Nk, MaxNk, deltak, Atten, Freq )

  ! Routine to read header from disk file
  ! This routine is essentially the same as at/misc/ReadHeader except that the range vector is replaced
  ! by a wavenumber vector.

  IMPLICIT NONE
  INTEGER,            PARAMETER   :: GRNFile = 20, PRTFile = 6
  INTEGER,            INTENT(IN ) :: MaxNsd, MaxNrd, MaxNk
  INTEGER,            INTENT(OUT) :: Nsd, Nrd, Nk
  INTEGER                         :: Ntheta, IOStat, LRecL
  REAL,               INTENT(OUT) :: k( MaxNk ), rd( MaxNrd ), sd( MaxNsd ), Atten, Freq, deltak
  CHARACTER (LEN=* ), INTENT(IN ) :: FileName
  CHARACTER (LEN=80), INTENT(OUT) :: PlotTitle

  OPEN( FILE = FileName, UNIT = GRNFile, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 10, IOSTAT = IOStat )
  IF ( IOStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'FIELDS:RDHead', 'Unable to open GRNFile' )
  READ( GRNFile, REC = 1 ) LRECL
  CLOSE( GRNFile )
  OPEN( FILE = FileName, UNIT = GRNFile, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4 * LRECL )

  ! Read data
  READ( GRNFile, REC = 1 ) LRECL, PlotTitle
  !READ( GRNFile, REC = 2 ) PlotType, XS, YS
  READ( GRNFile, REC = 3 ) Freq, Ntheta, Nsd, Nrd, Nk, Atten

  IF ( Nsd > MaxNsd ) CALL ERROUT( PRTFile, 'F', 'FIELDS:RDHead', 'Too many source depths  ' )
  IF ( Nrd > MaxNrd ) CALL ERROUT( PRTFile, 'F', 'FIELDS:RDHead', 'Too many receiver depths' )
  IF ( Nk  > MaxNk  ) CALL ERROUT( PRTFile, 'F', 'FIELDS:RDHead', 'Too many k-space points ' )

  READ( GRNFile, REC = 5 ) sd( 1 : Nsd )
  READ( GRNFile, REC = 6 ) rd( 1 : Nrd )
  READ( GRNFile, REC = 7 ) k(  1 : Nk  )  ! vector of k-space points
  deltak = k( 2 ) - k( 1 )

END SUBROUTINE RDHead
  !**********************************************************************C
  SUBROUTINE SHADE( FileRoot, k, Atten, Nk, G, NT, kInterp, deltakInterp, Ginterp, &
       &   R, Nrr, NrrLast, Nrrsubsample, Iratiodeltar, deltaR, &
       &   sd, Nsd, rd, Nrd, PlotTitle, Freq, AttInt, Option, P, Temp )

    ! Performs the transforms to convert the Green's function file to a shade file
    ! Expects
    !    k, Nk: The wavenumber data
    !    sd, rd: Source depth, receiver depth data
    ! Returns
    !    Nothing

    IMPLICIT NONE
    INTEGER, PARAMETER :: PRTFile = 6, GRNFile = 20, SHDFile = 25
    INTEGER            :: Isd, Ird, ISR, IRec, Nsd, Nrd, Nk, NT, Nrr, Ntheta, NrrLast, Nrrsubsample, Iratiodeltar
    REAL               :: R( Nrrsubsample ), kInterp( NT ), theta( 1 ), sd( Nsd ), rd( Nrd ), k( Nk ), &
                          Rmin, deltakInterp, Freq, AttInt, Atten, deltar
    COMPLEX            :: G( Nk ), Ginterp( NT ), P( Nrr ), Temp( NT )
    CHARACTER (LEN=80) :: PlotTitle, FileRoot
    CHARACTER (LEN=3 ) :: Option
    CHARACTER (LEN=10) :: PlotType

    theta( 1 ) = 0.0   ! dummy bearing anle
    Ntheta     = 1
    Rmin       = R( 1 )

    CALL WriteHeader( TRIM( FileRoot ) // '.shd', PlotTitle, theta, Ntheta, sd, Nsd, rd, Nrd, &
         &  R, Nrrsubsample, Freq, Atten, PlotType, 0.0, 0.0 )

    DO Isd = 1, Nsd   ! Loop over all source/rcvr combos
       WRITE( PRTFile, * ) 'Transform for source depth: ', sd( Isd )

       DO Ird  = 1, Nrd
          ISR  = ( Isd - 1 ) * Nrd + Ird   ! Index of source/receiver
          IRec = 7 + ISR
          READ( GRNFile, REC = IRec ) G( 1 : Nk )
          CALL INTERP( k, Atten, AttInt, G, Nk, kInterp, Ginterp, NT, Option )   ! Interpolate

          ! Evaluate either Fourier or Hankel transform
          IF ( Option(1:1) == 'X' ) THEN
             CALL FTS( NT, kInterp( 1 ), deltakInterp, AttInt, Rmin, Nrr, deltaR, Ginterp, Temp, P )
          ELSE
             CALL HTS( NT, kInterp( 1 ), deltakInterp, AttInt, Rmin, Nrr, deltaR, Ginterp, Temp, P, Option )
          ENDIF

          WRITE( SHDFile, REC = IRec ) P( 1 : NrrLast : Iratiodeltar )   ! Write out the field

       END DO   ! next receiver depth
    END DO   ! next source depth

  END SUBROUTINE SHADE

  !**********************************************************************C

  SUBROUTINE INTERP( k, Atten, AttInt, G, Nk, kInterp, Ginterp, NT, Option )

    ! Produces an evenly sampled kernel from input G

    IMPLICIT NONE
    INTEGER, PARAMETER   :: PRTFile = 6, ISize = 3
    COMPLEX, PARAMETER   :: i = ( 0.0, 1.0 )
    INTEGER                 ICenter, IRight, INew, Nk, Nt, It
    REAL,    INTENT(IN)  :: k( Nk ), kInterp( Nt ), Atten, AttInt
    COMPLEX              :: xinterp, X( ISize ), F( ISize ), PADE, PC
    COMPLEX, INTENT(IN)  :: G( Nk )
    COMPLEX, INTENT(OUT) :: Ginterp( Nt )
    CHARACTER   (LEN=80) :: ErrorMessage
    CHARACTER   (LEN=3 ) :: Option

    ! Initialize interpolation data
    ICenter = ISize / 2 + 1   ! counter to mark center of abscissas used for interpolation
    IRight  = ISize
    INew    = ISize

    X( 1 : ISize ) = ( k( 1 : ISize ) + i * Atten )**2   ! abscissa
    F( 1 : ISize ) = G( 1 : ISize )                      ! ordinate

    DO It = 1, NT ! Main loop: step through values in kInterp

       ! Desirable/possible to advance interpolation window?
       DO WHILE ( kInterp( It ) > k( ICenter ) .AND. IRight < Nk )
          ICenter = ICenter + 1
          IRight  = IRight  + 1
          INew    = MOD( INew, ISize ) + 1   ! this is the next open slot

          X( INew ) = ( k( IRight ) + i * Atten ) ** 2
          F( INew ) = G( IRight )
       END DO

       ! Interpolate or zero fill
       IF ( kInterp( It ) <= k( Nk ) ) THEN
          xinterp = ( kInterp( It ) + i * AttInt ) ** 2
          IF ( Option( 3 : 3 ) == 'O' ) THEN
             Ginterp( It ) = PC(   xinterp, X, F, ISize, ErrorMessage )   ! polynomial
          ELSE
             Ginterp( It ) = Pade( xinterp, X, F, ISize, ErrorMessage )   ! Pade
          ENDIF

          IF ( ErrorMessage( 1 : 6 ) /= '      ' ) WRITE( PRTFile, * ) ErrorMessage
       ELSE
          Ginterp( It ) = 0.0
       ENDIF
    END DO

  END SUBROUTINE INTERP
