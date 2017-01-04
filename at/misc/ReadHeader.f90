SUBROUTINE ReadHeader( SHDFile, FileName, Title, Freq, Atten, DeltaR, PlotType, xs, ys )

  ! Read header from disk file

  ! FileName is a SHDFIL for complex pressure or a GRNFIL for a Green's function
  ! Title  arbitrary title
  ! theta  vector of bearing lines,   theta( 1: Ntheta )
  ! sd     vector of source   depths, sd(1:Nsd)
  ! rd     vector of receiver depths, rd(1:Nrd)
  ! r      vector of receiver ranges, r(1:Nr)
  ! Freq   frequency
  ! Atten  stabilizing attenuation (which is only important for FFP runs. Use zero for other cases.)

  USE SdRdRMod
  IMPLICIT NONE
  INTEGER, PARAMETER                :: PRTFile = 6
  REAL (KIND=4),      INTENT( OUT ) :: DeltaR, Freq, xs, ys, Atten
  INTEGER                           :: IAllocStat, IOStat, LRecL, SHDFile
  CHARACTER (LEN=80)                :: Title, FileName
  CHARACTER (LEN=10), INTENT( OUT ) :: PlotType

  ! Open file, read header
  IF ( SHDFile == 0 ) SHDFile = 25
  IF ( FileName( 1 : 1 ) == ' ' ) FileName = 'SHDFIL'

  ! INQUIRE( FILE = FileName, RECL = IRECL )
  OPEN( UNIT = SHDFile,   FILE = FileName, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4, IOSTAT = IOStat )
  IF ( IOStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadHeader', 'Unable to open shade file' )
  READ( SHDFile, REC = 1 ) LRecl
  CLOSE( UNIT = SHDFile )
  OPEN(  UNIT = SHDFile,   FILE = FileName, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4 * LRecl )

  READ( SHDFile, REC = 1 ) LRecl, Title
  READ( SHDFile, REC = 2 ) PlotType, xs, ys
  READ( SHDFile, REC = 3 ) Freq, Ntheta, Nsd, Nrd, Nr, atten

  ALLOCATE( sd( Nsd ), rd( Nrd ), R( Nr ), theta( Ntheta ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'RDHEAD', 'Too many source/receiver combinations' )

  READ( SHDFile, REC = 4 ) theta
  READ( SHDFile, REC = 5 ) sd
  READ( SHDFile, REC = 6 ) rd
  READ( SHDFile, REC = 7 ) R

  xs = xs / 1000.0   ! convert to km
  ys = ys / 1000.0

  DeltaR = R( NR ) - R( NR - 1 )

END SUBROUTINE ReadHeader
