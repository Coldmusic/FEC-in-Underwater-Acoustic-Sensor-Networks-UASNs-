SUBROUTINE WriteHeader( FileName, Title, theta, Ntheta, sd, Nsd, rd, Nrd, r, Nr, Freq, Atten, PlotType, XS, YS )

  ! Write header to disk file

  ! FileName is a SHDFIL for complex pressure or a GRNFIL for a Green's function
  ! Title  arbitrary title
  ! theta  vector of bearing lines,   theta( 1: Ntheta )
  ! sd     vector of source   depths, sd(1:Nsd)
  ! rd     vector of receiver depths, rd(1:Nrd)
  ! r      vector of receiver ranges, r(1:Nr)
  ! Freq   frequency
  ! Atten  stabilizing attenuation (which is only important for FFP runs. Use zero for other cases.)

  IMPLICIT NONE
  INTEGER, PARAMETER          :: SHDFile = 25
  INTEGER,       INTENT( IN ) :: Ntheta, Nsd, Nrd, Nr
  REAL (KIND=4), INTENT( IN ) :: theta( Ntheta ), r( Nr ), rd( Nrd ), sd( Nsd ), Freq, XS, YS, Atten
  CHARACTER,     INTENT( IN ) :: FileName*( * ), Title*( * ), PlotType*( 10 )

  INTEGER LRecl

  LRecl = MAX( 40, Ntheta, Nsd, Nrd, 2 * Nr )   ! words/record

  OPEN ( FILE = FileName, UNIT = SHDFile, STATUS = 'UNKNOWN', ACCESS = 'DIRECT', RECL = 4 * LRecl, FORM = 'UNFORMATTED')
  WRITE( SHDFile, REC = 1 ) LRecl, Title( 1 : 80 )
  WRITE( SHDFile, REC = 2 ) PlotType, XS, YS
  WRITE( SHDFile, REC = 3 ) Freq, Ntheta, Nsd, Nrd, Nr, atten
  WRITE( SHDFile, REC = 4 ) theta
  WRITE( SHDFile, REC = 5 ) sd
  WRITE( SHDFile, REC = 6 ) rd
  WRITE( SHDFile, REC = 7 ) r

END SUBROUTINE WriteHeader

!**********************************************************************!

SUBROUTINE WriteField( P, Nrd, Nr, IRec )

  ! Write the field to disk

  IMPLICIT NONE
  INTEGER, PARAMETER    :: SHDFile = 25
  INTEGER               :: IREC, I, Nrd, Nr
  COMPLEX, INTENT( IN ) :: P( Nrd, Nr )

  DO I = 1, Nrd
     IRec = IRec + 1
     WRITE( SHDFile, REC = IRec ) P( I, : )
  END DO

END SUBROUTINE WriteField
