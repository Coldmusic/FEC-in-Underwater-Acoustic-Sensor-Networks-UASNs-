MODULE beampatternmod

  ! Loads a source beam pattern

  SAVE
  INTEGER, PARAMETER         :: SBPFile = 50
  INTEGER                    :: NSBPPts
  REAL (KIND=8), ALLOCATABLE :: SrcBmPat( :, : )

CONTAINS

  SUBROUTINE ReadPAT( FileRoot, SBP, PRTFile )

    IMPLICIT NONE
    INTEGER            :: PRTFile
    INTEGER            :: I, IAllocStat, IOStat
    CHARACTER (LEN=1)  :: SBP
    CHARACTER (LEN=80) :: FileRoot

    IF ( SBP == '*' ) THEN
       WRITE( PRTFile, * ) '*********************************'
       WRITE( PRTFile, * ) 'Using source beam pattern file'

       OPEN( UNIT = SBPFile,   FILE = TRIM( FileRoot ) // '.sbp', STATUS = 'OLD', IOSTAT = IOStat )
       IF ( IOstat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP-ReadPat', 'Unable to open source beampattern file' )
 
       READ(  SBPFile, * ) NSBPPts
       WRITE( PRTFile, * ) 'Number of source beam pattern points', NSBPPts

       ALLOCATE( SrcBmPat( NSBPPts, 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
            CALL ERROUT( PRTFile, 'F', 'BELLHOP-ReadPat', 'Insufficient memory for source beam pattern data: reduce # SBP points' )

       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) ' Angle (degrees)  Power (dB)'

       DO I = 1, NSBPPts
          READ(  SBPFile, * ) SrcBmPat( I, : )
          WRITE( PRTFile, FMT = "( 2G11.3 )" ) SrcBmPat( I, : )
       END DO

    ELSE   ! no pattern given, use omni source pattern
       NSBPPts = 2
       ALLOCATE( SrcBmPat( 2, 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP-ReadPat', 'Insufficient memory'  )
       SrcBmPat( 1, : ) = [ -180.0, 0.0 ]
       SrcBmPat( 2, : ) = [  180.0, 0.0 ]
    ENDIF

    SrcBmPat( :, 2 ) = 10 ** ( SrcBmPat( :, 2 ) / 20 )  ! convert dB to linear scale

  END SUBROUTINE ReadPAT

END MODULE beampatternmod
