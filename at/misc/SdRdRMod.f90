MODULE SdRdRMod
  ! Reads in source depths, receiver depths, and receiver ranges

  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER   :: Number_to_Echo = 21
  INTEGER              :: Nsd, Nrd, NR, Ntheta
  INTEGER, ALLOCATABLE :: Isd( : ), Ird( : )
  REAL,    ALLOCATABLE ::  sd( : ),  rd( : ), WS( : ), WR( : ), R( : ), theta( : )

CONTAINS

  SUBROUTINE SDRD( ENVFile, PRTFile, zMin, zMax )

    IMPLICIT NONE
    INTEGER ENVFile, PRTFile, IS, IR, IAllocStat
    REAL    zMin, zMax

    ! *** Read source depths ***

    READ(  ENVFile, * ) Nsd
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of sources   = ', Nsd
    WRITE( PRTFile, * ) 'Source depths (m)'

    IF ( Nsd <= 0 ) THEN
       WRITE( PRTFile, * ) ' Nsd = ', Nsd
       CALL ERROUT( PRTFile, 'F', 'SDRD', 'Number of sources must be positive'  )
    ENDIF

    IF ( ALLOCATED( sd ) ) DEALLOCATE( sd, WS, Isd )
    ALLOCATE( sd( MAX( 3, Nsd ) ), WS( Nsd ), Isd( Nsd ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) THEN
       WRITE( PRTFile, * ) 'Nsd = ', Nsd
       CALL ERROUT( PRTFile, 'F', 'SDRD', 'Too many sources'  )
    END IF

    sd( 3 ) = -999.9
    READ( ENVFile, * ) sd( 1 : Nsd )

    CALL SUBTAB( sd, Nsd )
    !CALL SORT(   sd, Nsd )

    IF ( Nsd >= 1  ) WRITE( PRTFile, "( 5G14.6 )" ) ( sd( IS ), IS = 1, MIN( Nsd, Number_to_Echo ) )
    IF ( Nsd > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', sd( Nsd )

    ! *** Read receiver depths ***

    READ(  ENVFile, * ) Nrd
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of receivers = ', Nrd
    WRITE( PRTFile, * ) 'Receiver depths (m)'

    IF ( Nrd <= 0 ) THEN
       WRITE( PRTFile, * ) ' Nrd = ', Nrd
       CALL ERROUT( PRTFile, 'F', 'SDRD', 'Number of receivers must be positive'  )
    ENDIF

    IF ( ALLOCATED( rd ) ) DEALLOCATE( rd, WR, Ird )
    ALLOCATE( rd( MAX( 3, Nrd ) ), WR( Nrd ), Ird( Nrd ), Stat = IAllocStat  )
    IF ( IAllocStat /= 0 ) THEN
       WRITE( PRTFile, * ) 'Nrd = ', Nrd
       CALL ERROUT( PRTFile, 'F', 'SDRD', 'Too many receivers'  )
    END IF

    rd( 3 ) = -999.9
    READ( ENVFile, * ) rd( 1 : Nrd )

    CALL SUBTAB( rd, Nrd )
    !CALL SORT(   rd, Nrd )

    IF ( Nrd >= 1 ) WRITE( PRTFile, "( 5G14.6 )" ) ( rd( IR ), IR = 1, MIN( Nrd, Number_to_Echo ) )
    IF ( Nrd > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', rd( Nrd )

    ! *** Check for sd/rd in upper or lower halfspace ***
    ! should be converted to use vector operations

    DO IS = 1, Nsd
       IF      ( sd( IS ) < zMin ) THEN
          sd( IS ) = zMin
          CALL ERROUT( PRTFile, 'W', 'SdRdRMod', 'Source above the top bdry has been moved down' ) 
       ELSE IF ( sd( IS ) > zMax ) THEN
          sd( IS ) = zMax
          CALL ERROUT( PRTFile, 'W', 'SdRdRMod', 'Source below the bottom bdry has been moved up' ) 
       ENDIF
    END DO

    DO IR = 1, Nrd
       IF      ( rd( IR ) < zMin ) THEN
          rd( IR ) = zMin
          CALL ERROUT( PRTFile, 'W', 'SdRdRMod', 'Receiver above the top bdry has been moved down' ) 
       ELSE IF ( rd( IR ) > zMax ) THEN
          rd( IR ) = zMax
          CALL ERROUT( PRTFile, 'W', 'SdRdRMod', 'Receiver below the bottom bdry has been moved up' ) 
       ENDIF
    END DO

    RETURN
  END SUBROUTINE SDRD

  !********************************************************************

  SUBROUTINE RANGES( ENVFile, PRTFile )

    ! Read receiver ranges

    IMPLICIT NONE
    INTEGER, PARAMETER   :: Number_to_Echo = 21
    INTEGER                 ENVFile, PRTFile, IR, IAllocStat

    READ(  ENVFile, * ) NR
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of ranges   = ', NR
    WRITE( PRTFile, * ) 'Receiver ranges (km)'

    IF ( ALLOCATED( R ) ) DEALLOCATE( R )
    ALLOCATE( R( MAX( 3, NR ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'RANGES', 'Too many range points' )

    R( 3 ) = -999.9
    READ( ENVFile, * ) R( 1 : NR )

    CALL SUBTAB( R, NR )
    CALL SORT(   R, NR )

    WRITE( PRTFile, "( 5G14.6 )" ) ( R( IR ), IR = 1, MIN( NR, Number_to_Echo ) )
    IF ( NR > Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', R( NR )

    R( 1 : NR ) = 1000.0 * R( 1 : NR )   ! Convert ranges to meters

    ! For a point source can't have receiver at origin
    ! IF ( OPT(1:1) == 'R' .AND. R( 1 ) <= 0.0 ) 
    ! IF ( R( 1 ) <= 0.0 ) R( 1 ) = MIN( 1.0, R( 2 ) )

    RETURN
  END SUBROUTINE ranges

END MODULE SdRdRMod
