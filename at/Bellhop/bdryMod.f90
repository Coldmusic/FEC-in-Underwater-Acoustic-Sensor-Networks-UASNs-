MODULE bdrymod

  ! Loads
  ! altimetry (top bdry) and bathymetry (bottom bdry) data
  !
  ! x = coordinate of boundary
  ! t = tangent for a segment
  ! n = normal  for a segment
  ! Len = length of tangent (temporary variable to normalize tangent)
  ! tangents are also calculated for the nodes, if the curvilinear option is used for the boundaries.

  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: ATIFile = 40, BTYFile = 41, Number_to_Echo = 21
  INTEGER            :: NATIPTS, NBTYPTS, ii, I, IOStat, IAllocStat
  CHARACTER  (LEN=1) :: atiType, btyType

  TYPE BdryPt
     REAL (KIND=8) :: x( 2 ), t( 2 ), n( 2 ), Nodet( 2 ), Noden( 2 ), Len, Kappa
  END TYPE

  TYPE(BdryPt), ALLOCATABLE :: Bot( : ), Top( : )

CONTAINS

  SUBROUTINE ReadATI( FileRoot, TopATI, DepthT, rBox, PRTFile )

    IMPLICIT NONE
    INTEGER                       PRTFile
    CHARACTER          (LEN=1) :: TopATI
    REAL (KIND=8)              :: DepthT, rBox, sss
    REAL (KIND=8), ALLOCATABLE :: phi( : )
    CHARACTER (LEN=80)         :: FileRoot

    IF ( TopATI == '*' ) THEN
       WRITE( PRTFile, * ) '*********************************'
       WRITE( PRTFile, * ) 'Using top-altimetry file'

       OPEN( UNIT = ATIFile,   FILE = TRIM( FileRoot ) // '.ati', STATUS = 'OLD', IOSTAT = IOStat )
       IF ( IOsTAT /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadATI', 'Unable to open altimetry file' )

       READ(  ATIFile, * ) atiType
       SELECT CASE ( atiType )
       CASE ( 'C' )
          WRITE( PRTFile, * ) 'Curvilinear Interpolation'
       CASE ( 'L' )
          WRITE( PRTFile, * ) 'Piecewise linear interpolation'
       CASE DEFAULT
          CALL ERROUT( PRTFile, 'F', 'ReadATI', 'Unknown option for selecting altimetry interpolation' )
       END SELECT

       READ(  ATIFile, * ) NatiPts
       WRITE( PRTFile, * ) 'Number of altimetry points', NatiPts
       NatiPts = NatiPts + 2   ! we'll be extending the altimetry to infinity to the left and right

       ALLOCATE( Top(  NatiPts ), phi( NatiPts ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
            CALL ERROUT( PRTFile, 'F', 'BELLHOP:ReadATI', 'Insufficient memory for altimetry data: reduce # ati points' )

       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) ' Range (km)  Depth (m)'

       DO I = 2, NatiPts - 1
          READ(  ATIFile, * ) Top( I )%x
          IF ( I < Number_to_Echo .OR. I == NatiPts ) THEN   ! echo some values
             WRITE( PRTFile, FMT = "(2G11.3)" ) Top( I )%x
          END IF
          IF ( Top( I )%x( 2 ) < DepthT ) THEN
             CALL ERROUT( PRTFile, 'F', 'BELLHOP:ReadATI', 'Altimetry rises above highest point in the sound speed profile' )
          END IF
       END DO

       CLOSE( ATIFile )

       Top( : )%x( 1 ) = 1000.0 * Top( : )%x( 1 )   ! Convert ranges in km to m

       ! extend the altimetry to +/- infinity in a piecewise constant fashion

       Top( 1 )%x( 1 )       = -sqrt( huge( Top( 1 )%x( 1 ) ) ) / 1.0d5
       Top( 1 )%x( 2 )       = Top( 2 )%x( 2 )
       Top( NatiPts )%x( 1 ) = +sqrt( huge( Top( 1 )%x( 1 ) ) ) / 1.0d5
       Top( NatiPts )%x( 2 ) = Top( NatiPts - 1 )%x( 2 )

       ! compute tangent and outward-pointing normal to top
       ! tToP( 1, : ) = xTop( 1, 2 : NatiPts ) - xTop( 1, 1 : NatiPts - 1 )
       ! tTop( 2, : ) = xTop( 2, 2 : NatiPts ) - xTop( 2, 1 : NatiPts - 1 )

       DO ii = 1, NatiPts - 1
          Top( ii )%t = Top( ii+1 )%x - Top( ii )%x
          Top( ii )%Len = SQRT( Top( ii )%t( 1 ) ** 2 + Top( ii )%t( 2 ) ** 2 )
          Top( ii )%t = Top( ii )%t / Top( ii )%Len

          Top( ii )%n( 1 ) =  Top( ii )%t( 2 )
          Top( ii )%n( 2 ) = -Top( ii )%t( 1 )
       END DO

       IF ( atiType == 'C' ) THEN ! curvilinear option: compute tangent and normal at node by averaging normals on adjacent segments

          DO ii = 2, NatiPts - 1
 
             sss = Top( ii - 1 )%Len / ( Top( ii - 1 )%Len + Top( ii )%Len )
             sss = 0.5
             Top( ii )%Nodet = ( 1.0 - sss ) * Top( ii - 1 )%t + sss * Top( ii )%t

             !tTopNode( :, ii ) = ( LenTop( ii ) * tTop( :, ii - 1 ) + LenTop( ii - 1 ) * tTop( :, ii ) ) / &
             !                    ( LenTop( ii )                     + LenTop( ii - 1 ) )
          END DO

          Top( 1       )%Nodet = [ 1.0, 0.0 ]   ! tangent left-end  node
          Top( NatiPts )%Nodet = [ 1.0, 0.0 ]   ! tangent right-end node
          Top( : )%Noden( 1 )  = +Top( : )%Nodet( 2 )
          Top( : )%Noden( 2 )  = -Top( : )%Nodet( 1 )

          ! compute curvature in each segment
          phi = atan2( Top( : )%Nodet( 2 ), Top( : )%Nodet( 1 ) )
          DO ii = 1, NatiPts - 2
             Top( ii )%kappa = ( phi( ii + 1 ) - phi( ii ) ) / Top( ii )%Len ! this is curvature = dphi/ds
          END DO 

       ELSE
          Top%kappa = 0
       END IF

    ELSE   ! no bathymetry given, use SSP depth for flat top
       NatiPts = 2
       atiType = 'L'
       ALLOCATE( Top( 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP', 'Insufficient memory'  )
       Top( 1 )%x = [ -rBox, DepthT ]
       Top( 2 )%x = [  rBox, DepthT ]

       Top( 1 )%t = [ 1.0,  0.0 ]   ! tangent to top
       Top( 1 )%n = [ 0.0, -1.0 ]   ! outward-pointing normal
       Top( 2 )%t = [ 1.0,  0.0 ]   ! tangent to top
       Top( 2 )%n = [ 0.0, -1.0 ]   ! outward-pointing normal

       Top%Len = sqrt( Top( 1 )%t( 1 )** 2 + Top( 1 )%t( 2 )** 2 )

       Top( 1 )%kappa = 0

    ENDIF

  END SUBROUTINE ReadATI

  ! **********************************************************************!

  SUBROUTINE ReadBTY( FileRoot, BotBTY, DepthB, rBox, PRTFile )

    ! Reads in the bottom bathymetry

    IMPLICIT NONE
    INTEGER                       PRTFile
    CHARACTER          (LEN=1) :: BotBTY
    REAL              (KIND=8) :: DepthB, rBox, sss
    REAL (KIND=8), ALLOCATABLE :: phi( : )
    CHARACTER         (LEN=80) :: FileRoot

    IF ( BotBTY == '*' ) THEN
       WRITE( PRTFile, * ) '*********************************'
       WRITE( PRTFile, * ) 'Using bottom-bathymetry file'

       OPEN( UNIT = BTYFile,   FILE = TRIM( FileRoot ) // '.bty', STATUS = 'OLD', IOSTAT = IOStat )
       IF ( IOsTAT /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadBTY', 'Unable to open bathymetry file' )
 
       READ( BTYFile, * ) btyType
 
       SELECT CASE ( btyType )
       CASE ( 'C' )
          WRITE( PRTFile, * ) 'Curvilinear Interpolation'
       CASE ( 'L' )
          WRITE( PRTFile, * ) 'Piecewise linear interpolation'
       CASE DEFAULT
          CALL ERROUT( PRTFile, 'F', 'ReadBTY', 'Unknown option for selecting bathymetry interpolation' )
       END SELECT

       READ(  BTYFile, * ) NbtyPts
       WRITE( PRTFile, * ) 'Number of bathymetry points', NbtyPts

       NbtyPts = NbtyPts + 2   ! we'll be extending the bathymetry to infinity on both sides
       ALLOCATE( Bot( NbtyPts ), phi( NbtyPts ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
            CALL ERROUT( PRTFile, 'F', 'BELLHOP:ReadBTY', 'Insufficient memory for bathymetry data: reduce # bty points' )

       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) ' Range (km)  Depth (m)'

       DO I = 2, NbtyPts - 1
          READ(  BTYFile, * ) Bot( I )%x
          IF ( I < Number_to_Echo .OR. I == NbtyPts ) THEN   ! echo some values
             WRITE( PRTFile, FMT = "(2G11.3)" ) Bot( I )%x
          END IF
          IF ( Bot( I )%x( 2 ) > DepthB ) THEN
             CALL ERROUT( PRTFile, 'F', 'BELLHOP:ReadBTY', 'Bathymetry drops below lowest point in the sound speed profile' )
          END IF
 
       END DO

       CLOSE( BTYFile )

       Bot( : )%x( 1 ) = 1000.0 * Bot( : )%x( 1 )   ! Convert ranges in km to m

       ! extend the bathymetry to +/- infinity in a piecewise constant fashion

       Bot( 1 )%x( 1 )       = -sqrt( huge( Bot( 1 )%x( 1 ) ) ) / 1.0d5
       Bot( 1 )%x( 2 )       = Bot( 2 )%x( 2 )
       Bot( NbtyPts )%x( 1 ) = +sqrt( huge( Bot( 1 )%x( 1 ) ) ) / 1.0d5
       Bot( NbtyPts )%x( 2 ) = Bot( NbtyPts - 1 )%x( 2 )

       ! compute tangent and outward-pointing normal to each bottom segment
       !tBot( 1, : ) = xBot( 1, 2:NbtyPts ) - xBot( 1, 1:NbtyPts - 1 )
       !tBot( 2, : ) = xBot( 2, 2:NbtyPts ) - xBot( 2, 1:NbtyPts - 1 )
       ! above caused compiler problems

       DO ii = 1, NbtyPts - 1
          Bot( ii )%t   = Bot( ii+1 )%x - Bot( ii )%x
          Bot( ii )%Len = SQRT( Bot( ii )%t( 1 ) ** 2 + Bot( ii )%t( 2 ) ** 2 )
          Bot( ii )%t   =  Bot( ii )%t / Bot( ii )%Len

          Bot( ii )%n( 1 ) = -Bot( ii )%t( 2 )
          Bot( ii )%n( 2 ) = +Bot( ii )%t( 1 )
       END DO

       IF ( btyType == 'C' ) THEN ! curvilinear option: compute tangent and normal at node by averaging normals on adjacent segments

          DO ii = 2, NbtyPts - 1
             sss = Bot( ii - 1 )%Len / ( Bot( ii - 1 )%Len + Bot( ii )%Len )
             sss = 0.5
             Bot( ii )%Nodet = ( 1.0 - sss ) * Bot( ii - 1 )%t + sss * Bot( ii )%t
          END DO

          Bot( 1       )%Nodet = [ 1.0, 0.0 ]   ! tangent left-end  node
          Bot( NbtyPts )%Nodet = [ 1.0, 0.0 ]   ! tangent right-end node
          Bot( :       )%Noden( 1 ) = -Bot( : )%Nodet( 2 )
          Bot( :       )%Noden( 2 ) = +Bot( : )%Nodet( 1 )

          ! compute curvature in each segment
          phi = atan2( Bot( : )%Nodet( 2 ), Bot( : )%Nodet( 1 ) )   ! this is the angle at each node
          DO ii = 1, NbtyPts - 2
             Bot( ii )%kappa = ( phi( ii + 1 ) - phi( ii ) ) / Bot( ii )%Len ! this is curvature = dphi/ds
          END DO

       ELSE
          Bot%kappa = 0
       END IF

    ELSE   ! no bathymetry given, use SSP depth for flat bottom
       btyType = 'L'
       NbtyPts = 2
       ALLOCATE( Bot( 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP', 'Insufficient memory'  )
       Bot( 1 )%x = [ -rBox, DepthB ]
       Bot( 2 )%x = [  rBox, DepthB ]

       Bot( 1 )%t = [ 1.0, 0.0 ]   ! tangent to bottom
       Bot( 1 )%n = [ 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 2 )%t = [ 1.0, 0.0 ]   ! tangent to bottom
       Bot( 2 )%n = [ 0.0, 1.0 ]   ! outward-pointing normal

       Bot%Len = sqrt( Bot( 1 )%t( 1 )** 2 + Bot( 1 )%t( 2 )** 2 )
       Bot( 1 )%kappa = 0

    ENDIF

  END SUBROUTINE ReadBTY

END MODULE bdrymod
