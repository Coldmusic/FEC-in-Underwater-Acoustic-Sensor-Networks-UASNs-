PROGRAM FIELD3D 

  ! Generate file of replica vectors
  ! Multiple receivers depths are handled inefficiently.

  USE SdRdRMod
  USE ElementMod
  IMPLICIT NONE
  INTEGER, PARAMETER   :: SHDFile = 25, MaxM = 200, MaxSet = 9000
  !  if you change MaxM, change it in eval3d.f, evalpdq.f, evalgbt.f too
  INTEGER              :: M( MaxSet ), Is, Ir, Itheta, Msource, NSets, IREC, IElementSource, ICorner, INode, Mlimit, &
                          IdentifySourceElement
  REAL                 :: Rmin, Rmax, xs, ys, freq, atten                             
  COMPLEX              :: k( MaxM, MaxSet ), PhiR( MaxM, MaxSet ), PhiST( MaxM, 3 )                   
  CHARACTER   (LEN=50) :: Option
  CHARACTER   (LEN=80) :: Title, TitleEnv, FileRoot, SHDFileName
  CHARACTER   (LEN=10) :: PlotType
  COMPLEX, ALLOCATABLE :: P( :, : ), PhiS( :, :, : )

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) '__________________________________________________'
  WRITE( PRTFile, * ) 'Running FIELD3D'
  WRITE( PRTFile, * ) 'Sums modes, producing pressure'
  WRITE( PRTFile, * )

  CALL READIN( FileRoot, Title, Option, Mlimit, XS, ys, Rmin, Rmax )  ! Read in all the user input
  ALLOCATE( PhiS( MaxM, Nsd, 3 ) )
  CALL GETADJ                                               ! Build AdjElt table
  IElementSource = IdentifySourceElement( XS, ys )                         ! Identify the source element

  ! Read modes at the source depths
  Msource = MaxM 
  DO ICorner = 1, 3 
     INode = Node( ICorner, IElementSource ) 
     CALL GetModes( ModeFileName( INode ), 1, MaxM, sd, Nsd, 'N', &
                  k( 1, ICorner ), PhiS( 1, 1, ICorner ), M( ICorner ), freq, TitleEnv )
     Msource = MIN( Msource, M( ICorner ) ) 
  END DO

  ! Write header
  SHDFileName  = TRIM( FileRoot ) // '.shd'
  atten    = 0.0
  CALL WriteHeader( SHDFileName, Title, theta, Ntheta, sd, Nsd, rd, Nrd, R, NR, freq, atten, PlotType, xs, ys )
  ALLOCATE( P( Ntheta, NR ) )

  ! MAIN Loop: over receiver depth

  DO Ir = 1, Nrd
     CALL SLURP( MaxSet, MaxM, rd( IR ), NSets, M, k, PhiR )   ! Get modes at the receiver depth   

     ! Loop over source depth
     DO Is = 1, Nsd 
        PhiST( 1:Msource, : ) = PhiS( 1:Msource, IS, : ) ! Get modes at the source depth

        ! Call the appropriate routine to evaluate the field
        SELECT CASE ( Option(1:3) )
        CASE ( 'STD' )  ! STANDARD (ignores hor. refraction)
           CALL EVAL3D(  k, PhiR, PhiST, M,       IElementSource, xs, ys, theta, Ntheta, Rmin, Rmax, NR, Mlimit, P )       
        CASE ( 'PDQ' )  ! Quick (ignores and ignores)
           CALL EVALPDQ( k, PhiR, PhiST, M,       IElementSource, xs, ys, theta, Ntheta, Rmin, Rmax, NR, Mlimit, P )       
        CASE ( 'GBT' )  ! GAUSSIAN BEAM (include hor. refraction)
           CALL EVALGB(  k, PhiR, PhiST, M, MaxM, IElementSource, xs, ys, theta, Ntheta, Rmin, Rmax, NR, Mlimit, Option, P,&
                         freq )       
        CASE DEFAULT
           CALL ERROUT( PRTFile, 'F', 'FIELD3D', 'Unknown option' ) 
        END SELECT

        ! Write out the field

        DO Itheta = 1, Ntheta
           IREC = 7 + ( Itheta - 1 ) * Nsd * Nrd + ( Is - 1 ) * Nrd + Ir
           WRITE( SHDFile, REC = IRec ) P( Itheta, 1 : Nr )
        END DO

     END DO   ! Next source   depth 

  END DO   ! Next receiver depth

END PROGRAM FIELD3D

!**********************************************************************C

SUBROUTINE READIN( FileRoot, Title, Option, Mlimit, xs, ys, Rmin, Rmax )

  ! Reads in all the user input

  USE SdRdRMod
  USE ElementMod

  IMPLICIT NONE
  INTEGER, PARAMETER :: FLPFile = 5, MaxM = 200
  INTEGER            :: I, IAllocStat, Mlimit, IElt, itheta, iostat
  REAL               :: xs, ys, Rmin, Rmax
  CHARACTER (LEN=50) :: Option
  CHARACTER (LEN=80) :: Title, FileRoot

! Open the field parameter file
  OPEN( UNIT = FLPFile, FILE = TRIM( FileRoot ) // '.flp', STATUS = 'OLD', IOSTAT = iostat )
  IF ( IOSTAT /= 0 ) THEN   ! successful open?
     WRITE( PRTFile, * ) 'FLPFile = ', TRIM( FileRoot ) // '.flp'
     CALL ERROUT( PrtFile, 'F', 'FIELD3D - READIN', 'Unable to open the field parameter file' )
  END IF

  READ(  FLPFile, * ) Title 
  WRITE( PRTFile, * ) Title 

  READ(  FLPFile, * ) Option 
  WRITE( PRTFile, * ) 'Option = ', Option

  READ(  FLPFile, * ) Mlimit 
  WRITE( PRTFile, * ) 'Number of modes = ', Mlimit 
  Mlimit = MIN( Mlimit, MaxM ) 

  ! Read source/receiver information
  READ(  FLPFile, * ) xs, ys 
  WRITE( PRTFile, * ) 'Coords. of source = ', xs, ys 
  xs = 1000.0 * xs   ! convert km to m
  ys = 1000.0 * ys 

  CALL SDRD(   FLPFile, PRTFile, 0.0, 1.0E6 )  ! Read source/rcvr depths
  CALL RANGES( FLPFile, PRTFile )              ! Read receiver ranges

  Rmin = R(  1 )
  Rmax = R( NR )

  ! Read angles for radials
  READ(  FLPFile, * ) Ntheta
  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Number of bearings = ', Ntheta
  ALLOCATE( theta( MAX( 3, Ntheta ) ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) THEN
     WRITE( PRTFile, * ) 'Ntheta = ', Ntheta
     CALL ERROUT( PRTFile, 'F', 'SDRD', 'Too many bearings'  )
  END IF

  theta( 3 ) = -999.9 
  READ( FLPFile, * ) theta( 1 : Ntheta )
  CALL SUBTAB( theta, Ntheta ) 
  IF ( Ntheta >= 1 ) WRITE( PRTFile, "( 5G14.6 )" ) ( theta( itheta ), itheta = 1, MIN( Ntheta, Number_to_Echo ) )
  IF ( Ntheta >  Number_to_Echo ) WRITE( PRTFile, * ) ' ... ', theta( Ntheta )

  ! Read nodal coordinates
  READ(  FLPFile, * ) NNodes 
  WRITE( PRTFile, * ) 'NNodes = ', NNodes
  ALLOCATE( X( NNodes ), Y( NNodes ), ModeFileName( NNodes ), Iset( NNodes ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) THEN 
     CALL ERROUT( PRTFile, 'F', 'FIELD3D-READIN', 'Too many nodes' ) 
  ENDIF

  DO I = 1, NNodes 
     READ( FLPFile, * ) X( I ), Y( I ), ModeFileName( I ) 
     X( I ) = 1000.0 * X( I )
     Y( I ) = 1000.0 * Y( I )
  END DO

  ! Read in element definitions
  READ(  FLPFile, * ) NElts 
  WRITE( PRTFile, * ) 'NElts  = ', NElts 
  ALLOCATE( AdjElt( 3, NElts ), Node( 3, NElts ) )

  DO IElt = 1, NElts 
     READ( FLPFile, * ) Node( 1 : 3, IELt ) 
  END DO

  IF ( Option(4:4) == 'T' ) THEN 
     WRITE( PRTFile, * ) 'Performing a Tesselation check' 
     CALL TESCHK
     WRITE( PRTFile, * ) 'Passed the Tesselation check' 
  ENDIF

END SUBROUTINE READIN
!**********************************************************************C
SUBROUTINE GETADJ

  ! Constructs a table AdjElt(IElt, iside) which gives the
  ! element number that shares iside with element number IElt

  USE ElementMod
  IMPLICIT NONE
  INTEGER ICorner( 3, 2 ), IElt, IEltT, iside, isideT, Node1, Node1T, Node2, Node2T, J
  DATA ((ICorner(iside, J), iside=1, 3), J=1, 2) /1, 2, 3, 2, 3, 1/ 

  ! ICorner maps a side (1, 2 or 3) and a local node (1 or 2) to a
  ! corner (1, 2, OR 3) of the triangle

  AdjElt = 0

  DO IElt = 1, NElts  ! Loop over each triangle
     
     SIDE: DO iside = 1, 3 ! Loop over triangle sides

        IF ( AdjElt( iside, IElt ) == 0 ) THEN 

           Node1 = Node( ICorner( iside, 1 ), IElt ) 
           Node2 = Node( ICorner( iside, 2 ), IElt ) 

           ! Search through other elements to find common side
           DO IEltT = 1, NElts 
              IF ( IEltT /= IElt ) THEN 
                 DO isideT = 1, 3 
                    Node1T = Node( ICorner( isideT, 1 ), IEltT ) 
                    Node2T = Node( ICorner( isideT, 2 ), IEltT ) 

                    ! Do IElt and IEltT share this side?
                    IF ( ( Node1 == Node1T .AND. Node2 == Node2T ) .OR. &
                         & ( Node1 == Node2T .AND. Node2 == Node1T ) ) THEN           
                       AdjElt( iside,  IElt  ) = IEltT 
                       AdjElt( isideT, IEltT ) = IElt 
                       CYCLE SIDE
                    ENDIF

                 END DO   ! next side
              ENDIF
           END DO   ! next element
        ENDIF

     END DO SIDE   ! next side

  END DO   ! next element

END SUBROUTINE GETADJ

!**********************************************************************C

FUNCTION IdentifySourceElement( xs, ys ) 

  !     Identifies the element containing the source at (xs, ys)          

  !     We define a function ENCL( xs, ys ) which is 1 when (xs, ys)     
  !     is inside a given triangle and decreases from 1 the further the   
  !     source moves away from the triangle.

  !     The element with the highest value of ENCL is identified as the source element.
  !     If several elements enclose, the highest numbered element is given posession.

  USE ElementMod
  IMPLICIT NONE
  INTEGER :: IELT, Node1, Node2, Node3, IdentifySourceElement
  REAL    :: xs, ys, X1, X2, X3, Y1, Y2, Y3, A1, A2, A3, Delta, ENCL, ENCLMax

  ENCLMax = 0.0 

  DO IElt = 1, NElts 
     Node1 = Node( 1, IElt )
     Node2 = Node( 2, IElt ) 
     Node3 = Node( 3, IElt ) 

     X1 = X( Node1 )   ;   Y1 = Y( Node1 ) 
     X2 = X( Node2 )   ;   Y2 = Y( Node2 ) 
     X3 = X( Node3 )   ;   Y3 = Y( Node3 ) 

     ! Compute areas of triangles
     Delta = ( X2*Y3 - Y2*X3 ) - ( X1*Y3 - Y1*X3 ) + ( X1*Y2 - Y1*X2 )
     A1    = ( X2*Y3 - Y2*X3 ) - ( xs*Y3 - ys*X3 ) + ( xs*Y2 - ys*X2 )
     A2    = ( xs*Y3 - ys*X3 ) - ( X1*Y3 - Y1*X3 ) + ( X1*ys - Y1*xs )
     A3    = ( X2*ys - Y2*xs ) - ( X1*ys - Y1*xs ) + ( X1*Y2 - Y1*X2 )

     ENCL = ABS( DEltA ) / ( ABS( A1 ) + ABS( A2 ) + ABS( A3 ) ) 

     IF ( ENCL > ENCLMax ) THEN 
        IdentifySourceElement  = IElt 
        ENCLMax = ENCL 
     ENDIF

  END DO

END FUNCTION IdentifySourceElement
!**********************************************************************C
SUBROUTINE TESCHK 

  !     Checks to see that triangulation is a tesselation                 
  !     that is, that there are no overlapping triangles                  
  !     (holes or triangles inside triangles are still possible)          

  USE ElementMod
  IMPLICIT NONE
  INTEGER  :: ICorner( 3, 2 ), IElt1, IElt2, iside1, iside2, iside, J
  DATA ((ICorner(iside, J), iside=1, 3), J=1, 2) /1, 2, 3, 2, 3, 1/
  REAL     :: X1, X2, X3, X4, Y1, Y2, Y3, Y4, UX, UY, VX, VY, WX, WY, DEL, S1, S2

  ! ICorner maps a side (1, 2 or 3) and a local node (1 or 2) to a
  ! corner (1, 2, OR 3) of the triangle

  DO IElt1 = 1, NElts-1 
     DO iside1 = 1, 3 
        X1 = X( Node( ICorner( iside1, 1 ), IElt1 ) ) 
        Y1 = Y( Node( ICorner( iside1, 1 ), IElt1 ) ) 
        X2 = X( Node( ICorner( iside1, 2 ), IElt1 ) ) 
        Y2 = Y( Node( ICorner( iside1, 2 ), IElt1 ) ) 
        UX = X2 - X1   ;   UY = Y2 - Y1 

        DO IElt2 = IElt1 + 1, NElts 
           DO iside2 = 1, 3 
              X3 = X( Node( ICorner( iside2, 1 ), IElt2 ) ) 
              Y3 = Y( Node( ICorner( iside2, 1 ), IElt2 ) ) 
              X4 = X( Node( ICorner( iside2, 2 ), IElt2 ) ) 
              Y4 = Y( Node( ICorner( iside2, 2 ), IElt2 ) ) 
              VX = X4 - X3
              VY = Y4 - Y3 
              WX = X4 - X2
              WY = Y4 - Y2 
              DEL = UX * VY - UY * VX 

              IF ( ABS( DEL ) > MAX( ABS( UX*UY ), ABS( VX*VY ) ) / 10000.0 ) THEN 
                 S1 = ( UX * WY - UY * WX ) / DEL 
                 S2 = ( VX * WY - VY * WX ) / DEL 
                 IF ( 0.001 < S1 .AND. S1 < 0.999 .AND.       &
                  &   0.001 < S2 .AND. S2 < 0.999 ) THEN      
                    WRITE( PRTFile, * ) 'Sides cross for elts = ', IElt1, IElt2      
                    WRITE( PRTFile, * ) 'SIDE 1:' 
                    WRITE( PRTFile, * ) '   (', X1, Y1, ')' 
                    WRITE( PRTFile, * ) '   (', X2, Y2, ')' 
                    WRITE( PRTFile, * ) '   S = ', S1 

                    WRITE( PRTFile, * ) 'SIDE 2:' 
                    WRITE( PRTFile, * ) '   (', X3, Y3, ')' 
                    WRITE( PRTFile, * ) '   (', X4, Y4, ')' 
                    WRITE( PRTFile, * ) '   S = ', S2 
                    STOP 
                 ENDIF
              ENDIF

           END DO
        END DO

     END DO
  END DO

END SUBROUTINE TESCHK

!**********************************************************************C

SUBROUTINE SLURP( MaxSet, MaxM, rd, NSets, M, k, PhiR )      

  ! Reads in the values of the modes at the given receiver depth
  ! for every node in the triangulation

  USE ElementMod
  IMPLICIT NONE
  INTEGER, INTENT( IN  ) :: MaxSet, MaxM
  INTEGER, INTENT( OUT ) :: M( * )
  INTEGER                :: INode, JNode, NSets
  REAL                   :: freq, rd
  COMPLEX, INTENT( OUT ) :: k( MaxM, * ), PhiR( MaxM, * ) 
  CHARACTER     (LEN=80) :: TitleEnv

  NSets = 0 

  NodeLoop: DO INode = 1, NNodes 

     ! Check if the modes have already been read
     IF ( INode >= 2 ) THEN 
        DO JNode = 1, INode-1 
           IF ( ModeFileName( INode ) == ModeFileName( JNode ) ) THEN 
              Iset( INode ) = Iset( JNode )  ! Copy previously read modes
              CYCLE NodeLoop  ! process next node 
           ENDIF
        END DO
     ENDIF

     NSets = NSets + 1 
     IF ( NSets > MaxSet ) THEN 
        WRITE( PRTFile, * ) 'MaxSet = ', MaxSet 
        CALL ERROUT( PRTFile, 'F', 'FIELD3D-SLURP', 'Too many mode sets' )               
     ENDIF
     Iset( INode ) = NSets 

     ! Check for 'DUMMY' elts (acoustic absorbers)
     IF ( ModeFileName( INode )(1:5) == 'DUMMY' ) THEN 
        M( NSets ) = 0 
     ELSE  ! Necessary to read in modes
        !WRITE( *, * ) 'Reading INode', INode 
        CALL GetModes( ModeFileName( INode ), 1, MaxM, rd, 1, 'N', k( 1, NSets ), PhiR( 1, NSets ), M( NSets ), freq, TitleEnv )
     ENDIF

  ENDDO NodeLoop

END SUBROUTINE SLURP
