SUBROUTINE EVAL3D( k, PhiR, PhiS, M, IElementSource, xs, ys, theta, Ntheta, RminM, RmaxM, NR, MACT, P )

  ! Computes 3-D pressure field using adiabatic mode theory           
  ! Normalized to pressure of point source at 1 meter                 
  ! Note RminM must be zero

  USE ElementMod                                          
  IMPLICIT NONE
  INTEGER, PARAMETER :: KBARFile = 55, ZBARFile = 56, maxM = 200
  REAL,    PARAMETER :: PI = 3.141592, DegRad = PI / 180.0
  COMPLEX, PARAMETER :: i  = ( 0.0, 1.0 )
  INTEGER            :: M( * ), Ielt, IElementSource, IR, Itheta, L, MACT, MProp, NEWELT, NR, Ntheta, Outside
  REAL               :: theta( Ntheta ), alpha, delta_r, Rin, Rout, RM, RminM, RmaxM, xs, ys, tsx, tsy
  COMPLEX            :: PhiR( maxM, * ), PhiIn( maxM ), PhiOut( maxM ), PhiS( maxM, * ), const( maxM ), &
                        P( Ntheta, * ), PhiInt
  COMPLEX            :: k(   maxM, * ),  kIn( maxM ),  kOut( maxM ), T, kInt
  COMPLEX, ALLOCATABLE :: sumk( : )

  ! Open file for eigenfunctions                                  
  OPEN ( FILE = 'ZBARFile', UNIT = ZBARFile, STATUS = 'UNKNOWN' ) 

  !  *** Loop over angle ***                                           

  delta_r = ( RmaxM - RminM ) / ( NR - 1 ) 

  DO Itheta = 1, Ntheta 
     tsx = COS( DegRad * theta( Itheta ) )
     tsy = SIN( DegRad * theta( Itheta ) ) 

     ! Get modal values                                           
     Ielt = IElementSource 
     WRITE( *, * ) 'Tracing bearing ', Itheta, ATAN2( tsy, tsx ) / DegRad 

     CALL FIRST( Ielt, Outside, RIn, ROut, xS, yS, tsx, tsy,        &
          MProp, M, maxM, k, PhiR, PhiS, const, kIn, PhiIn, kOut, PhiOut )                          

     MProp = MIN( MACT, MProp ) 

     ! write modes at first range                                 
     IF ( Itheta == 1 ) THEN 
        WRITE( ZBARFile, * ) MProp 
        WRITE( ZBARFile, * ) const( 1 : MProp ) 
     END IF

     const( 1 : MProp ) = i * SQRT( 2.0 * PI ) * EXP( i * PI / 4.0 ) * const( 1 : MProp )

     IF ( Itheta == 1 ) ALLOCATE( sumk( MProp ) )
     sumk = 0.0

     ! *** March forward in range ***                                 

     DO IR = 1, NR 
        RM = RminM + ( IR - 1 ) * delta_r 
        IF ( RM == 0.0 ) RM = MIN( 1.0, delta_r ) 

        ! Crossing into new element?                                  
        DO WHILE ( RM > ROut ) 

           ! Copy outside info to inside                          
           NEWELT = AdjElt( Outside, Ielt ) 
           RIn    = ROut 

           kIn(   1 : MProp ) =   kOut( 1 : MProp ) 
           PhiIn( 1 : MProp ) = PhiOut( 1 : MProp ) 

           ! Get new outside info                                 
           CALL OUT( Ielt, NEWELT, Outside, ROut, xS, yS, tsx, tsy, MProp, M, maxM, k, PhiR, kOut, PhiOut )             
           Ielt = NEWELT 

        END DO

        ! *** Compute modal contribution at this range ***            

        T = 0.0 

        IF ( RIn /= ROut ) THEN 
           alpha = ( RM - RIn ) / ( ROut - RIn ) 
           alpha = MIN( MAX( alpha,  0.0 ), 1.0 ) 
        ELSE 
           alpha = 0.0 
        ENDIF

        IF ( IR == NR ) WRITE( ZBARFile, * ) MProp 
        DO L = 1, MProp 
           kInt      =   kIn( L ) + alpha   * (   kOut( L ) -   kIn( L ) ) 
           PhiInt    = PhiIn( L ) + alpha * ( PhiOut( L ) - PhiIn( L ) ) 
           sumk( L ) =  sumk( L ) + delta_r * kInt 
           T         = T + PhiInt * const( L ) * EXP( -i * sumk( L ) ) / SQRT( kInt )    
           IF ( IR == NR ) WRITE( ZBARFile, * ) PhiInt  ! write mode at last range     
        END DO

        P( Itheta, IR ) = CMPLX( T / SQRT( RM ) )

     END DO   ! Next range

     ! Write average wavenumber to file                           
     OPEN ( FILE = 'KBARFile', UNIT = KBARFile, STATUS = 'UNKNOWN' ) 
     WRITE( KBARFile, * ) MProp, RM 
     WRITE( KBARFile, * ) sumk( 1:MProp ) / RM

  END DO ! Next bearing

END SUBROUTINE EVAL3D
!**********************************************************************C
SUBROUTINE FIRST( Ielt, Outside, RIn, ROut, xS, yS, tsx, tsy,     &
     MProp, M, maxM, k, PhiR, PhiS, const, kIn, PhiIn, kOut, PhiOut )

  ! Given an element number for the source                            
  ! Computes                                                          
  !    Outside    the side through which path exits                   
  !    ROut       the range at which path exits                       
  !    PhiIn      Mode values at entrance
  !    PhiOut     Mode values at exit
  !    const      Interpolated mode excitation coeffs                            

  USE ElementMod
  IMPLICIT NONE
  INTEGER                :: Inside, Outside, M( * ), ICor( 3, 2 ), IBad, IGood1, IGood2, &
                            ICor1, ICor2, ICor3, ICor4, Ielt, &
                            ISide, IS, Iset1, Iset2, J, L, maxM, MProp, node1, node2 
  REAL                   :: R,  RIn, ROut, RV( 3 ), SV( 3 ), RVC( 3 ), SIn, SOut, xS, yS, x1S, y1S, xCenter, yCenter, &
                            Delta, Tx, Ty, TxC, TyC, alpha, tsx, tsy
  COMPLEX, INTENT( IN  ) :: PhiS( maxM, * ), PhiR( maxM, * )
  COMPLEX, INTENT( IN  ) :: k(   maxM, * ),  kIn( maxM ),  kOut( maxM )
  COMPLEX, INTENT( OUT ) :: PhiIn( * ), PhiOut( * ), Const( * )          
  DATA ( (ICor( IS, J), IS = 1, 3), J = 1, 2) /1, 2, 3, 2, 3, 1/ 

  ! ICor maps a side (1, 2 OR 3) and a local node (1 OR 2) to a       
  !      corner (1, 2, or 3) of the triangle                        

  ! WRITE( *, * ) 'Ielt = ', Ielt 
  MProp = HUGE( MProp )

  ! Coordinates of the centroid of the source element                 

  xCenter = SUM( x( node( 1 : 3, Ielt ) ) ) / 3.0
  yCenter = SUM( y( node( 1 : 3, Ielt ) ) ) / 3.0 

  DO ISide = 1, 3 

     node1 = node( ICor( ISide, 1 ), Ielt )
     node2 = node( ICor( ISide, 2 ), Ielt )
     Iset1 = Iset( node1 )
     Iset2 = Iset( node2 ) 
     MProp = MIN( MProp, M( Iset1 ), M( Iset2 ) ) 

     x1S = x( node1 ) - xS
     y1S = y( node1 ) - yS
     TxC = x( node1 ) - xCenter
     TyC = y( node1 ) - yCenter                                 
     Tx  = x( node2 ) - x( node1 )
     Ty  = y( node2 ) - y( node1 ) 

     Delta = tsx * Ty - tsy * Tx 

     ! *** Radial parallel to side? ***                               

     IF ( Delta == 0.0 ) THEN 
        SV(  ISide )    = HUGE( SV( ISide ) )
     ELSE 
        RVC( ISide ) = ( TxC * Ty  - TyC * Tx  ) / Delta 
        RV(  ISide ) = ( x1S * Ty  - y1S * Tx  ) / Delta 
        SV(  ISide ) = ( x1S * tsy - y1S * tsx ) / Delta 
     ENDIF
  END DO

  ! Identify two good sides and one bad side based on the intercept point

  IBad = 1 
  IF ( ABS( SV( 2 ) - 0.5 ) > ABS( SV( IBad ) - 0.5 ) ) IBad = 2                                                       
  IF ( ABS( SV( 3 ) - 0.5 ) > ABS( SV( IBad ) - 0.5 ) ) IBad = 3                                                     

  IGood1 = 1
  IGood2 = 2
  IF ( IBad == 1 ) THEN 
     IGood1 = 3 
  ELSE IF ( IBad == 2 ) THEN 
     IGood2 = 3 
  ENDIF

  ! The side with the lesser RVC is the inside                    

  IF ( RVC( IGood1 ) < RVC( IGood2 ) ) THEN 
     Inside  = IGood1
     Outside = IGood2 
  ELSE 
     Inside  = IGood2
     Outside = IGood1 
  ENDIF

  SIn      = SV( Inside ) 
  SIn      = MIN( MAX( SIn,  0.0 ), 1.0 ) 
  RIn      = RV( Inside ) 

  SOut     = SV( Outside ) 
  SOut     = MIN( MAX( SOut, 0.0 ), 1.0 ) 
  ROut     = RV( Outside ) 

  ! Get values of modes at Z = SD and (x, y) = intercept points   

  ICor1 = ICor(  Inside, 1 )
  ICor2 = ICor(  Inside, 2 ) 
  ICor3 = ICor( Outside, 1 )
  ICor4 = ICor( Outside, 2 ) 

  ! Interpolate to get modal values at source                     
  R = 0.0 

  IF ( RIn /= ROut ) THEN 
     alpha = ( R - RIn ) / ( ROut - RIn ) 
     alpha = MIN( MAX( alpha,  0.0 ), 1.0 ) 
  ELSE 
     alpha = 0.0 
  ENDIF

  DO L = 1, MProp 
     PhiIn(  L ) = PhiS( L, ICor1 ) + SIn   * ( PhiS( L, ICor2 ) - PhiS( L, ICor1 ) )   
     PhiOut( L ) = PhiS( L, ICor3 ) + SOut  * ( PhiS( L, ICor4 ) - PhiS( L, ICor3 ) )   
     const(  L ) = PhiIn( L )       + alpha * ( PhiOut( L )      - PhiIn( L ) ) 
  END DO

  ! Obtain values of modes at z = Rd and (x, y)=intercept points      

  CALL INTER1( Ielt, Inside,  SIn,  MProp, M, maxM, k, PhiR, kIn,  PhiIn  )          
  CALL INTER1( Ielt, Outside, SOut, MProp, M, maxM, k, PhiR, kOut, PhiOut )          
 
END SUBROUTINE FIRST
!**********************************************************************C
SUBROUTINE OUT( Ielt, NEWELT, Outside, ROut, xS, yS, &
     &   tsx, tsy, MProp, M, maxM, k, PhiR, kOut, PhiOut)          

  ! Given an element number and a side through which prop path enters 
  ! Computes                                                          
  !    Outside    the side through which path exits                   
  !    ROut       the range at which path exits                       
  !    Mode values                                                    

  USE ElementMod
  IMPLICIT NONE
  INTEGER          :: Outside, M( * ), ICor( 3, 2 ), Ielt, ISide, IS, J, NEWELT, maxM, MProp, node1, node2
  REAL             :: ROut, ROutT, SOut, ST, x1S, y1S, xS, yS, Delta, Tx, Ty, tsx, tsy
  COMPLEX          :: PhiOut( * ), PhiR( maxM, * )
  COMPLEX          :: k(   maxM, * ), kOut( maxM )
  DATA ( (ICor( IS, J), IS = 1, 3), J = 1, 2) /1, 2, 3, 2, 3, 1/ 

  ! ICor maps a side (1, 2 OR 3) and a local node (1 OR 2) to a       
  !        corner (1, 2, or 3) of the triangle                        
  ! WRITE( *, * ) '   Crossing into new element = ', NEWELT           

  ! If no adj elt then exit. Previous values of KOut, PhiOut retained 

  IF ( NEWELT == 0 ) THEN 
     ROut = HUGE( ROut )
     RETURN 
  ENDIF

  ! *** Loop over sides to find outside ***                           

  SOut = HUGE( SOut )

  DO ISide = 1, 3 

     ! Don't try an outside the same as the inside             
     IF ( AdjElt( ISide, NEWELT ) /= Ielt ) THEN 

        node1 = node( ICor( ISide, 1 ), NEWELT ) 
        node2 = node( ICor( ISide, 2 ), NEWELT ) 

        x1S = x( node1 ) - xS
        y1S = y( node1 ) - yS                                               
        Tx  = x( node2 ) - x( node1 )
        Ty  = y( node2 ) - y( node1 ) 

        Delta = tsx * Ty - tsy * Tx 

        ! *** Radial parallel to side? ***                            

        IF ( Delta == 0.0 ) THEN 
           ST    = HUGE( ST )
        ELSE 
           ROutT = ( x1S * Ty  - y1S * Tx  ) / Delta 
           ST    = ( x1S * tsy - y1S * tsx ) / Delta 
        ENDIF
        ! If intercept is in segment, store the side number       
        IF ( ABS( ST - 0.5 ) < ABS( SOut - 0.5) ) THEN 
           Outside = ISide 
           SOut    = ST 
           ROut    = ROutT 
        ENDIF
     ENDIF
  END DO

  ! Obtain values of modes at intercept points                        

  CALL INTER1( NEWELT, Outside, SOut, MProp, M, maxM, k, PhiR, kOut, PhiOut )          
 
END SUBROUTINE OUT
!**********************************************************************C
SUBROUTINE INTER1( Ielt, ISide, S, MProp, M, maxM, k, PhiR, kInt, PhiInt )          

  ! Given                                                             
  !    Ielt         the element                                       
  !    ISide        the side                                          
  !    S            the proportional distance along the side          

  ! Returns                                                           
  !    kInt, PhiInt  interpolated modal values                       
  !    MProp         number of propagating modes                     

  USE ElementMod
  IMPLICIT NONE
  INTEGER, INTENT( IN  ) :: Ielt, ISide
  INTEGER, INTENT( OUT ) :: MProp
  INTEGER                :: M( * ), ICor( 3, 2 ), Iset1, Iset2, I, IS, J, maxM, node1, node2 
  REAL,     INTENT( IN ) :: S
  REAL                   :: St
  COMPLEX                :: PhiR( maxM, * )
  COMPLEX, INTENT( OUT ) :: PhiInt( * )
  COMPLEX                :: k(   maxM, * )
  COMPLEX, INTENT( OUT ) :: kInt( maxM )

  DATA ( (ICor( IS, J), IS = 1, 3), J = 1, 2) /1, 2, 3, 2, 3, 1/ 

  ! ICor maps a side   (1, 2, or 3) and a local node (1 or 2) to
  !           a corner (1, 2, or 3) of the triangle                        

  node1 = node( ICor( ISide, 1 ), Ielt )
  node2 = node( ICor( ISide, 2 ), Ielt )                                         
  Iset1 = Iset( node1 )
  Iset2 = Iset( node2 ) 
  MProp = MIN( MProp, M( Iset1 ), M( Iset2 ) ) 

  ! Extrapolation is blocked by making sure s is in [0, 1]                                   

  St = MIN( MAX( S, 0.0 ), 1.0 ) 

  DO I = 1, MProp 
     kInt(   I ) =    k( I, Iset1 ) + St   * (    k( I, Iset2 ) -    k( I, Iset1 ) )  
     PhiInt( I ) = PhiR( I, Iset1 ) + St * ( PhiR( I, Iset2 ) - PhiR( I, Iset1 ) )  
  END DO

END SUBROUTINE INTER1
