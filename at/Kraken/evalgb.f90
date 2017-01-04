SUBROUTINE EVALGB( k, phi, phiR, M, maxM, IelementSource, xS, yS, theta, Ntheta, Rmin, Rmax, NR, MActive, Option, P, freq )

  ! Computes 3-D pressure field using adiabatic mode theory.
  ! Uses gaussian beam tracing for horizontal refraction.           
  ! Note beam curvature change at interfaces still needed.          
  ! Phase may be off by e(i pi / 4 ) but TL is correct.
  ! Note rays here should be complex; they are projected here onto the real plane
  ! This happens wherever we use DBLE( CA ) instead of the complex CA

  USE ElementMod
  IMPLICIT NONE
  INTEGER, PARAMETER :: FLPFIL = 5, RAYFIL = 21
  REAL,    PARAMETER :: pi = 3.141592, DegRad = pi / 180.0, c0 = 1500.0
  COMPLEX, PARAMETER :: i = ( 0.0, 1.0 )
  LOGICAL            :: EXITED
  INTEGER            :: M( * ), Ialpha, IElement, IelementSource, Iset1, Iset2, Iset3, Iside, Istep, KMAHA, KMAHB, &
                        MActive, maxM, mode, Mprop, Nsteps, NR, Nalpha, Ntheta
  REAL               :: theta( * ),  freq, Rmin, Rmax, xs, ys
  REAL, ALLOCATABLE  :: RadVec( :, : )
  REAL (KIND=8), ALLOCATABLE :: xV( : ), yV( : )
  REAL    (KIND=8)   :: alpha, Dalpha, A1, A2, A3, D12, D13, D23, DB1, DB2, DB3, delta, &
                        EpsilonMult, qAT, qBT, tsx, tsy, xa, xb, ya, yb, x1, y1, x2, y2, x3, y3, &
                        xia, xib, etaa, etab, alpha1, alpha2, step, stepG, Hwidth
  COMPLEX            :: k( maxM, * ), phi( maxM, * ), phiR( maxM, 3 ), P( Ntheta, * )
  COMPLEX (KIND=8)   :: const, phiA, phiB, pA, pB, qA, qB, tauA, tauB, CA, CB, eps, EpsOpt, Cx, Cy
  CHARACTER          :: Option*( * )

  IF ( Option(6:6) == 'R' ) OPEN( FILE = 'RAYFIL', UNIT = RAYFIL ) ! optional ray file

  ! following should only be done on first call
  READ( FLPFIL, * ) alpha1, alpha2, Nalpha 
  READ( FLPFIL, * ) stepG, Nsteps 
  READ( FLPFIL, * ) EpsilonMult
  
  ! automatic calculation of number of beams
  IF ( Nalpha == 0 ) THEN
     Nalpha = MAX( INT( 0.3 * Rmax * Freq / c0 * ( alpha2 - alpha1 ) / 360 ), 300 )
  END IF

  WRITE( PRTFile, * ) 'alpha1, alpha2    = ', alpha1, alpha2
  WRITE( PRTFile, * ) 'Number of beams   = ', Nalpha
  WRITE( PRTFile, * ) 'stepsize          = ', stepG
  WRITE( PRTFile, * ) 'Number of steps   = ', Nsteps
  WRITE( PRTFile, * ) 'Epsilon multplier = ', EpsilonMult

  alpha1 = DegRad * alpha1 
  alpha2 = DegRad * alpha2 
  Dalpha = ( alpha2 - alpha1 ) / ( Nalpha - 1 )

  ! Radial vector for each bearing line
  ALLOCATE( RadVec( 2, Ntheta ), xV( Nsteps ), yV( Nsteps ) )
  RadVec( 1, : )        = COS( DegRad * theta( 1:Ntheta ) )
  RadVec( 2, : )        = SIN( DegRad * theta( 1:Ntheta ) )
  P( 1:Ntheta, 1:NR )   = 0.0

  ! Loop over beam takeoff angle
  ANGLE: DO Ialpha = 1, Nalpha                                  
     alpha = alpha1 + ( Ialpha - 1 ) * Dalpha 
     tsx   = COS( alpha ) 
     tsy   = SIN( alpha ) 
     WRITE( *, * ) 'Ialpha, tsx, tsy', Ialpha, tsx, tsy 

     ! Loop over modes
     DO mode = 1, MActive 

        ! Get mode excitation coef and initialize
        IElement = IelementSource 
        CALL NEWELT( IElement, k, mode, M, maxM, Iset1, Iset2, Iset3, x1, y1, x2, y2, x3, y3, D12, D13, D23, delta, Cx, Cy, Mprop )
        IF ( mode > Mprop ) CYCLE ANGLE

        ! Evaluate modes at source depth                          
        DB1 = xS * y1 - yS * x1 
        DB2 = xS * y2 - yS * x2 
        DB3 = xS * y3 - yS * x3 

        ! Compute areas of triangles                              
        A1   = (  D23 - DB3 + DB2 ) / delta 
        A2   = (  DB3 - D13 - DB1 ) / delta 
        A3   = ( -DB2 + DB1 + D12 ) / delta 
        CA   = A1 / k( mode, Iset1 ) + A2 / k( mode, Iset2 ) + A3 / k( mode, Iset3 )
        phiA = A1 * phiR( mode, 1  ) + A2 * phiR( mode, 2  ) + A3 * phiR( mode, 3  )

        IF ( Option(5:5) == 'F' ) THEN ! Space filling in far field
           Hwidth = 2.0 / ( ( 1.0 / DBLE( CA ) ) * Dalpha ) 
           EpsOpt = 0.5 * Hwidth ** 2 
        ENDIF

        IF ( Option(5:5) == 'M' ) THEN ! Minimum width at Rmax
           Hwidth = SQRT( 2.0 * DBLE( CA ) * Rmax ) 
           EpsOpt = 0.5 * Hwidth ** 2 
        ENDIF

        eps   = EpsilonMult * i * EpsOpt 
        const = phiA * SQRT( eps / CA ) * Dalpha 
        xA    = xS 
        yA    = yS 
        xiA   = tsx / DBLE( CA )
        etaA  = tsy / DBLE( CA ) 
        pA    = 1.0 
        qA    = eps 
        tauA  = 0.0 
        KMAHA = 1 

        ! Evaluate modes at rcvr depth                            
        phiA = A1 * phi( mode, Iset1 ) + A2 * phi( mode, Iset2 ) + A3 * phi( mode, Iset3 )        

        ! March forward in range
        EXITED = .FALSE. 
        step = stepG 

        DO Istep = 1, Nsteps 
           xB   = xA   + step * DBLE( CA ) * xiA 
           yB   = yA   + step * DBLE( CA ) * etaA 
           xiB  = xiA  - step * DBLE( Cx / ( CA * CA ) )
           etaB = etaA - step * DBLE( Cy / ( CA * CA ) ) 
           pB   = pA 
           qB   = qA   + step * DBLE( CA ) * pA
           tauB = tauA + step / DBLE( CA ) 

           step = stepG ! Update step size

           IF ( Option(6:6) == 'R' ) THEN 
              xV( Istep ) = xA 
              yV( Istep ) = yA 
           ENDIF

           ! Compute KMAH index                                   
           KMAHB = KMAHA 
           IF ( REAL( qB ) < 0.0 ) THEN 
              qAT = AIMAG( qA ) 
              qBT = AIMAG( qB ) 
              IF ( ( qAT < 0.0 .AND. qBT >= 0.0 ) .OR. &
                   ( qAT > 0.0 .AND. qBT <= 0.0 ) ) KMAHB = -KMAHA
           ENDIF

           ! Mode interpolation coefs                             
           IF ( .NOT. EXITED ) THEN 
              DB1 = xB * y1 - yB * x1 
              DB2 = xB * y2 - yB * x2 
              DB3 = xB * y3 - yB * x3 

              ! Compute areas of triangles                        
              A1 = (  D23 - DB3 + DB2 ) / delta 
              A2 = (  DB3 - D13 - DB1 ) / delta 
              A3 = ( -DB2 + DB1 + D12 ) / delta 

              ! Crossing into new element?                            
              ! SHOULD CHECK FOR A STEP WHICH CROSSES TWO ELTS!       
              IF ( ABS( A1 ) + ABS( A2 ) + ABS( A3 ) > 1.0 + EPSILON( A1 ) ) THEN 

                 ! Identify the side through which exitting       
                 IF ( A1 < 0.0 ) Iside = 2 
                 IF ( A2 < 0.0 ) Iside = 3 
                 IF ( A3 < 0.0 ) Iside = 1
                 IElement = ADJELT( Iside, IElement ) 

                 ! Normal elt transition or exit into free space? 
                 IF ( IElement == 0 ) THEN 
                    EXITED = .TRUE. 
                    Cx = 0.0 
                    Cy = 0.0 
                 ELSE 
                    CALL NEWELT( IElement, k, mode, M, maxM, Iset1, Iset2, Iset3, &
                         x1, y1, x2, y2, x3, y3, D12, D13, D23, delta, Cx, Cy, Mprop )
                    IF ( Mprop < mode ) EXIT ! If mode cuts off, skip to next mode
                 ENDIF
              ENDIF

              ! Evaluate modes at the rcvr depth                  
              CB   = A1 /   k( mode, Iset1 ) + A2 /   k( mode, Iset2 ) + A3 /   k( mode, Iset3 )          
              phiB = A1 * phi( mode, Iset1 ) + A2 * phi( mode, Iset2 ) + A3 * phi( mode, Iset3 )      
           ENDIF

           ! Compute beam influence
           CALL INFLU( xA-xS, yA-yS, xiA, etaA, pA, qA, tauA, CA, KMAHA, phiA, &
                       xB-xS, yB-yS, xiB, etaB, pB, qB, tauB, CB, KMAHB, phiB, RadVec, Ntheta, Rmin, Rmax, NR, const, P )

           xA    = xB
           yA    = yB
           xiA   = xiB
           etaA  = etaB
           pA    = pB
           qA    = qB
           tauA  = tauB
           CA    = CB
           KMAHA = KMAHB
           phiA  = phiB
                                                    
        ENDDO ! Next step

        ! Optionally dump rays to disk
        IF ( Option(6:6) == 'R' ) THEN 
           WRITE( RAYFIL, * ) Istep - 1 
           DO Istep = 1, Istep - 1 
              WRITE( RAYFIL, * ) xV( Istep ), yV( Istep ) 
           ENDDO
        ENDIF

     ENDDO   ! next mode 

  END DO ANGLE   ! next beam take-off angle

  RETURN 
END SUBROUTINE EVALGB
!**********************************************************************C
SUBROUTINE NEWELT( IElement, k, mode, M, maxM, Iset1, Iset2, Iset3, x1, y1, x2, y2, x3, y3, D12, D13, D23, delta, Cx, Cy, Mprop )

  ! Given elt number, returns info which is constant in elt           

  USE ElementMod
  IMPLICIT NONE  
  INTEGER          :: M( * ), node1, node2, node3, Iset1, Iset2, Iset3, Mprop, mode, maxM, IElement
  REAL (KIND=8)    :: D12, D13, D23, delta, x1, y1, x2, y2, x3, y3, A1x, A2x, A3x, A1y, A2y, A3y
  COMPLEX          :: k( maxM, * )
  COMPLEX (KIND=8) :: Cx, Cy

  node1 = node( 1, IElement ) 
  node2 = node( 2, IElement ) 
  node3 = node( 3, IElement ) 

  Iset1 = Iset( node1 ) 
  Iset2 = Iset( node2 ) 
  Iset3 = Iset( node3 ) 

  ! If mode cuts off, return to origin and do next mode           
  Mprop = MIN( M( Iset1 ), M( Iset2 ), M( Iset3 ) ) 
  IF ( mode > Mprop ) RETURN 

  x1 = x( node1 ) 
  y1 = y( node1 ) 
  x2 = x( node2 ) 
  y2 = y( node2 ) 
  x3 = x( node3 ) 
  y3 = y( node3 ) 

  D12   = x1 * y2 - y1 * x2 
  D13   = x1 * y3 - y1 * x3 
  D23   = x2 * y3 - y2 * x3 
  delta = D23 - D13 + D12 

  ! Gradient                                                      
  A1x = -y3 + y2 
  A2x =  y3 - y1 
  A3x = -y2 + y1 
  A1y =  x3 - x2 
  A2y = -x3 + x1 
  A3y =  x2 - x1 

  Cx = ( A1x / k( mode, Iset1 ) + A2x / k( mode, Iset2 ) + A3x / k( mode, Iset3 ) ) / delta        
  Cy = ( A1y / k( mode, Iset1 ) + A2y / k( mode, Iset2 ) + A3y / k( mode, Iset3 ) ) / delta        

  RETURN 
END SUBROUTINE NEWELT
!**********************************************************************C
SUBROUTINE INFLU( xA, yA, xiA, etaA, pA, qA, tauA, CA, KMAHA, phiA, &
                  xB, yB, xiB, etaB, pB, qB, tauB, CB, KMAHB, phiB, &
                  RadVec, Ntheta, Rmin, Rmax, NR, const, P )

  ! Computes contribution to receivers assuming beam cannot           
  ! contribute to a radial in an incoming sense                       

  IMPLICIT NONE
  REAL,    PARAMETER    :: BeamWindow = 5
  COMPLEX, PARAMETER    :: i = (0.0, 1.0)

  INTEGER,          INTENT( IN ) :: KMAHA, KMAHB, Ntheta, NR
  REAL,             INTENT( IN ) :: RadVec( 2, * ), Rmin, Rmax
  REAL    (KIND=8), INTENT( IN ) :: xA, yA, xiA, etaA, xB, yB, xiB, etaB
  COMPLEX (KIND=8), INTENT( IN ) :: const, phiA, phiB, pA, pB, qA, qB, tauA, tauB, CA, CB

  INTEGER            :: KMAH, Itheta, ir1, ir2, ir
  REAL    (KIND=8)   :: NA, NB, Nsquared, R, RA, RB, deltar, deltaA, deltaB, W, qAT, qBT 
  COMPLEX            :: P( Ntheta, * )
  COMPLEX (KIND=8)   :: phiMid, pMid, qMid, tauMid, CMid, contrib

  deltaR = ( Rmax - Rmin ) / ( NR - 1 ) 

  ! Loop over radials of receiver line and compute contribution

  DO Itheta = 1, Ntheta

     ! Compute intercept range, ra, & index preceding rcvr        
     deltaA = RadVec( 1, Itheta ) * xiA + RadVec( 2, Itheta ) * etaA 
     deltaB = RadVec( 1, Itheta ) * xiB + RadVec( 2, Itheta ) * etaB

     IF ( ABS( deltaA ) < TINY( deltaA ) ) CYCLE
     IF ( ABS( deltaB ) < TINY( deltaB ) ) CYCLE

     RA  = ( yA * etaA + xA * xiA ) / deltaA
     RB  = ( yB * etaB + xB * xiB ) / deltaB

     ir1 = MAX( MIN( INT( ( RA - Rmin ) / deltaR ) + 1, Nr ), 0 )
     ir2 = MAX( MIN( INT( ( RB - Rmin ) / deltaR ) + 1, Nr ), 1 )

     ! If a receiver is bracketted, compute influence
     IF ( ir2 > ir1 .AND. deltaA * deltaB > 0 ) THEN

        ! Normal distance                                         
        NA = ( xA * RadVec( 2, Itheta ) - yA * RadVec( 1, Itheta ) ) / ( DBLE( CA ) * deltaA )     
        NB = ( xB * RadVec( 2, Itheta ) - yB * RadVec( 1, Itheta ) ) / ( DBLE( CB ) * deltaB )     

        DO ir = ir1+1, ir2 
           R         = Rmin + ( ir - 1 ) * deltaR 
           W         = ( R - RA ) / ( RB - RA ) 
           pMid      =   pA + W * ( pB - pA ) 
           qMid      =   qA + W * ( qB - qA ) 
           Nsquared  = ( NA + W * ( NB - NA ) ) ** 2

           ! Within beam window?                               
           IF ( -0.5 * AIMAG( pMid / qMid ) * Nsquared < BeamWindow ) THEN
              CMid   = CA   + W * ( CB   - CA   ) 
              tauMid = tauA + W * ( tauB - tauA ) 
              phiMid = phiA + W * ( phiB - phiA )

              ! Compute KMAH index                                
              KMAH = KMAHA 
              IF ( REAL( qMid ) < 0.0 ) THEN 
                 qAT = AIMAG( qA ) 
                 qBT = AIMAG( qMid ) 
                 IF ( ( qAT < 0.0 .AND. qBT >= 0.0 ) .OR. ( qAT > 0.0 .AND. qBT <= 0.0 ) ) KMAH = -KMAH                
              END IF

              contrib = const * phiMid * SQRT( CMid / qMid ) * EXP( -i * ( tauMid + 0.5 * pMid / qMid * Nsquared ) )           
              IF ( KMAH < 0 ) contrib = -contrib
              P( Itheta, ir ) = P( Itheta, ir ) + CMPLX( contrib ) 
           END IF
        END DO ! Next receiver on radial
     END IF

  END DO   ! Next radial

  RETURN 
END SUBROUTINE INFLU
