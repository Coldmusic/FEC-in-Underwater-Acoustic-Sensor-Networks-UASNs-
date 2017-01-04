SUBROUTINE EVALCM( FileRoot, RProf, NProf, phiS, phi, rd, Nrd, R, Nr, k, M, Opt, P )                                         

  ! Computes pressure field using coupled mode theory                 
  ! Normalized to pressure of point source at 1 meter

  ! Opt = X     Cartesian   (X, z) coordinates                        
  ! Opt = R     Cylindrical (R, z) coordinates

  ! Note number of propagating modes is reset after first segment.    
  ! Thus M restricts the number of modes in the source field but
  ! thereafter energy is allowed to couple into higher-order modes.

  IMPLICIT NONE
  INTEGER, PARAMETER :: MaxM = 20000
  REAL,    PARAMETER :: pi = 3.1415926
  COMPLEX, PARAMETER :: i = ( 0.0, 1.0 )
  LOGICAL            :: first
  INTEGER            :: NProf, iProf, M1, M, Nr, Nrd, ir, ird
  REAL               :: rd( Nrd ), RProf( NProf + 1 ), r( Nr ), freq
  COMPLEX            :: phiS( * ), phi( MaxM, * ), A( MaxM ), sum, P( Nrd, * ), k( MaxM )
  CHARACTER (LEN=50) :: Opt
  CHARACTER (LEN=80) :: Title, FileRoot
  SAVE first
  DATA first /.TRUE./

  ! Compute ranges (in meters) where new profiles are used    
  IF ( first ) THEN  

     DO iProf = NProf, 2, -1 
        RProf( iProf ) = 500.0 * ( RProf( iProf ) + RProf( iProf-1 ) )
     END DO

     RProf(  NProf + 1 ) = HUGE( RProf( NProf ) )
     first = .FALSE.
  ENDIF

  ! Evaluate mode excitation coefficients, A(mode)            
  iProf = 1

  ! THIS HAS ALREADY BEEN READ IN PLOTTLR!                            
  CALL GetModes( FileRoot, iProf, MaxM, rd, Nrd, 'N', k, phi,  M1, freq, Title )              

  M = MIN( M, M1 )                                                    
  IF ( Opt( 1 : 1 ) == 'X' ) THEN   ! Cartesian coordinates
     A( 1 : M ) =     SQRT( 2.0 * pi ) * EXP( i * pi / 4.0 ) * phiS( 1 : M ) /       k( 1 : M )
  ELSE                          ! Cylindrical coordinates
     A( 1 : M ) = i * SQRT( 2.0 * pi ) * EXP( i * pi / 4.0 ) * phiS( 1 : M ) / SQRT( k( 1 : M ) )
  ENDIF

  ! March forward in range                                    
  DO ir = 1, Nr 
     ! WRITE( *, * ) 'Range = ', r( ir )
     IF ( r( ir ) > RProf( iProf + 1 ) ) THEN  ! Crossing into new range segment?

        iProf = iProf + 1 

        ! Advance to interface
        IF ( ir == 1 ) THEN   ! first range
           A( 1 : M ) = A( 1 : M ) * EXP( -i * k( 1 : M ) *   RProf( iProf ) )
        ELSE
           A( 1 : M ) = A( 1 : M ) * EXP( -i * k( 1 : M ) * ( RProf( iProf ) - r( ir - 1 ) ) )
        ENDIF

        ! Here's where we cross over
        IF ( iProf  <= NProf ) THEN 
           CALL NewProfile( FileRoot, iProf, k, phi, M, rd, Nrd, A )
           WRITE( *, * ) 'New profile read at receiver range', r(ir ) / 1000.0, iProf, ' #modes = ', M
        ENDIF

        ! Are there other segments to cross? Advance phase through each segment
        DO WHILE ( r( ir ) > RProf( iProf + 1 ) )
           iProf    = iProf + 1
           A( 1 : M ) = A( 1 : M ) * EXP( -i * k( 1 : M ) * ( RProf( iProf ) - RProf( iProf - 1 ) ) )

           IF ( iProf <= NProf ) CALL NewProfile( FileRoot, iProf, k, phi, M, rd, Nrd, A )
        END DO

        ! Advance the remaining distance past the last interface
        A( 1 : M ) = A( 1 : M ) * EXP( -i * k( 1 : M ) * ( r( ir ) - RProf( iProf ) ) )

     ELSE  ! no new segment, just advance the phase based on the range step
        IF ( ir == 1 ) THEN   ! first range
           A( 1 : M ) = A( 1 : M ) * EXP( -i * k( 1 : M ) *   r( ir )  )
        ELSE
           A( 1 : M ) = A( 1 : M ) * EXP( -i * k( 1 : M ) * ( r( ir ) - r( ir - 1 ) ) )
        ENDIF

     ENDIF

     ! For each rcvr add up modal contributions
     DO ird = 1, Nrd              
        IF ( Opt( 1 : 1 ) == 'R' .AND. r( ir ) /= 0.0 ) THEN
           P( ird, ir ) = SUM( A( 1 : M ) * phi( 1 : M, ird ) ) / SQRT( r( ir ) )
        ELSE
           P( ird, ir ) = SUM( A( 1 : M ) * phi( 1 : M, ird ) ) 
        ENDIF
     END DO

  END DO    ! next range step 

END SUBROUTINE EVALCM

!**********************************************************************C

SUBROUTINE NewProfile( FileRoot, iProf, k, phiR, MR, rd, Nrd, A ) 

  ! For a given profil number:                                        
  !     read in modes for current segment                             
  !     project the pressure field onto the new modes
  !     extract values of the modes at rcvr depths

  IMPLICIT NONE
  INTEGER, PARAMETER   :: MaxM = 20000,  MaxN = 16001, ModeFile = 30
  INTEGER              :: ird( Nrd ), Mode, Nr, Nrd, NTot, ML, MR, iProf, ir, iz, IRecProfile
  REAL                 :: z(  MaxN ), rd( * ), W( Nrd ), depthTL, depthBL, depthTR, depthBR, rhoBR, rhoTR, zt
  COMPLEX              :: P(  MaxN ), k( * ), phiR( MaxM, * ), A( * ), sum1, &
                          gamTL( MaxM ), gamBL( MaxM ), phiTL( MaxM ), phiBL( MaxM ),    &
                          gamTR, gamBR, phiTR, phiBR, KTop2R, KBot2R, tail
  COMPLEX, ALLOCATABLE :: phi( : ), phiTmp( : )
  COMPLEX     (KIND=8) :: PekerisRoot, gamma2
  CHARACTER   (LEN= 1) :: BCBotR, BCTopR
  CHARACTER   (LEN=80) :: FileRoot
  SAVE IRecProfile

  ! Compute pressure along the left of the interface
  CALL PLEFT( FileRoot, iProf, IRecProfile, A, k, z, MR, P, NR, NTot,                &
       &   BCTopR, rhoTR, KTop2R, depthTR, BCBotR, rhoBR, KBot2R, depthBR, &
       &   gamTL, gamBL, depthTL, depthBL, phiTL, phiBL, ML, MaxM )

  !  Read in eigenfunctions and extract receiver values       
  CALL WEIGHT( z, NTot, rd, Nrd, W, ird ) ! Compute weights for mode interpolation at rcvr depths
  ALLOCATE( phi( NR ), phiTmp( NTot ) )

  DO mode = 1, MR
     READ( ModeFile, REC = IRecProfile + 5 + mode ) ( phi( iz ), iz = 1, Nr )   ! read in the mode

     IF ( BCTopR == 'A' ) THEN 
        phiTR  = phi( 1 )
        gamma2 = k( mode ) ** 2 - KTop2R
        gamTR  = CMPLX( PekerisRoot( gamma2 ) )
     ENDIF

     IF ( BCBotR == 'A' ) THEN 
        phiBR  = phi( Nr )
        gamma2 = k( mode ) ** 2 - KBot2R
        gamBR  = CMPLX( PekerisRoot( gamma2 ) )
     ENDIF

     ! tabulate the new mode on the grid from the previous segment
     DO iz = 1, NTot
        zT = z( iz ) 
        IF      ( zT > depthBR ) THEN
           IF ( BCBotR == 'A' ) phiTmp( iz ) = phiBR * EXP( -gamBR * ( zT - depthBR ) )
        ELSE IF ( zT < depthTR ) THEN
           IF ( BCTopR == 'A' ) phiTmp( iz ) = phiTR * EXP( -gamTR * ( depthTR - zT ) )
        ELSE
           phiTmp( iz ) = phi( iz )
        ENDIF

     END DO

     ! Compute new amplitudes:
     !       A = Integral[ P( z ) * phi( z ) dz ]
     !       (Repeat for each mode phi to produce excitation coef. A)
     !       Integral is done using trapezoidal rule

     sum1 = SUM( P( 1 : NTot ) * phiTmp )

     IF ( BCTopR == 'A' ) THEN  ! contribution from upper halfspace
        sum1   = sum1 + tail( z( 1    ), phiTL, gamTL, depthTL, ML, phiTR / rhoTR, gamTR, depthTR )
     ENDIF

     IF ( BCBotR == 'A' ) THEN  ! contribution from lower halfspace
        sum1   = sum1 + tail( z( NTot ), phiBL, gamBL, depthBL, ML, phiBR / rhoBR, gamBR, depthBR )
     ENDIF

     IF ( mode > ML ) A( mode ) = 0.0 
     A( mode ) = sum1

     ! Subtabulate modes at receiver depths
     DO ir = 1, Nrd
        phiR( mode, ir ) = 0.0 
        IF      ( rd( ir ) < depthTR ) THEN ! Rcvr in upper halfspace
           IF ( BCTopR == 'A' ) phiR( mode, ir ) = phiTR * EXP( -gamTR * ( depthTR   - rd( ir  ) ) )
        ELSE IF ( rd( ir ) > depthBR ) THEN ! Rcvr in lower halfspace
           IF ( BCBotR == 'A' ) phiR( mode, ir ) = phiBR * EXP( -gamBR * ( rd( ir )  - depthBR   ) )
        ELSE
           iz = ird( ir ) 
           phiR( mode, ir ) = phi( iz ) + W( ir ) * ( phi( iz + 1 ) - phi( iz ) )
        ENDIF
     END DO

  END DO    ! next mode

  WRITE( *, * ) 'depth-averaged power:', SUM( ABS( A( 1 : MR ) ) ** 2 )

END SUBROUTINE NewProfile

!**********************************************************************C

SUBROUTINE PLEFT( FileRoot, iProf, IRecProfile, A, k, z, M, P, NR, NTot,           &
     &   BCTop, rhoT, KTop2, depthT, &
     &   BCBot, rhoB, KBot2, depthB,        &
     &   gamTL, gamBL, depthTL, depthBL, phiTL, phiBL, ML, MaxM )

  ! Computes the pressure field along the interface sampled on the grid of the new segment                 
  ! Also returns information needed for the tails in the halfspaces   

  IMPLICIT NONE
  INTEGER,   PARAMETER :: MaxN = 16001, Maxmed = 50, ModeFile = 30
  INTEGER              :: N( Maxmed ), IRecProfile, MaxM, M, ML, NL, NR, NTot, NMat, med, mode, NMedia, iz, izL, LRecL, iProf
  REAL                 :: freq, zL( MaxN ), z( * ), zt, rhoL( Maxmed ), rho( Maxmed ), depthL( Maxmed ), depth( Maxmed ), &
                          depthT, depthB, depthTL, depthBL, DBelow, rhoT, rhoB, rhoMed, rhoBel, h
  COMPLEX              :: k( * ),  P( * ), A( * ), gamTL( * ), gamBL( * ), phiTL( * ), phiBL( * ), KTop2, KBot2
  COMPLEX, ALLOCATABLE :: phi( : ), PL( : )
  COMPLEX     (KIND=8) :: PekerisRoot, gamma2
  CHARACTER   (LEN=80) :: Title, FileRoot
  CHARACTER   (LEN= 8) :: Material( Maxmed )
  CHARACTER   (LEN= 1) :: BCBot, BCTop

  ! Read modal info at end of last segment                    
  CALL ModeHeader( FileRoot, iProf - 1, IRecProfile, LRecl, Title, freq, Nmedia, NL, NMat, N, Material, depthL, rhoL, &
       BCTop, rhoT, depthTL,  BCBot, rhoB, depthBL, ML, MaxM, zL, k, KTop2, KBot2 )

  ML = M   ! we only know M amplitudes
!! should be min( ML, M )

  ! Compute pressure at the interface
  ALLOCATE( phi( NL ), PL( NL ) )
  PL = 0.0

  DO mode = 1, ML 
     READ( ModeFile, REC = IRecProfile + 5 + mode ) phi

     PL = PL + A( mode ) * phi

     ! Halfspace information
     phiTL( mode ) = A( mode ) * phi( 1  )
     phiBL( mode ) = A( mode ) * phi( NL )                       
     gamTL( mode ) = 0.0
     gamBL( mode ) = 0.0

     IF ( BCTop == 'A' ) THEN   ! top halfspace
        gamma2        = k( mode ) ** 2 - KTop2
        gamTL( mode ) = CMPLX( PekerisRoot( gamma2 ) )
     END IF

     IF ( BCBot == 'A' ) THEN   ! bottom halfpsace
        gamma2        = k( mode ) ** 2 - KBot2
        gamBL( mode ) = CMPLX( PekerisRoot( gamma2 ) )
     END IF

  END DO

  ! set record pointer to beginning of next mode set
  IRecProfile = IRecProfile + 6 + M + 1 + ( 2 * M - 1 ) / LRecL

  ! Read modal data in new segment
  CALL ModeHeader( FileRoot, iProf, IRecProfile, LRecl, Title, freq, Nmedia, NR, NMat, N, Material, depth, rho, &
       &   BCTop, rhoT, depthT, BCBot, rhoB, depthB, M, MaxM, z, k, KTop2, KBot2 )

  IF ( z( 1 ) /= depthT .OR. z( NR ) /= depthB ) THEN
     WRITE( *, * ) 'Fatal Error: modes must be tabulated throughout the ocean and sediment to compute the coupling coefs.'
     STOP
  END IF

  ! Upslope? Extend the z vector with data from zL
  ! This code should be generalized to include bottom and top cases
  NTot = NR 
  DO izL = 1, NL 
     IF ( zL( izL ) > z( NTot ) ) THEN
        NTot      = NTot + 1
        z( NTot ) = zL( izL )
     ENDIF
  END DO

  ! Retabulate the pressure on the new grid
  izL    = 1 
  med    = 1
  rhomed = rho( 1 )

  ! Depth of next interface below this one
  IF ( med < Nmedia ) THEN 
     DBelow = depth( med + 1 )
     rhoBel = rho(   med + 1 )
  ELSE
     DBelow = depthB 
     rhoBel = rhoB
  ENDIF

  DO iz = 1, NTot 
     zT = z( iz )

     ! Get medium density
     IF      ( zT < depthT ) THEN
        rhomed = rhoT
     ELSE IF ( zT > depthB ) THEN
        rhomed = rhoB
     ELSE IF ( med < Nmedia ) THEN
        IF   ( zT > DBelow ) THEN
           med    = med + 1
           rhomed = rho( med )
        ENDIF
     ENDIF

     ! depth of next interface Below this one
     IF ( med < Nmedia ) THEN
        DBelow = depth( med + 1 )
        rhoBel = rho(   med + 1 )
     ELSE
        DBelow = depthB
        rhoBel = rhoB 
     ENDIF

     DO WHILE ( zT > zL( izL + 1 ) .AND. izL < NL - 1 ) 
        izL = izL + 1
     END DO

     ! Calculate P at that depth
     IF      ( zT > depthBL ) THEN       ! lower halfspace
        IF ( BCBot == 'A' ) THEN
           P( iz ) = SUM( phiBL( 1 : ML ) * EXP( -gamBL( 1 : ML ) * ( zT - depthBL ) ) )
        ENDIF

     ELSE IF ( zT < depthTL ) THEN       ! upper halfpace
        IF ( BCTop == 'A' ) THEN
           P( iz ) = SUM( phiTL( 1 : ML ) * EXP( -gamTL( 1 : ML ) * ( depthTL - zT ) ) )
        ENDIF
     ELSE 
        P( iz ) = PL( izL ) +                                       &
             ( zT - zL( izL ) ) / ( zL( izL + 1 ) - zL( izL ) ) *     &
                                  ( PL( izL + 1 ) - PL( izL ) )
     ENDIF

     ! compute mesh width, h
     IF ( iz == 1 ) THEN         ! first point
        h = 0.5 * ( z(  2   ) - z(   1      ) ) / rhomed
     ELSE IF ( iz == NTot ) THEN ! Last point
        h = 0.5 * ( z( NTot ) - z( NTot - 1 ) ) / rhomed
     ELSE                        ! Point just above or below the interface
        IF ( z( iz - 1 ) < DBelow .AND. z( iz + 1 ) > DBelow ) THEN                         
           h = 0.5 * ( z( iz + 1 ) / rhoBel - z( iz - 1 ) / rhomed &
                &         - DBelow / rhoBel +      DBelow / rhomed )
        ELSE
           h = 0.5 * ( z( iz + 1 ) - z( iz   - 1 ) ) / rhomed
        ENDIF
     ENDIF

     P( iz ) = h * P( iz ) 

  END DO   ! next depth iz

END SUBROUTINE PLEFT
!**********************************************************************C
COMPLEX FUNCTION tail( D, phiL, gamL, DL, ML, phiR, gamR, DR )

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ML
  REAL,    INTENT(IN) :: D, DL, DR
  COMPLEX, INTENT(IN) :: gamL( ML ), phiL( ML ), gamR, phiR
  COMPLEX             :: FR

  FR =  phiR * EXP( -gamR * ( D - DR ) ) 

  IF ( D == DL ) THEN 
     tail = FR * SUM( phiL                             / ( gamL + gamR ) )
  ELSE
     tail = FR * SUM( phiL * EXP( -gamL * ( D - DL ) ) / ( gamL + gamR ) )
  ENDIF

END FUNCTION tail
