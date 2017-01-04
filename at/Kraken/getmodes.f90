SUBROUTINE GetModes( FileRoot, IProf, MaxM, rd, Nrd, Comp, k, PhiR, M, Freq, Title )                                     

  ! Read in modes and extract values at rcvr depths                   

  ! INPUT:
  !    FileRoot for the mode set
  !    IProf    unit number  for the mode set
  !    FileName name of file for the mode set
  !    MaxM     row dimension of PhiR in calling program
  !    rd       vector of receiver depths where modes are to be evaluate
  !    Nrd      number of such receivers
  !    Comp     component (vertical, horizontal,...) for elastic problem
  !                 ignored for purely acoustic problems

  ! OUTPUT:                                                           
  !    k      vector of eigenvalues                                   
  !    PhiR   matrix of tabulated modes                               
  !    M      number of modes
  !    Freq   frequency                                               
  !    Title  title in mode file

  ! If IProf = 1 then the file is opened and IRecProfile is set to the first record          
 
  IMPLICIT NONE
  INTEGER, PARAMETER :: PRTFile = 6, MaxN = 16001, MaxMedium = 51
  INTEGER            :: N( MaxMedium ), Ird( Nrd ), ir, iz, IProf, IRecProfile, LRecl, M, MaxM, Mode, NMat, NMedia, Nrd, NTot
  REAL               :: rd( Nrd ), Z( MaxN ), W( Nrd ), Depth( MaxMedium ), rho( MaxMedium ), Tol, freq, &
                        rhoT, rhoB, csT, csB, depthT, depthB, wT
  COMPLEX            :: PhiR( MaxM, Nrd ), PhiT( MaxN ), k( * ), kTop2, kBot2
  CHARACTER          :: Title*( * ), Comp*( * )
  CHARACTER (LEN=*)  :: FileRoot
  CHARACTER (LEN=8)  :: Material( MaxMedium )
  CHARACTER (LEN=1)  :: BCTop, BCBot

  SAVE iRecProfile

  ! Read the header data from the mode file
  CALL ModeHeader( FileRoot, IProf, IRecProfile, LRecl, Title, Freq, NMedia, NTot, NMat, N, Material, Depth, rho,  &
       BCTop, rhoT, DepthT, BCBot, rhoB, DepthB, M, MaxM, Z, k, kTop2, kBot2 )

  IF ( M <= 0 ) RETURN  ! no modes? quick return.
  CALL WEIGHT( Z, NTot, rd, Nrd, W, Ird )   ! Locate indices of receiver points

  ! Loop over receiver depths to check for safe interpolation
  ! Receivers must be within a fraction of a wavelength
  ! of tabulated pts. We accept one wavelength at 1500 m/s

  Tol = 1500.0 / Freq 

  DO ir = 1, Nrd 

     iz = Ird( ir ) 
     WT = ABS( MIN( W( ir ), 1.0 - W( ir ) ) ) 

     IF ( rd( ir ) < DepthT ) THEN        ! Rcvr in upper halfspace
        cST = 0.0   ! should be passed by ModeHeader
        IF ( cST /= 0.0 .OR. BCTop(1:1) /= 'A' ) THEN 
           WRITE( PRTFile, * ) 'Receiver depth: ', rd( ir )
           WRITE( PRTFile, * ) 'Highest valid depth: ', DepthT
           CALL ERROUT( PRTFile, 'F', 'GetMode', 'Rcvr above highest valid depth' )
        ENDIF

     ELSE IF ( rd( ir ) > DepthB ) THEN   ! Rcvr in lower halfspace
        cSB = 0.0   ! should be passed by ModeHeader
        IF ( cSB /= 0.0 .OR. BCBot(1:1) /= 'A' ) THEN 
           WRITE( PRTFile, * ) 'Receiver depth: ', rd( ir )
           WRITE( PRTFile, * ) 'Lowest valid depth: ', DepthB 
           CALL ERROUT( PRTFile, 'F', 'GetMode', 'Rcvr below lowest valid depth' )
        ENDIF

     ELSE IF ( NTot > 1 ) THEN            ! Rcvr between two grid points or large extrapolation
        IF ( WT * ( Z( iz + 1 ) - Z( iz ) ) > Tol ) THEN
           WRITE( PRTFile, * ) 'Receiver depth: ', rd( ir )
           WRITE( PRTFile, * ) 'Nearest depths: ', Z( iz ), Z( iz + 1 )
           WRITE( PRTFile, * ) 'Tolerance: ', Tol 
           CALL ERROUT( PRTFile, 'F', 'GetMode', 'Modes not tabulated near requested pt.' )
        ENDIF
     ELSE                                 ! Rcvr near a single grid point
        IF ( ABS( rd( ir ) - Z( iz ) ) > Tol ) THEN 
           WRITE( PRTFile, * ) 'Rd, Tabulation depth ', rd( ir), Z( iz )
           WRITE( PRTFile, * ) 'Tolerance: ', Tol 
           CALL ERROUT( PRTFile, 'F', 'GetMode', 'Modes not tabulated near requested pt.' )
        ENDIF
     ENDIF

  END DO   ! next receiver depth

  ! Read in the modes
  DO Mode = 1, M
     CALL GETONE( Mode, IRecProfile, NTot, NMat, W, Ird, N, Material, NMedia, Comp, &
          kTop2, DepthT, BCTop, kBot2, DepthB, BCBot, rd, Nrd, k, PhiT )
     PhiR( Mode, 1:Nrd ) = PhiT( 1:Nrd )
  END DO   ! next mode

  IRecProfile = IRecProfile + 6 + M + 1 + ( 2 * M - 1 ) / LRecL   ! advance to next profile

END SUBROUTINE GetModes

!**********************************************************************C

SUBROUTINE ModeHeader( FileRoot, IProf, IRecProfile, LRecl, Title, Freq, NMedia, NTot, NMat, N, Material, Depth, rho, &
     &   BCTop, rhoT, DepthT, BCBot, rhoB, DepthB, M, MaxM, Z, k, kTop2, kBot2 )

  ! Reads the header information from ModeFile                          
  ! Note T suffix means top
  !      B suffix means bottom                                        

  ! IProf     is a profile number
  ! FileRoot  is the user-provided file name                             
  ! FileNameT is the temporary name we build
  ! These have to be two separate variables to ensure there
  ! is space allocated to construct the file name even when           
  ! the user calls us with FileName = ' ' 

  ! IRecProfile must point to the first record of the profile                             

  IMPLICIT NONE
  INTEGER, PARAMETER :: PRTFile = 6, ModeFile = 30 
  REAL,    PARAMETER :: pi = 3.141592
  LOGICAL            :: OpenFlag
  INTEGER            :: N( * ) , NMedia, NTot, NMat, M, MaxM, IFirst, ILast, iProf, &
                        iostat, IRecProfile, IRec, LRecL, Medium
  REAL               :: Depth( * ), Z( * ), rho( * ), Freq, rhoT, DepthT, rhoB, DepthB 
  COMPLEX            :: k( * ), cPT, cST, cPB, cSB, kTop2, kBot2
  CHARACTER          :: Title*( * )
  CHARACTER (LEN=80) :: FileNameT
  CHARACTER (LEN=* ) :: FileRoot
  CHARACTER (LEN=8 ) :: Material( * )
  CHARACTER (LEN=1 ) :: BCTop, BCBot

  ! open ModeFile
  FileNameT = TRIM( FileRoot ) // '.mod'

  INQUIRE( FILE = FileNameT, OPENED = OpenFlag )
  IF ( .NOT. OpenFlag ) THEN
     OPEN( UNIT = ModeFile, FILE = FileNameT, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', &
          RECL = 100, IOSTAT = iostat )
     IF ( IOSTAT /= 0 ) THEN
        WRITE( PRTFile, * ) 'Mode file = ', FileNameT
        CALL ERROUT( PrtFile, 'F', 'GetMode - ModeHeader', 'Unable to open the mode file' )
     END IF

     READ( ModeFile, REC = 1 ) LRecL
     CLOSE( UNIT = ModeFile )
     OPEN(  UNIT = ModeFile, FILE = FileNameT, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', &
          RECL = 4 * LRecL, IOSTAT = iostat )
  END IF

  write( *, * ) 'iprof = ', iprof
  ! If this is the first profile, reset the record counter to the beginning of the file
  IF ( IProf == 1 ) THEN
     IRecProfile = 1   ! set counter pointing to the first record to read
  END IF

  ! Read header info
  READ( ModeFile, REC = IRecProfile     ) LRecL, Title( 1 : 80 ), Freq, NMedia, NTot, NMat  

  write( *, * ) Title                
  READ( ModeFile, REC = IRecProfile + 1 ) ( N( Medium ), Material( Medium ), Medium = 1, NMedia)
  READ( ModeFile, REC = IRecProfile + 2 ) BCTop(1:1), cPT, cST, rhoT, DepthT, BCBot(1:1), cPB, cSB, rhoB, DepthB
  READ( ModeFile, REC = IRecProfile + 3 ) ( Depth( Medium ), rho( Medium ), Medium = 1, NMedia )
  READ( ModeFile, REC = IRecProfile + 4 ) M, LRecL

  IF ( M > MaxM ) THEN 
     WRITE( PRTFile, * ) 'M = ', M, '   MaxM = ', MaxM
     CALL ERROUT( PRTFile, 'F', 'ModeHeader', 'Insufficient storage to read all the modes: increase MaxM' )
  ENDIF

  IF ( M == 0 ) RETURN 
  READ( ModeFile, REC = IRecProfile + 5 ) Z( 1 : NTot ) 

  ! Read in eigenvalues, k( I )
  ! They come at the end of the ModeFile, because they're calculated after the eigenvectors
  IFirst = 1 
  DO IRec = 1, 1 + ( 2 * M - 1 ) / LRecL
     ILast  = MIN( M, IFirst + LRecL / 2 - 1 )
     READ( ModeFile, REC = IRecProfile + 5 + M + IRec ) k( IFirst : ILast )      
     IFirst = ILast + 1 
  END DO

  IF ( BCTop( 1 : 1 ) == 'A' ) kTop2 = ( 2.0 * pi * Freq / cPT ) **2 
  IF ( BCBot( 1 : 1 ) == 'A' ) kBot2 = ( 2.0 * pi * Freq / cPB ) **2

END SUBROUTINE ModeHeader

!**********************************************************************C

SUBROUTINE GETONE( Mode, IRecProfile, NTot, NMat, W, Ird, N, Material, NMedia, Comp, &
     &   kTop2, DepthT, BCTop, kBot2, DepthB, BCBot, rd, Nrd, k, PhiR )

  ! Read in a single eigenfunction and extract receiver values
  ! Results are returned in PhiR

  IMPLICIT NONE
  INTEGER,PARAMETER :: ModeFile = 30
  LOGICAL           :: TufLuk 
  INTEGER           :: N( * ), Ird( * ), ir, iz, IRecProfile, j, Mode, NTot, NMat, NMedia, Nrd
  REAL              :: rd( * ), W( * ), DepthT, DepthB
  COMPLEX           :: PhiR( * ), k( * ), Phi( NMat ), gammaT, gammaB, kTop2, kBot2
  COMPLEX  (KIND=8) :: PekerisRoot, gamma2
  CHARACTER         :: Comp*( *)
  CHARACTER (LEN=8) :: Material( * )
  CHARACTER (LEN=1) :: BCTop, BCBot

  READ( ModeFile, REC = IRecProfile + 5 + Mode ) ( Phi( j ), j = 1, NMat ) 

  ! Is there an elastic medium in the problem?
  TufLuk = .FALSE. 
  IF ( ANY( Material( 1:NMedia ) == 'ELASTIC' ) ) TufLuk = .TRUE.

  ! Extract the component specified by 'Comp'
  IF ( TufLuk ) CALL EXTRACT( Phi, N, Material, NMedia, Comp ) 

  ! Extract values at receiver depths
  gammaT = 0.0
  gammaB = 0.0 

  ! n.b. should be using real( k(mode) ) for KRAKEN
  IF ( BCTop(1:1) == 'A' ) THEN 
     gamma2 = k( Mode ) ** 2 - kTop2
     gammaT = CMPLX( PekerisRoot( gamma2 ) )
  END IF

  IF ( BCBot(1:1) == 'A' ) THEN 
     gamma2 = k( Mode ) ** 2 - kBot2
     gammaB = CMPLX( PekerisRoot( gamma2 ) )
  END IF

  DO ir = 1, Nrd 
     IF ( rd( ir ) < DepthT ) THEN      ! Rcvr in upper halfspace
        PhiR( ir ) = Phi( 1    ) * EXP( -gammaT * ( DepthT - rd( ir  ) ) )

     ELSE IF ( rd( ir ) > DepthB ) THEN ! Rcvr in lower halfspace
        PhiR( ir ) = Phi( NTot ) * EXP( -gammaB * ( rd( ir ) - DepthB ) )

     ELSE IF ( NTot > 1 ) THEN 
        iz         = Ird( ir )
        PhiR( ir ) = Phi( iz ) + W( ir ) * ( Phi( iz + 1 ) - Phi( iz ) )

     ELSE                               ! mode is tabulated at only one depth
        iz = Ird( ir ) 
        PhiR( ir ) = Phi( iz ) 
     ENDIF

  END DO

END SUBROUTINE GETONE
!**********************************************************************C
SUBROUTINE EXTRACT( Phi, N, Material, NMedia, Comp ) 

  ! For elastic problems where a stress-displacement vector is output,
  ! extracts the desired component                                    

  IMPLICIT NONE
  INTEGER,           INTENT( IN    ) :: N( * ), NMedia
  INTEGER                            :: i, j, k, Medium
  COMPLEX,           INTENT( INOUT ) :: Phi( * ) 
  CHARACTER,         INTENT( IN    ) :: Comp*( *)
  CHARACTER (LEN=8), INTENT( IN    ) :: Material( * )
 
  j = 1 
  k = 1 

  DO Medium = 1, NMedia 

     DO i = 1, N( Medium ) + 1
        SELECT CASE ( Material( Medium ) )
        CASE ( 'ACOUSTIC' )
           Phi( j ) = Phi( k )
           k = k + 1
        CASE ( 'ELASTIC' )
           SELECT CASE ( Comp )
           CASE ( 'H' )
              Phi( j ) = Phi( k     )
           CASE ( 'V' )
              Phi( j ) = Phi( k + 1 )
           CASE ( 'T' )
              Phi( j ) = Phi( k + 2 )
           CASE ( 'N' )
              Phi( j ) = Phi( k + 3 )
           END SELECT

           k = k + 4
        END SELECT   ! Comp
        j = j + 1 
     END DO   ! Material

  END DO   ! next medium

END SUBROUTINE EXTRACT
