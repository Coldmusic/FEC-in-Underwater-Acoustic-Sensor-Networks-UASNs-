PROGRAM MODBIN 

  ! Converts asciii mode file to binary format                        

  INTEGER,     PARAMETER :: BINFIL = 20, ASCFIL = 5
  COMPLEX                   CPT, CST, CPB, CSB
  CHARACTER     (LEN=80) :: TITLE
  CHARACTER     (LEN=1 ) :: BCTOP, BCBOT
  INTEGER,   ALLOCATABLE :: N( : ) 
  REAL,      ALLOCATABLE :: Z(: ), RHO( : ), depth( : )
  COMPLEX,   ALLOCATABLE :: PHI( : ), k( : )
  CHARACTER (LEN=8), ALLOCATABLE :: Material( : )

  ! Read header

  READ( ASCFIL, * ) LRECL 
  READ( ASCFIL, * ) TITLE 
  READ( ASCFIL, * ) FREQ, Nmedia, NTOT, NMAT, M 

  ALLOCATE( N( Nmedia), Material( Nmedia), depth( Nmedia ), rho( Nmedia ), Z( NTOT ), PHI( NTOT ), k( M ) )

  DO medium = 1, Nmedia 
     READ( ASCFIL, * ) N( medium ), depth( medium ), rho( medium ), Material( medium )                     
  ENDDO

  READ( ASCFIL, * ) BCTOP(1:1), CPT, CST, RHOT, depthT 
  READ( ASCFIL, * ) BCBOT(1:1), CPB, CSB, rhoB, depthB 

  DO IZ = 1, NTOT 
     READ( ASCFIL, * ) Z( IZ ) 
  END DO

  ! Write header

  OPEN ( FILE = 'MODFIL', UNIT = BINFIL, STATUS = 'NEW', ACCESS = 'DIRECT', RECL = 4 * LRECL, FORM = 'UNFORMATTED' )    

  WRITE( BINFIL, REC = 1 ) LRECL, TITLE(1:80), FREQ, Nmedia, NTOT, NMAT                  
  WRITE( BINFIL, REC = 2 ) ( N( medium ), Material( medium ), medium = 1, Nmedia)                
  WRITE( BINFIL, REC = 3 ) BCTOP(1:1), CPT, CST, rhoT, depthT, BCBOT(1:1), CPB, CSB, rhoB, depthB        
  WRITE( BINFIL, REC = 4 ) ( depth( medium ), rho( medium ), medium = 1, Nmedia )               
  WRITE( BINFIL, REC = 5 ) M, LRECL 
  WRITE( BINFIL, REC = 6 ) ( Z(   IZ ), IZ = 1, NTOT ) 

  ! Read in eigenvalues, k( I )

  READ( ASCFIL, * ) 
  DO MODE = 1, M 
     READ( ASCFIL, * ) kR, kI 
     k( MODE ) = CMPLX( kR, kI ) 
  ENDDO

  ! Write out eigenvalues

  IFIRST = 1 
  DO IREC = 1, 1 + ( 2 * M - 1 ) / LRECL 
     ILAST  = MIN( M, IFIRST + LRECL / 2 - 1 ) 

     WRITE( BINFIL, REC = 6 + M + IREC ) k( IFIRST : ILAST )   

     IFIRST = ILAST + 1 
  END DO

  ! Loop to read and write the modes

  DO MODE = 1, M 
     READ( ASCFIL, * ) 
     DO IZ = 1, NTOT 
        READ(  ASCFIL, * ) PHIR, PHII 
        PHI( IZ ) = CMPLX( PHIR, PHII ) 
     ENDDO
     WRITE( BINFIL, REC = 6 + MODE ) PHI( 1 : NTOT ) 
  END DO

  STOP 
END PROGRAM MODBIN
