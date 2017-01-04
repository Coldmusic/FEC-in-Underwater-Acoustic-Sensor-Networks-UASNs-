PROGRAM COVAR 

!!! this version has not been tested since converting to using the module SRDRD
  !     Takes replica shade files and computes synthesized                
  !     data vectors with Gaussian random noise added                     
  !     Averages data vector to from covariance matrix                    

  !     Data vectors are produced for the phone positions in              
  !     the shade file (not those of the mode file).                      

  USE SdRdRMod                                              
  INTEGER, PARAMETER :: INPFIL = 5, PRTFIL = 6, SHDFIL = 20, DATFIL = 21, COVFIL = 22, &
       MaxM = 1500, MaxNR = 5000, MaxPHONES = 101
  REAL,    PARAMETER :: PI = 3.141592                                                         
  INTEGER   ISEED( 4 ), MIN_LOC( 1 )
  REAL      COVISO( MaxPHONES, MaxPHONES )
  COMPLEX   P( MaxNR ), PHIS( MaxM ), PHIR( MaxM, MaxPHONES ), k( MaxM )
  COMPLEX, ALLOCATABLE :: D( : ), DNOISY( : ), NOISE( : ), COV( :, : )
  CHARACTER TITLE*80, TITLE2*80, PLTTYP*10, NOISETYPE*1 

  SAVE ISEED 
  DATA ISEED /2301, 4562, 3219, 1823/ 

  ! *** Read data file header ***                                     

  CALL ReadHeader( SHDFIL, ' ', Title, Freq, Atten, DELTAR, PLTTYP, xs, ys, theta )

  OPEN ( FILE = 'DATFIL', UNIT = DATFIL, STATUS = 'NEW' ) 
  OPEN ( FILE = 'COVFIL', UNIT = COVFIL, STATUS = 'NEW' ) 

  ! Read user specification for source position and compute indices   

  READ( INPFIL, * ) TITLE 
  READ( INPFIL, * ) SDEPTH, SRANGE 
  READ( INPFIL, * ) NOISETYPE, SNRDB 
  READ( INPFIL, * ) NTIMES, ISEED( 4 ) 
  !     NTIMES = 30000  ! **********************************              
  !     SNRDB  = 200.0  ! **********************************             
  ! Convert km to meters          
  SRANGE = 1000.0 * SRANGE 
  SNR    = 10.0 ** ( SNRDB / 10.0 ) 

  ! *** Calculate source depth index ***                              

  MIN_LOC = MINLOC( ABS( SD(1:NSD) - SDEPTH ) )
  ISDEP = MIN_LOC( 1 )

  ! Check to see if within tolerance                              

  IF ( FREQ == 0.0 ) THEN 
     TOL = 0.001 
  ELSE
     TOL = 1500.0 / ( 6.0 * FREQ )   ! one wavelength
  END IF

  IF ( ABS( SD( ISDEP ) - SDEPTH ) > TOL ) THEN 
     WRITE( PRTFIL, * ) 'Sources exist at the following depths:' 
     WRITE( PRTFIL, * ) ( SD( I ), I = 1, NSD ) 
     CALL ERROUT( PRTFIL, 'F', 'COVAR', 'No source at specified depth' )                  
  END IF

  ! *** Calculate source range index ***                              

  MIN_LOC = MINLOC( ABS( R(1:NR) - SRANGE ) )
  ISRAN = MIN_LOC( 1 )

  ! Check to see if within tolerance                              

  IF ( ABS( R( ISRAN ) - SRANGE ) > TOL ) THEN 
     WRITE( PRTFIL, * ) 'Sources exist at the following ranges:' 
     WRITE( PRTFIL, * ) R( 1 : NR ) 
     CALL ERROUT( PRTFIL, 'F', 'COVAR', 'No source at specified range' )                  
  END IF

  ! extract data vector                                           

  WRITE( PRTFIL, * ) 'Data vector being extracted for SDEP, SRAN', SD( ISDEP), R( ISRAN )
  ALLOCATE( D( NRD ), DNOISY( NRD ) )

  IREC = ( ISDEP - 1 ) * NRD + 6 
  DO IRD1 = 1, NRD
     IREC = IREC + 1 
     READ ( SHDFIL, REC = IREC ) P( 1 : NR ) 
     D( IRD1 ) = P( ISRAN )
  END DO

  ! Normalize data vector so that array power = 1                 

  D = D / SQRT( DOT_PRODUCT( D, D ) )

  PHNPOWER = 1.0 / NRD 

  ! Zero out covariance matrix                                    

  ALLOCATE( NOISE( NRD ), COV( NRD, NRD ) )
  NOISE = 0.0
  COV   = 0.0

  ! *** If making colored noise, need modal data ***                  

  IF ( NOISETYPE == 'C' ) THEN 
     IPROF   = 1 
     ! noise sources are at 1 m                 
     SD( 1 ) = 1.0 
     NSD     = 1 

     CALL GETMOD( IPROF, ' ', MaxM, SD, NSD, 'N', k, PHIS, M, FREQ, TITLE2 )                                 
     CALL GETMOD( IPROF, ' ', MaxM, RD, NRD, 'N', k, PHIR, M, FREQ, TITLE2 )                                 

     IF ( ANY( AIMAG( k( 1:M ) ) >= 0.0 ) ) THEN 
        CALL ERROUT( PRTFIL, 'F', 'COVAR', 'Modes present with no loss: Noise field is infinite ' )   
     END IF

     ! Construct PHI and POWER
     DO MODE = 1, M 
        PHIR( MODE, 1:NRD ) = PHIS( MODE ) * PHIR( MODE, 1:NRD ) / &
             SQRT( -REAL( k( MODE ) ) * AIMAG( k( MODE ) ) )
        POWER = SUM( PHIR( MODE, 1:NRD ) * CONJG( PHIR( MODE, 1:NRD ) ) )
     END DO
  END IF

  ! header for data vector file                                   

  WRITE( DATFIL, * ) '''', TITLE(1:60), '''' 
  WRITE( DATFIL, * ) FREQ 
  WRITE( DATFIL, * ) NTIMES, NRD 

  DO IRD1 = 1, NRD
     WRITE( DATFIL, * ) RD( IRD1 )
  END DO

  ! *** Main loop: generate realizations of noise field ***           

  AVGSIGNAL = 0.0
  AVGNOISE  = 0.0
  SIGMA = SQRT( PHNPOWER / SNR )   ! Noise level on a phone

  DO ITIME = 1, NTIMES 

     IF (      NOISETYPE == 'W' ) THEN 
        CALL WHITE( ISEED, SIGMA, NOISE, NRD ) 
     ELSE IF ( NOISETYPE == 'C' ) THEN 
        CALL COLOR( ISEED, PHIR, M, NRD, SIGMA, POWER, NOISE ) 
     ELSE 
        STOP 'Unknown noise type' 
     END IF

     DNOISY = D + NOISE

     AVGSIGNAL = SUM( ABS( D     ) ** 2 )
     AVGNOISE  = SUM( ABS( NOISE ) ** 2 )

     ! write noise vector to DATFIL                               
     !         WRITE( DATFIL, 1000 ) ( ITIME, IRD, DNOISY( IRD ) , IRD = 1, NRD )               
1000 FORMAT( I3, I3, ' (', F11.7, ',', F11.7, ' )' ) 

     ! Compute cross-sensor correlation matrix R                  

     DO I = 1, NRD
        COV( I, : ) = COV( I, : ) + DNOISY( I ) * CONJG( DNOISY( : ) )
     END DO

  END DO

  AVGSIGNAL = AVGSIGNAL / ( NTIMES * NRD ) 
  AVGNOISE  = AVGNOISE  / ( NTIMES * NRD ) 

  WRITE( PRTFIL, * ) 'AVGSIGNAL = ', AVGSIGNAL 
  WRITE( PRTFIL, * ) 'AVGNOISE  = ', AVGNOISE 

  WRITE( PRTFIL, * ) 'Average signal (dB) = ', 10.0 * LOG10( AVGSIGNAL )                                      
  WRITE( PRTFIL, * ) 'Average noise  (dB) = ', 10.0 * LOG10( AVGNOISE  )                                      
  WRITE( PRTFIL, * ) 'Average SNR    (dB) = ', 10.0 * LOG10( AVGSIGNAL / AVGNOISE )                           

  ! *** Analytic formula for covariance matrix in isovelocity ocean

  COV( 1:NRD, 1:NRD ) = 0.0

  DO I = 1, NRD
     DO J = 1, NRD
        DEPTH = 100.0 
        DO MODE = 1, M 
           GAMMA = ( MODE - .5 ) * PI / DEPTH 
           COVISO( I, J ) = COVISO( I, J ) + SIN( GAMMA * SD( 1 ) ) ** 2 * &
                SIN( GAMMA * RD( I ) ) * SIN( GAMMA * RD( J ) )
        END DO
     END DO
  END DO

  ! *** Write normalized cross-sensor correlation matrix column-wise  

  WRITE( COVFIL, * ) '''', TITLE(1:60), '''' 
  WRITE( COVFIL, * ) FREQ 
  WRITE( COVFIL, * ) NRD 
  DO IRD1 = 1, NRD
     WRITE( COVFIL, * ) RD( IRD1 )
  END DO

  DO I = 1, NRD 
     DO J = 1, NRD 
        COV( I, J ) = NRD * COV( I,J ) / ( ( 1.0 + SIGMA**2 ) * NTIMES )                     
        WRITE( COVFIL, 1000 ) I, J, COV( I, J ) 
        !               ABS( COVISO( I, J ) * COV( 1, 1 ) / COVISO( 1, 1 ) ) 
     END DO
  END DO

  STOP 
END PROGRAM COVAR

!----------------------------------------------------------------------C

SUBROUTINE WHITE( ISEED, SIGMA, NOISE, NRD ) 

  ! Generates Gaussian random white noise                             

  COMPLEX i 
  PARAMETER ( i = (0.0, 1.0), PI = 3.1415926 ) 

  INTEGER   ISEED( 4 ) 
  COMPLEX   NOISE( * ) 

  DO IRD = 1, NRD
     CALL RANDOM_NUMBER( X )
     CALL RANDOM_NUMBER( Y ) 
     X = X + 0.000001
     NOISE( IRD ) = SIGMA * SQRT( -LOG( X ) ) * EXP( i * 2.0 * PI * Y )
  END DO

  ! Note power = NRD * SIGMA ** 2                                 

  RETURN 
END SUBROUTINE WHITE

!----------------------------------------------------------------------C

SUBROUTINE COLOR( ISEED, PHI, M, NRD, SIGMA, POWER, NOISE ) 

  ! Generates colored noise vector                                    

  COMPLEX i
  PARAMETER ( MaxM = 1500, i = (0.0, 1.0), PI = 3.1415926 ) 

  INTEGER   ISEED( 4 ) 
  COMPLEX PHI( MaxM, NRD ), NOISE( NRD ), COEF

  NOISE( 1:NRD ) = 0.0  ! Zero out noise vector

  ! loop over modes                                               
  DO MODE = 1, M 

     ! compute the random coefficient                             
     CALL RANDOM_NUMBER( X )
     CALL RANDOM_NUMBER( Y ) 
     X = X + 0.000001

     COEF = SQRT( -LOG( X ) ) * EXP( i * 2.0 * PI * Y ) 
     NOISE = NOISE + COEF * PHI( MODE, : )
  END DO

  ! Normalize so that power = NRD * SIGMA**2                      

  NOISE = SIGMA * NOISE * SQRT( NRD / POWER )

  RETURN 
END SUBROUTINE COLOR
