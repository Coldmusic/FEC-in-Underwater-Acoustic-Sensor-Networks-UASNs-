SUBROUTINE EVALPDQ( k, phiR, phiS, M, IelementSC, XS, YS, theta, Ntheta, RminM, RmaxM, NR, MACT, P )

  ! Computes 3-D pressure field using adiabatic mode theory.           
  ! Normalized to pressure of point source at 1 meter.                 
  ! Note RminM must be zero.

  USE ElementMod

  INTEGER, PARAMETER :: maxM = 200
  REAL,    PARAMETER :: pi = 3.141592, DEGRAD = pi/180.0
  COMPLEX, PARAMETER :: i = ( 0.0, 1.0 )
  INTEGER            :: M( * ), outside
  REAL               :: theta( * )
  COMPLEX            :: phiR( maxM, * ), phiin( maxM ), phiout( maxM ),  phiS( maxM, * ), coef( maxM ), &
       &  k(   maxM, * ), kin( maxM ), kout( maxM ),  P( Ntheta, * )
  COMPLEX, ALLOCATABLE :: PhaseInc( : ), kAVG( : )
  
  deltar = ( RmaxM - RminM ) / ( NR - 1 )

  DO Itheta = 1, Ntheta    ! Loop over angle
     tsx = COS( DEGRAD * theta( Itheta ) )
     tsy = SIN( DEGRAD * theta( Itheta ) )

     ! Get modal values                                          
     Ielement = IelementSC
     WRITE( *, * )
     WRITE( *, * ) 'Itheta, tsx, tsy', Itheta, tsx, tsy

     CALL FIRST( Ielement, outside, Rin, Rout, XS, YS, tsx, tsy,        &
          &      Mprop, M, maxM, k, phiR, phiS, coef, kin, phiin, kout, phiout )
     Mprop = MIN( MACT, Mprop )

     IF ( Itheta == 1 ) ALLOCATE( PhaseInc( Mprop ), kAVG( Mprop ) )
     coef(     1:Mprop ) = SQRT( 2.0 * pi ) * EXP( i * pi / 4.0 ) * coef( 1:Mprop )                         
     kAVG(     1:Mprop ) = 0.5 * ( kin( 1:Mprop ) + kout( 1:Mprop ) ) 
     PhaseInc( 1:Mprop ) = EXP( -i * kAVG( 1:Mprop ) * deltar )

     DO IR = 1, NR ! March forward in range 
        RM = RminM + ( IR - 1 ) * deltar 
        IF ( RM == 0.0 ) RM = deltar 

        ! Crossing into new element?                                  
        DO WHILE ( RM > Rout )

           ! Copy outside info to inside                          
           NEWelement = ADJELT( outside, Ielement ) 
           Rin = Rout
           kin(      1:Mprop ) =  kout(  1:Mprop ) 
           phiin(    1:Mprop ) = phiout( 1:Mprop ) 
           kAVG(     1:Mprop ) = 0.5 * ( kin( 1:Mprop ) + kout( 1:Mprop ) ) 
           PhaseInc( 1:Mprop ) = EXP( -i * kAVG( 1:Mprop ) * deltar )

           ! Get new outside info                                 
           CALL OUT( Ielement, NEWelement, outside, Rout, XS, YS, tsx, tsy, Mprop, M, maxM, k, phiR, kout, phiout )
           Ielement = NEWelement 

        END DO   ! Check again

        ! Compute modal contribution at this range

        coef( 1:Mprop ) = coef( 1:Mprop ) * PhaseInc( 1:Mprop )
        P( Itheta, IR ) = SUM( coef( 1:Mprop ) * phiin( 1:Mprop ) ) / SQRT( RM * kin( 1 ) ) 

     END DO ! next range
  END DO ! next bearing

END SUBROUTINE EVALPDQ
