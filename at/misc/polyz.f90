FUNCTION PZ( x0, x, F, N )

  ! Polynomial approximant of order N at x0

  IMPLICIT NONE
  INTEGER N, i, j
  COMPLEX ( KIND=8 ) PZ, x0, x( N ), f( N ), fT( N ), h( N )

  ! Initialize arrays
  h  = x - x0
  fT = f

  ! Recursion for solution
  IF ( N >= 2 ) THEN
     DO i = 1, N - 1
        DO j = 1, N - i
           fT( j ) = ( h( j + i ) * fT( j ) - h( j ) * fT( j + 1 ) ) / ( h( j + i ) - h( j ) )
        END DO
     END DO
  ENDIF
  PZ = fT( 1 )

END FUNCTION PZ
