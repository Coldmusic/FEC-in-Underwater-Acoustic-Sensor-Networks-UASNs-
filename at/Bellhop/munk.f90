SUBROUTINE ANALYT( x, c, gradc, crr, crz, czz )

  IMPLICIT NONE
  REAL, INTENT( IN  ) :: x( 2 )
  REAL, INTENT( OUT ) :: c, gradc( 2 ), crr, crz, czz
  REAL                :: c0, cr, cz, DxtDz, xt

  c0 = 1500.0

  IF ( x( 2 ) < 5000.0 ) THEN
     xt    = 2.0 * ( x( 2 ) - 1300.0 ) / 1300.0
     DxtDz = 2.0 / 1300.0
     c     = C0 * ( 1.0 + 0.00737*( xt - 1.0 + EXP( -xt ) ) )
     cz    = C0 * 0.00737 * ( 1.0 - EXP( -xt ) ) * DxtDz
     czz   = C0 * 0.00737 * EXP( -xt ) * DxtDz ** 2
  ELSE
     ! Homogeneous half-space
     xt   = 2.0 * ( 5000.0 - 1300.0 ) / 1300.0
     c    = C0 * ( 1.0 + 0.00737 * ( xt - 1.0 + EXP( -xt ) ) )
     cz   = 0.0
     czz  = 0.0
  ENDIF

  cr = 0.0
  gradc = [ cr, cz ]
  crz = 0.0
  crr = 0.0

  RETURN
END SUBROUTINE ANALYT
