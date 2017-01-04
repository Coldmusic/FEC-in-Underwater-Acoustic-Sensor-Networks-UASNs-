MODULE spamod

  SAVE
  INTEGER,       PARAMETER :: ENVFile = 5, PRTFile = 6, GRNFile = 25, RTSFile = 35, MaxN = 17000, MaxMedium = 500, MaxIt = 1000000
  COMPLEX,       PARAMETER :: i = ( 0.0, 1.0 )
  REAL (KIND=8), PARAMETER :: PI = 3.1415926535898D0
  INTEGER              :: LOC( MaxMedium ), N( MaxMedium ), NMedia, NSig, Nk, NTot1, Nrr, NTout
  REAL                 :: C2R( MaxN ), C2I( MaxN ), Z( MaxN ), &
                          H( MaxMedium ), CMin, CLow, CHigh, Omega2, DeltaK, &
                          Deltat, CrossT, CMax, TStart, V, TMult, alpha, beta, FMin, FMax
  REAL        (KIND=8) :: Depth( MaxMedium ), Rho( MaxN ), sigma( MaxMedium )
  CHARACTER (LEN=8 )   :: Material( MaxMedium ), TopOpt, BotOpt
  CHARACTER (LEN=4 )   :: Pulse
  CHARACTER (LEN=80)   :: Title
  REAL,    ALLOCATABLE :: k( : ), Tout( : ), RTSrd( :, : ), RTSrr( :, : )
  COMPLEX, ALLOCATABLE :: Green( :, :, : )

END MODULE spamod
