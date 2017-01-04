MODULE krakcmod

  SAVE

  INTEGER,          PARAMETER   :: ENVFile = 5, PRTFile = 6, MODFile = 20, EVMFile = 22, MaxM = 20000, MaxMedium = 500, NSets = 5
  REAL    (KIND=8), PARAMETER   :: pi = 3.1415926535898D0, DegRad = pi / 180.0D0
  COMPLEX (KIND=8), PARAMETER   :: i = ( 0.0D0, 1.0D0 )

  INTEGER                       :: FirstAcoustic, LastAcoustic, NMedia, NV( NSets ), ISet, M, &
                                   LRecordLength, IRecProfile = 1, ModeCount, Mode, IProf
  REAL      (KIND=8)            :: ET( NSets ), hV( NSets ), VG( MaxM ), cMin, cLow, cHigh, Freq, omega2, RMax
  COMPLEX   (KIND=8)            :: k( MaxM ), EVMat( NSets, MaxM ), Extrap( NSets, MaxM )
  ! Halfspace properties
  TYPE HSInfo
     CHARACTER (LEN=1)          :: BC                          ! Boundary condition type
     COMPLEX (KIND=8)           :: cP, cS                      ! P-wave, S-wave speeds
     REAL    (KIND=8)           :: rho, BumpDensity, eta, xi   ! density, boss parameters
  END TYPE HSInfo

  TYPE( HSInfo )                :: HSTop, HSBot

  ! media properties
  INTEGER                       :: Loc(   MaxMedium ), NG( MaxMedium ), N(     MaxMedium )
  REAL      (KIND=8)            :: Depth( MaxMedium ), H(  MaxMedium ), sigma( MaxMedium )
  CHARACTER (LEN= 8)            :: Material( MaxMedium ), TopOpt, BotOpt

  ! storage for finite-difference equations

  CHARACTER (LEN=80)            :: Title
  REAL    (KIND=8), ALLOCATABLE :: rho( : )
  COMPLEX (KIND=8), ALLOCATABLE :: B1( : ), B2( : ), B3( : ), B4( : )

END MODULE krakcmod
