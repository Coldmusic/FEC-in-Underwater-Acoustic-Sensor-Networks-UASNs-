MODULE scomod

   SAVE
   INTEGER,          PARAMETER :: ENVFile = 5, PRTFile = 6, GRNFile = 25, IRCFile = 12, MaxMedium = 500
   REAL    (KIND=8), PARAMETER :: pi = 3.1415926535898D0
   COMPLEX (KIND=8), PARAMETER :: i  = ( 0.0, 1.0 )

   INTEGER           :: N( MaxMedium ), Loc( MaxMedium ), FirstAcoustic, LastAcoustic, NMedia, Nk
   REAL     (KIND=8) :: Depth( MaxMedium ), H( MaxMedium ), sigma( MaxMedium ), Cmin, Clow, CHigh, omega2, &
                        RMax, Atten
   CHARACTER (LEN=8) :: Material( MaxMedium ), TopOpt, BotOpt

  ! Halfspace properties
  TYPE HSInfo
     CHARACTER (LEN=1)         :: BC                          ! Boundary condition type
     COMPLEX (KIND=8)          :: cP, cS                      ! P-wave, S-wave speeds
     REAL    (KIND=8)          :: rho, BumpDensity, eta, xi   ! density, boss parameters
  END TYPE HSInfo

  TYPE( HSInfo )               :: HSTop, HSBot

END MODULE scomod
