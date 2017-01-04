MODULE bellmod

   INTEGER,          PARAMETER :: PRTFile = 6, MaxN = 100000
   REAL    (KIND=8), PARAMETER :: pi = 3.1415926535898D0, RadDeg = 180.0D0 / pi
   COMPLEX (KIND=8), PARAMETER :: i  = ( 0.0D0, 1.0D0 )

   ! Reduce MaxN (= max # of steps along a ray) to reduce storage
   ! Note space is wasted in NumTopBnc, NumBotBnc ...

   INTEGER            :: NumTopBnc( MaxN), NumBotBnc( MaxN )
   INTEGER            :: NSteps, Nrd_per_range
   REAL    ( KIND=8 ) :: omega
   CHARACTER (LEN=6 ) :: TopOpt

   ! Halfspace properties
   TYPE HSInfo
      CHARACTER (LEN=1) :: BC       ! Boundary condition type
      COMPLEX (KIND=8)  :: cP, cS   ! P-wave, S-wave speeds
      REAL    (KIND=8)  :: rho, Depth  ! density, depth
      REAL    (KIND=8)  :: BumpDensity, eta, xi   ! Twersky boss parameters
   END TYPE

   TYPE( HSInfo )       :: HSTop, HSBot

   TYPE rayPt
       REAL (KIND=8 ) :: x( 2 ), t( 2), p( 2 ), q( 2 ), tau, c, Amp, Phase
   END TYPE

   TYPE( rayPt )      :: ray( MaxN )

END MODULE bellmod
