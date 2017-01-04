MODULE sspmod

   ! holds SSP input by user and associated variables

   IMPLICIT NONE
   SAVE
   INTEGER, PARAMETER :: MaxSSP = 2001
   INTEGER            :: NSSP, Layer
   REAL  ( KIND = 8 ) :: zSSPV( MaxSSP ), cSSPV( MaxSSP ), cZV( MaxSSP ), n2V( MaxSSP ), n2ZV( MaxSSP ), &
                         cVS( 4, MaxSSP )

END MODULE sspmod
