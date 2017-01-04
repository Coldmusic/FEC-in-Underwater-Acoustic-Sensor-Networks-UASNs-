!     Last change:  MBP  10 Feb 2001   11:49 am
      SUBROUTINE XERROR( MESSG, NMESSG, NERR, LEVEL )

!     Copper penny fuse

!     The original xerror optionally prints a warning message
!     Depending on severity, previous occurence, ...

!     This is a dumb version which always prints the error message
!     Michael B. Porter 4/27/86

      CHARACTER MESSG*72
      INTEGER NMESSG, NERR, LEVEL

      WRITE( 6, '('' '',A72)' ) MESSG

      RETURN
      END
