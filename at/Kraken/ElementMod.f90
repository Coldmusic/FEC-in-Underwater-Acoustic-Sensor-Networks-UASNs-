MODULE ElementMod

  SAVE

  INTEGER,              PARAMETER :: PRTFile = 6
  INTEGER                         :: NElts, NNodes
  INTEGER,            ALLOCATABLE :: Node( :, : ), AdjElt( :, : ), Iset( : )
  REAL,               ALLOCATABLE :: x( : ), y( : )                   
  CHARACTER (LEN=50), ALLOCATABLE :: ModeFileName( : )

END MODULE ElementMod
