MODULE ElementMod

  IMPLICIT NONE
  SAVE
  INTEGER,              PARAMETER :: PRTFile = 6
  INTEGER                         :: NElts, NNodes, iCorner( 3, 2 ), i1, j1
  INTEGER,            ALLOCATABLE :: Node( :, : ), AdjacentElement( :, : ), Iset( : )
  REAL,               ALLOCATABLE :: x( : ), y( : )                   
  CHARACTER (LEN=80), ALLOCATABLE :: ModeFileName( : )
  DATA ( ( iCorner( i1, j1 ), i1 = 1, 3), j1 = 1, 2 ) /1, 2, 3, 2, 3, 1/ 
  ! ICor maps a side (1, 2 OR 3) and a local node (1 OR 2) to a       
  !      corner (1, 2, or 3) of the triangle                        

END MODULE ElementMod
