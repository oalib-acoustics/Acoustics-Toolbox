MODULE SubTabulate

  ! If x( 3 ) = -999.9 then subtabulation is performed
  ! i.e., a vector is generated with Nx points in [ x( 1 ), x( 2 ) ]
  ! If x( 2 ) = -999.9 then x( 1 ) is repeated into x( 2 )

  ! mbp 1/2015

  IMPLICIT NONE
  INTEGER, PRIVATE :: ix

  INTERFACE SubTab
     MODULE PROCEDURE SubTab_sngl, SubTab_dble
  END INTERFACE SubTab

CONTAINS

  SUBROUTINE SubTab_sngl( x, Nx )

    INTEGER, INTENT( IN )    :: Nx
    REAL,    INTENT( INOUT ) :: x( Nx )
    REAL                     :: deltax

    IF ( Nx >= 3 ) THEN
       IF ( x( 3 ) == -999.9 ) THEN   ! testing for equality here is dangerous
          IF ( x( 2 ) == -999.9 ) x( 2 ) = x( 1 )
          deltax      = ( x( 2 ) - x( 1 ) ) / ( Nx - 1 )
          x( 1 : Nx ) = x( 1 ) + [ ( ix, ix = 0, Nx - 1 ) ] * deltax
       END IF
    END IF

  END SUBROUTINE SubTab_sngl

  SUBROUTINE SubTab_dble( x, Nx )

    INTEGER,       INTENT( IN )    :: Nx
    REAL (KIND=8), INTENT( INOUT ) :: x( Nx )
    REAL (KIND=8)                  :: deltax

    IF ( Nx >= 3 ) THEN
       IF ( x( 3 ) == -999.9 ) THEN   ! testing for equality here is dangerous
          IF ( x( 2 ) == -999.9 ) x( 2 ) = x( 1 )
          deltax      = ( x( 2 ) - x( 1 ) ) / ( Nx - 1 )
          x( 1 : Nx ) = x( 1 ) + [ ( ix, ix = 0, Nx - 1 ) ] * deltax
       END IF
    END IF

  END SUBROUTINE SubTab_dble

END MODULE SubTabulate
