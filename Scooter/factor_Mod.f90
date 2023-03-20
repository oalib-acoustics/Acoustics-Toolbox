MODULE Factor_Mod

  ! MICHAEL B. PORTER 7/1/85
  !
  ! Based on TINVIT in EISPACK
  ! Gaussian elimination to factor a symmetric tridiagonal linear system
  !
  ! On input
  !    N is the order of the matrix
  !    d contains the diagonal elements of the input matrix
  !    e              subdiagonal
  !      in its last N-1 positions.  e(1) is arbitrary
  ! On output
  !    Rv1, Rv2, Rv3 and Rv4 contain the factored matrix

  IMPLICIT NONE
  INTEGER, PRIVATE :: i

  INTERFACE Factor
     MODULE PROCEDURE Factor_sngl, Factor_dble
  END INTERFACE Factor

CONTAINS

  SUBROUTINE Factor_sngl( N, d, e, Rv1, Rv2, Rv4 )

    INTEGER, INTENT( IN  ) :: N
    COMPLEX, INTENT( IN  ) :: d( N ), e( N )
    COMPLEX, INTENT( OUT ) :: Rv1( N ), Rv2( N ), Rv4( N )
    COMPLEX  :: U, V, XU

    ! LU decomposition without interchanges

    U = d( 1 )
    V = e( 2 )

    DO i = 2, N-1
       XU         = e( i ) / U
       Rv4( i   ) = XU
       Rv1( I-1 ) = 1.0 / U
       Rv2( I-1 ) = V
       U          = d( i ) - XU * V
       V          = e( I+1 )
    END DO

    ! following is same as above with i = N except V=e(N+1) is not used
    XU         = e( N ) / U
    Rv4( N   ) = XU
    Rv1( N-1 ) = 1.0 / U
    Rv2( N-1 ) = V
    U          = d( N ) - XU * V

    IF ( U == 0.0 ) WRITE( *, * ) 'Singular matrix'
    Rv1( N ) = 1.0 / U
    Rv2( N ) = 0.0

  END SUBROUTINE Factor_sngl

  !------------------------------------------------------------------------------
  
  Subroutine Factor_dble( N, d, e, Rv1, Rv2, Rv4 )

    INTEGER,          INTENT( IN  ) :: N
    COMPLEX (KIND=8), INTENT( IN  ) :: d( N ), e( N )
    COMPLEX (KIND=8), INTENT( OUT ) :: Rv1( N ), Rv2( N ), Rv4( N )
    COMPLEX (KIND=8) :: U, V, XU

    ! LU decomposition without interchanges

    U = d( 1 )
    V = e( 2 )

    DO i = 2, N-1
       XU         = e( i ) / U
       Rv4( i   ) = XU
       Rv1( I-1 ) = 1.0 / U
       Rv2( I-1 ) = V
       U          = d( i ) - XU * V
       V          = e( I+1 )
    END DO

    ! following is same as above with i = N except V=e(N+1) is not used
    XU         = e( N ) / U
    Rv4( N   ) = XU
    Rv1( N-1 ) = 1.0 / U
    Rv2( N-1 ) = V
    U          = d( N ) - XU * V

    IF ( U == 0.0 ) WRITE( *, * ) 'Singular matrix'
    Rv1( N ) = 1.0 / U
    Rv2( N ) = 0.0

  END SUBROUTINE Factor_dble

END MODULE Factor_Mod
