MODULE BackSub_Mod

  ! MICHAEL B. PORTER 7/1/85

  ! Based on TINVIT in EISPACK, once upon a time ...
  ! performs back-substitution for a symmetric tridiagonal linear system
  !
  ! On input
  !    N is the order of the matrix
  !    RV1,2,3,4 contain the LU factorization of A
  !    B contains the right hand side (Ax=b)
  ! On output
  !    b contains the solution

  IMPLICIT NONE
  INTEGER, PRIVATE :: i

  INTERFACE BackSub
     MODULE PROCEDURE BackSub_sngl, BackSub_dble
  END INTERFACE BackSub

CONTAINS

  SUBROUTINE BackSub_sngl( N, Rv1, Rv2, Rv4, b )

    INTEGER, INTENT( IN    ) :: N
    COMPLEX, INTENT( IN    ) :: Rv1( N ), Rv2( N ), Rv4( N )
    COMPLEX, INTENT( INOUT ) :: b( N )

    ! Forward elimination

    DO i = 2, N
       b( i ) = b( i ) - Rv4( i ) * b( i - 1 )
       ! prevent underflows from propagating:
       !IF ( ABS(  REAL( b( i ) ) ) < SMALL .AND. ABS( AIMAG( b( i ) ) ) < SMALL ) b( i ) = 0.0
    END DO

    ! Back-substitution (result in b)
    b( N ) = b( N ) * Rv1( N )

    IF ( N >= 2 ) THEN
       DO i = N - 1, 1, -1
          b( i ) = ( b( i ) - b( i + 1 ) * Rv2( i ) ) * Rv1( i )
       END DO
    END IF

  END SUBROUTINE BackSub_sngl

  !------------------------------------------------------------------------------

  SUBROUTINE BackSub_dble( N, Rv1, Rv2, Rv4, b )

    INTEGER,          INTENT( IN    ) :: N
    COMPLEX (KIND=8), INTENT( IN    ) :: Rv1( N ), Rv2( N ), Rv4( N )
    COMPLEX (KIND=8), INTENT( INOUT ) :: b( N )

    ! Forward elimination

    DO i = 2, N
       b( i ) = b( i ) - Rv4( i ) * b( i - 1 )
       ! prevent underflows from propagating:
       !IF ( ABS(  REAL( b( i ) ) ) < SMALL .AND. ABS( AIMAG( b( i ) ) ) < SMALL ) b( i ) = 0.0
    END DO

    ! Back-substitution (result in b)
    b( N ) = b( N ) * Rv1( N )

    IF ( N >= 2 ) THEN
       DO i = N - 1, 1, -1
          b( i ) = ( b( i ) - b( i + 1 ) * Rv2( i ) ) * Rv1( i )
       END DO
    END IF

  END SUBROUTINE BackSub_dble

END MODULE BackSub_Mod
