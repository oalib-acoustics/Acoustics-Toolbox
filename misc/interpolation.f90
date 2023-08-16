MODULE interpolation

CONTAINS
  SUBROUTINE interp1( x, y, xi, yi )

    ! Given a set of points (x, y), computes the interpolated value at xi
    ! using piecewise linear interpolation

    ! Assumes x is monotonically increasing

    IMPLICIT NONE
    REAL (KIND=8), INTENT( IN  ) :: x( : ), y( : ), xi( : )
    REAL (KIND=8), INTENT( OUT ) :: yi( : )
    INTEGER                      :: N, Ni, I, iseg
    REAL (KIND=8)                :: R

    N    = SIZE( x  )
    Ni   = SIZE( xi )
    iseg = 1

    ! loop over the interpolation points
    DO I = 1, Ni
       ! search for the bracketing pair of tabulated values
       DO WHILE ( xi( I ) > x( iseg + 1 ) )  ! is the xi point to the right of the current segment?
          IF ( iseg < N - 2 ) THEN
             iseg = iseg + 1
          END IF
       END DO

       DO WHILE ( xi( I ) < x( iseg     ) )  ! is the xi point to the left  of the current segment?
          IF ( iseg > 1 ) THEN
             iseg = iseg - 1
          END IF
       END DO

       ! proportional distance between points
       R       = ( xi( I ) - x( iseg ) ) / ( x( iseg + 1 ) - x( iseg ) )
       yi( I ) = ( 1.0 - R ) * y( iseg ) + R * y( iseg + 1 )
    END DO

  END SUBROUTINE interp1
END MODULE interpolation
