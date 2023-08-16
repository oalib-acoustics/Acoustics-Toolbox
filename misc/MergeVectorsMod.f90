MODULE MergeVectorsMod

  ! Merge two vectors x, y into z (duplicates are removed)
  ! x and y are assumed to be monotonically increasing

  IMPLICIT NONE

  INTEGER :: ix, iy, iz

  INTERFACE MergeVectors
     MODULE PROCEDURE MergeVectors_sngl, MergeVectors_dble
  END INTERFACE MergeVectors

CONTAINS
  SUBROUTINE MergeVectors_sngl( x, Nx, y, Ny, z, Nz )
    INTEGER, INTENT( IN  ) :: Nx, Ny
    INTEGER, INTENT( OUT ) :: Nz
    REAL,    INTENT( IN  ) :: x( Nx ), y( Ny )
    REAL,    INTENT( OUT ) :: z( * )

    ix = 1
    iy = 1
    iz = 1

    ! Loop to successively take values from vector x or y

    DO WHILE ( ix <= Nx .OR. iy <= Ny )

       IF      ( iy > Ny ) THEN             ! No more in vector y, take from x 
          z( iz ) = x( ix )
          ix      = ix + 1
          iz      = iz + 1
       ELSE IF ( ix > Nx ) THEN             ! No more in vector x, take from y
          z( iz ) = y( iy )
          iy      = iy + 1
          iz      = iz + 1
       ELSE IF ( x( ix ) <= y( iy ) ) THEN  ! Vector x has the next largest
          z( iz ) = x( ix )
          ix      = ix + 1
          iz      = iz + 1
       ELSE                                 ! Vector y has the next largest
          z( iz ) = y( iy )
          iz      = iz + 1
          iy      = iy + 1
       END IF

       ! Check that we didn't just add a duplicate
       IF ( iz > 2 ) THEN
          ! following is a test for approximately equal values; not a very good test
          IF ( ABS( z( iz - 1 ) - z( iz - 2 ) ) < 100. * EPSILON( z ) ) THEN
             iz = iz - 1
          END IF
       END IF
    END DO

    Nz = iz - 1

  END SUBROUTINE MergeVectors_sngl

    SUBROUTINE MergeVectors_dble( x, Nx, y, Ny, z, Nz )
    INTEGER,       INTENT( IN  ) :: Nx, Ny
    INTEGER,       INTENT( OUT ) :: Nz
    REAL (KIND=8), INTENT( IN  ) :: x( Nx ), y( Ny )
    REAL (KIND=8), INTENT( OUT ) :: z( * )

    ix = 1
    iy = 1
    iz = 1

    ! Loop to successively take values from vector x or y

    DO WHILE ( ix <= Nx .OR. iy <= Ny )

       IF      ( iy > Ny ) THEN             ! No more in vector y, take from x 
          z( iz ) = x( ix )
          ix      = ix + 1
          iz      = iz + 1
       ELSE IF ( ix > Nx ) THEN             ! No more in vector x, take from y
          z( iz ) = y( iy )
          iy      = iy + 1
          iz      = iz + 1
       ELSE IF ( x( ix ) <= y( iy ) ) THEN  ! Vector x has the next largest
          z( iz ) = x( ix )
          ix      = ix + 1
          iz      = iz + 1
       ELSE                                 ! Vector y has the next largest
          z( iz ) = y( iy )
          iz      = iz + 1
          iy      = iy + 1
       END IF

       ! Check that we didn't just add a duplicate
       IF ( iz > 2 ) THEN
          ! following is a test for approximately equal values; not a very good test
          IF ( ABS( z( iz - 1 ) - z( iz - 2 ) ) < 100. * EPSILON( z ) ) THEN
             iz = iz - 1
          END IF
       END IF
    END DO

    Nz = iz - 1

  END SUBROUTINE MergeVectors_dble

END MODULE MergeVectorsMod
