MODULE RayNormals

   IMPLICIT NONE
   REAL (KIND=8)                :: RL                 ! length of part of the tangent vector

CONTAINS

  SUBROUTINE RayNormal( t, phi, c, e1, e2 )

    ! Computes the ray normals

    REAL (KIND=8), INTENT( IN  ) :: t( 3 )             ! tangent vector (NOT normalized)
    REAL (KIND=8), INTENT( IN  ) :: phi                ! torsion
    REAL (KIND=8), INTENT( IN  ) :: c                  ! sound speed
    REAL (KIND=8), INTENT( OUT ) :: e1( 3 ), e2( 3 )   ! ray unit normals

    RL = NORM2( t( 1 : 2 ) )

    IF ( phi /= 0 ) THEN
       !  e1
       e1( 1 ) = ( c * t( 1 ) * t( 3 ) * COS( phi ) + t( 2 ) * SIN( phi ) ) / RL
       e1( 2 ) = ( c * t( 2 ) * t( 3 ) * COS( phi ) - t( 1 ) * SIN( phi ) ) / RL
       e1( 3 ) = -c * RL * COS( phi )

       !  e2
       e2( 1 ) = ( c * t( 1 ) * t( 3 ) * SIN( phi ) - t( 2 ) * COS( phi ) ) / RL
       e2( 2 ) = ( c * t( 2 ) * t( 3 ) * SIN( phi ) + t( 1 ) * COS( phi ) ) / RL
       e2( 3 ) = -c * RL * SIN( phi )
    ELSE
       e1 = [ c * t( 1 ) * t( 3  ) / RL, c * t( 2 ) * t( 3 ) / RL, -c * RL ]
       e2 = [ -t( 2 ) / RL,              t( 1 ) / RL,              0.0D0   ]
    END IF

    RETURN
  END SUBROUTINE RayNormal

  !**********************************************************************!

  SUBROUTINE RayNormal_unit( t, phi, e1, e2 )

    ! Computes the ray normals
    ! Same as routine RayNormal except this version assumes t is already normalized

    REAL (KIND=8), INTENT( IN  ) :: t( 3 )             ! tangent vector (normalized)
    REAL (KIND=8), INTENT( IN  ) :: phi                ! torsion
    REAL (KIND=8), INTENT( OUT ) :: e1( 3 ), e2( 3 )   ! ray normals

    RL = NORM2( t( 1 : 2 ) )

    !  e1
    e1( 1 ) = ( t( 1 ) * t( 3 ) * COS( phi ) + t( 2 ) * SIN( phi ) ) / RL
    e1( 2 ) = ( t( 2 ) * t( 3 ) * COS( phi ) - t( 1 ) * SIN( phi ) ) / RL
    e1( 3 ) = -RL * COS( phi )

    !  e2
    e2( 1 ) = ( t( 1 ) * t( 3 ) * SIN( phi ) - t( 2 ) * COS( phi ) ) / RL
    e2( 2 ) = ( t( 2 ) * t( 3 ) * SIN( phi ) + t( 1 ) * COS( phi ) ) / RL
    e2( 3 ) = -RL * SIN( phi )

    RETURN
  END SUBROUTINE RayNormal_unit

END MODULE RayNormals

