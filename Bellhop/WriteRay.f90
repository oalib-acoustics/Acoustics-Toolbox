MODULE WriteRay

  ! Compress the ray data keeping every iSkip point, points near surface or bottom, and last point.
  ! (The subsampling is currently disabled by setting iSkip=1, but the logic is left in for optional use.)
  ! Write to RAYFile.

  ! During an eigenray calculation, subsets of the full ray may be passed
  ! These have lengths Nsteps1 vs. Nsteps for the entire ray

  USE BellhopMod
  USE sspMod
  IMPLICIT NONE
  INTEGER, PRIVATE :: MaxNRayPoints = 500000   ! this is the maximum length of the ray vector that is written out
  INTEGER, PRIVATE :: is, N2, iSkip

CONTAINS

  SUBROUTINE WriteRay2D( alpha0, Nsteps1 )

    ! The 2D version is for ray traces in (r,z) coordinates

    INTEGER,       INTENT( IN ) :: Nsteps1
    REAL (KIND=8), INTENT( IN ) :: alpha0   ! take-off angle of this ray

    ! compression

    N2    = 1
    iSkip = MAX( Nsteps1 / MaxNRayPoints, 1 )

    Stepping: DO is = 2, Nsteps1
       ! ensure that we always write ray points near bdry reflections (works only for flat bdry)
       IF ( MIN( Bdry%Bot%HS%Depth - ray2D( is )%x( 2 ),  ray2D( is )%x( 2 ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
            MOD( is, iSkip ) == 0 .OR. is == Nsteps1 ) THEN
          N2 = N2 + 1
          ray2D( N2 )%x = ray2D( is )%x
       END IF
    END DO Stepping

    ! write to ray file

    WRITE( RAYFile, * ) alpha0
    WRITE( RAYFile, * ) N2, ray2D( Nsteps1 )%NumTopBnc, ray2D( Nsteps1 )%NumBotBnc

    DO is = 1, N2
       WRITE( RAYFile, * ) ray2D( is )%x
    END DO

  END SUBROUTINE WriteRay2D

  ! **********************************************************************!

  SUBROUTINE WriteRay3D( alpha0, beta0, Nsteps1 )

    ! The 3D version is for ray traces in (x,y,z) coordinates

    INTEGER,       INTENT( IN ) :: Nsteps1
    REAL (KIND=8), INTENT( IN ) :: alpha0, beta0   ! take-off angle of this ray

    ! if Nx2D run, copy r-z rays to x-y-z rays

    IF ( Beam%RunType( 6 : 6 ) == '2' ) THEN
       ray3D%x( 1 )    = xs_3D( 1 ) + ray2D%x( 1 ) * COS( beta0 )
       ray3D%x( 2 )    = xs_3D( 2 ) + ray2D%x( 1 ) * SIN( beta0 )
       ray3D%x( 3 )    = ray2D%x( 2 )
       ray3D%NumTopBnc = ray2D%NumTopBnc
       ray3D%NumBotBnc = ray2D%NumBotBnc
    END IF

    ! compression

    N2    = 1
    iSkip = MAX( Nsteps1 / MaxNRayPoints, 1 )

    Stepping: DO is = 2, Nsteps1
       ! ensure that we always write ray points near boundary reflections
       IF ( MIN( Bdry%Bot%HS%Depth - ray3D( is )%x( 3 ),  ray3D( is )%x( 3 ) - Bdry%Top%HS%Depth ) < 0.2 .OR. &
            MOD( is, iSkip ) == 0 .OR. is == Nsteps1 ) THEN
          N2 = N2 + 1
          ray3D( N2 )%x = ray3D( is )%x
       END IF
    END DO Stepping

    ! write to ray file

    WRITE( RAYFile, * ) alpha0
    WRITE( RAYFile, * ) N2, ray3D( Nsteps1 )%NumTopBnc, ray3D( Nsteps1 )%NumBotBnc

    DO is = 1, N2
       WRITE( RAYFile, * ) ray3D( is )%x
    END DO

  END SUBROUTINE WriteRay3D

END MODULE WriteRay
