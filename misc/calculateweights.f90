MODULE CalculateWeights

  ! mbp 12/2018, based on much older subroutine
  ! Given 
  !    x(*)    abscissas
  !    xTab(*) points for tabulation
  !    Nx      number of x    points
  !    NxTab   number of xTab points

  ! Compute
  !    w(*)    weights for linear interpolation
  !    Ix(*)   indices for    "         "
  !
  ! If xTab is outside the domain of x, assumes extrapolation will be done

  IMPLICIT NONE

  INTERFACE Weight
     MODULE PROCEDURE Weight_sngl, Weight_dble
  END INTERFACE Weight

CONTAINS
  SUBROUTINE Weight_sngl( x, Nx, xTab, NxTab, w, Ix )

    INTEGER, INTENT( IN  ) :: Nx, NxTab
    INTEGER, INTENT( OUT ) :: Ix( NxTab )
    REAL,    INTENT( IN  ) :: x( Nx ), xTab( NxTab )
    REAL,    INTENT( OUT ) :: w( NxTab )
    INTEGER                :: L, IxTab

    ! Quick return if just one X value for interpolation
    IF ( Nx == 1 ) THEN
       w(  1 ) = 0.0
       Ix( 1 ) = 1
       RETURN
    ENDIF

    L = 1

    DO IxTab = 1, NxTab   ! Loop over each point for which the weights are needed

       ! search for index, L, such that [ x( L ), x( L+1 ) ] brackets rcvr depth
       DO WHILE ( xTab( IxTab ) > x( L + 1 ) .AND. L < Nx - 1 )
          L = L + 1
       END DO

       ! make note of index, L, and associated weight for interpolation
       Ix( IxTab ) = L
       w(  IxTab ) = ( xTab( IxTab ) - x( L ) ) / ( x( L + 1 ) - x( L ) )

       ! special code for case of a bottom following receiver
!!$     IF ( w( IxTab ) /= 0.0 ) THEN
!!$        Ix( IxTab ) = MIN( 2, NxTab )   ! use second depth in the mode file, or first depth if only one
!!$        w( IxTab )  = 0.0               ! assert that it's a perfect match to the requested depth
!!$     END IF
    END DO

  END SUBROUTINE Weight_sngl

  SUBROUTINE Weight_dble( x, Nx, xTab, NxTab, w, Ix )

    INTEGER,       INTENT( IN  ) :: Nx, NxTab
    INTEGER,       INTENT( OUT ) :: Ix( NxTab )
    REAL (KIND=8), INTENT( IN  ) :: x( Nx ), xTab( NxTab )
    REAL (KIND=8), INTENT( OUT ) :: w( NxTab )
    INTEGER                      :: L, IxTab

    ! Quick return if just one X value for interpolation
    IF ( Nx == 1 ) THEN
       w(  1 ) = 0.0
       Ix( 1 ) = 1
       RETURN
    ENDIF

    L = 1

    DO IxTab = 1, NxTab   ! Loop over each point for which the weights are needed

       ! search for index, L, such that [ x( L ), x( L+1 ) ] brackets rcvr depth
       DO WHILE ( xTab( IxTab ) > x( L + 1 ) .AND. L < Nx - 1 )
          L = L + 1
       END DO

       ! make note of index, L, and associated weight for interpolation
       Ix( IxTab ) = L
       w(  IxTab ) = ( xTab( IxTab ) - x( L ) ) / ( x( L + 1 ) - x( L ) )

       ! special code for case of a bottom following receiver
!!$     IF ( w( IxTab ) /= 0.0 ) THEN
!!$        Ix( IxTab ) = MIN( 2, NxTab )   ! use second depth in the mode file, or first depth if only one
!!$        w( IxTab )  = 0.0               ! assert that it's a perfect match to the requested depth
!!$     END IF
    END DO

  END SUBROUTINE Weight_dble
END MODULE CalculateWeights
