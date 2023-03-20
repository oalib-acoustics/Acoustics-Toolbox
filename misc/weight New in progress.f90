SUBROUTINE Weight( x, Nx, xTab, NxTab, w, Ix )

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
  ! x is assumed to be monotically increasing
  ! xTab can be in any order
  !
  ! If x is not monotonically increasing then this routine will still work iff xTab = x( L ) for some L

  IMPLICIT NONE
  INTEGER, INTENT( IN  ) :: Nx, NxTab
  INTEGER                :: L, IxTab, IxTemp( 1 )
  INTEGER, INTENT( OUT ) :: Ix( NxTab )
  REAL,    INTENT( IN  ) :: x( Nx ), xTab( NxTab )
  REAL,    INTENT( OUT ) :: w( NxTab )

  ! Quick return if just one X value for interpolation ***
  IF ( Nx == 1 ) THEN
     w(  1 ) = 0.0
     Ix( 1 ) = 1
     RETURN
  ENDIF

  DO IxTab = 1, NxTab   ! Loop over each point for which the weights are needed

     ! get index, L, of the nearest point in the table
     IxTemp = MINLOC( ABS( x - xTab( IxTab ) ) )
     IF (  xTab( IxTab ) < x( IxTemp ) )
        L = IxTemp - 1

     !
     IF ( ABS( x( IxTemp - 1 ) - xTab( IxTab ) ) < ABS( x( IxTemp + 1 ) - xTab( IxTab ) ) )
        L = IxTemp - 1
     ELSE
        L = IxTemp + 1
     END IF

!!$     L = 1
!!$     DO WHILE ( xTab( IxTab ) >= x( L + 1 ) .AND. L < Nx - 1 )
!!$        L = L + 1
!!$     END DO

     ! make note of index, L, and associated weight for interpolation
     ! L = Ix( IxTab )
     Ix( IxTab ) = L
     w(  IxTab ) = ( xTab( IxTab ) - x( L ) ) / ( x( L + 1 ) - x( L ) )

     ! write( *, * ) IxTab, L, w( IxTab )
  END DO

END SUBROUTINE WEIGHT
