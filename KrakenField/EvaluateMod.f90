MODULE EvaluateMod
  IMPLICIT NONE

CONTAINS
  SUBROUTINE Evaluate( C, phi, Nz, R, Nr, ro, k, M, Option, P )

    ! Given modes and wavenumbers, compute pressure field on a grid Nz x Nr
    ! Normalized to pressure of point source at 1 meter
    ! ro is a vector of range offsets that describes array tilt

    ! Option:
    !    X Cartesian   (x, z) coordinates
    !    S scaled cylindrical (sqrt(r) removed)
    !    R Cylindrical (r, z) coordinates

    INTEGER, PARAMETER    :: MaxM = 20000, MinExp = -100
    REAL,    PARAMETER    :: pi = 3.1415926
    COMPLEX, PARAMETER    :: i = ( 0.0, 1.0 )
    INTEGER, INTENT( IN ) :: M, Nr, Nz
    REAL (KIND=8), INTENT( IN ) :: ro( Nz ), r( Nr )
    COMPLEX, INTENT( IN ) :: C( M ), phi( MaxM, Nz ), k( MaxM )
    COMPLEX, INTENT( OUT) :: P( Nz, Nr )
    CHARACTER (LEN=50), INTENT( IN ) :: Option
    INTEGER               :: ir, iz
    COMPLEX               :: Hank( M ), ik( M ), const( M ), Cmat( M, Nz ), factor

    ! If no modes, return vanishing pressure
    IF ( M <= 0 ) THEN
       P = 0.0
       RETURN
    END IF

    ! Initialization
    factor = i * SQRT( 2.0 * pi ) * EXP( i * pi / 4.0 )

    IF ( Option( 1 : 1 ) == 'X' ) THEN   ! line source
       const( 1 : M ) = factor * C( 1 : M ) / k( 1 : M )
    ELSE                                 ! point source or scaled cylindrical
       const( 1 : M ) = factor * C( 1 : M ) / SQRT( k( 1 : M ) )
    END IF

    ik( 1 : M ) = -i * k( 1 : M )   ! use e{ i( w t - k r ) } form
    IF ( Option( 4 : 4 ) == 'I' ) ik = REAL( ik )   ! Incoherent case

    ! phase changes due to array range-offsets representing array tilt
    DO iz = 1, Nz
       Cmat( :, iz ) = CMPLX( const( : ) * phi( 1 : M, iz ) * EXP( ik( : ) * ro( iz ) ) )
    END DO

    Ranges: DO ir = 1, Nr
       ! eliminate underflows (can raise CPU time)
       !WHERE (  REAL( ik * r( ir ) ) > MinExp )
       Hank = CMPLX( EXP( ik * r( ir ) ) )
       !ELSEWHERE
       !   Hank = 0.0
       !END WHERE

       IF ( Option( 4 : 4 ) /= 'I' )  THEN     ! coherent   case
          !!! This should be written as a simple matrix multiply
          Depths: DO iz = 1, Nz
             P( iz, ir ) =      SUM(   Cmat( :, iz ) * Hank( : ) )
          END DO Depths
       ELSE                                    ! incoherent case
          DepthsInc: DO iz = 1, Nz
             ! P( iz, ir ) = SQRT( SUM(  ABS( Cmat( :, iz ) * Hank( : ) ) ** 2 ) )
             P( iz, ir ) = SQRT( SUM( ( Cmat( :, iz ) * Hank( : ) ) ** 2 ) )
          END DO DepthsInc
       END IF

       ! Optionally include cylindrical spreading
       IF ( Option( 1 : 1 ) == 'R' ) THEN
          WHERE ( ABS( r( ir ) + ro( : ) ) > TINY( R( 1 ) ) )
             P( :, ir ) = CMPLX( P( :, ir ) / SQRT( r( ir ) + ro( : ) ) )
          END WHERE
       END IF

    END DO Ranges

  END SUBROUTINE Evaluate
END MODULE EvaluateMod
