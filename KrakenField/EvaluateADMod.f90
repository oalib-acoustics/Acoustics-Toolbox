MODULE EvaluateADMod

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: MaxM = 20000, MaxNfreq = 1000
  INTEGER, PRIVATE            :: Nfreq
  INTEGER, PRIVATE            :: ir, iProf
  REAL      (KIND=8)          :: rLeft, rMid, wMid
  REAL      (KIND=8), PRIVATE :: freqVec( MaxNfreq )
  CHARACTER (LEN=80), PRIVATE :: Title

CONTAINS
  SUBROUTINE EvaluateAD( FileRoot, rProf, NProf, ifreq, phiS, phiR, rd, Nrd, r, Nr, M, Option, P )

    ! Computes pressure field using adiabatic mode theory
    ! Normalized to pressure at 1 meter

    ! Option:
    !   X     Cartesian   (x, z) coordinates
    !   T     Translationally invariant ocean
    !   R     Cylindrical (r, z) coordinates
    !   S     Scaled cylindrical coordinates ( 1/r fall-off removed )

    USE ReadModes
    REAL,       PARAMETER       :: pi = 3.1415926
    COMPLEX,    PARAMETER       :: i = ( 0.0, 1.0 )
    INTEGER,       INTENT( IN    ) :: Nr, Nrd, NProf, ifreq     ! number of receiver ranges, depths, profiles
    REAL,          INTENT( IN    ) :: rd( Nrd )                 ! receiver depths
    REAL (KIND=8), INTENT( IN    ) :: r( Nr )                   ! receiver ranges
    INTEGER,       INTENT( INOUT ) :: M                         ! number of modes
    COMPLEX,       INTENT( IN    ) :: phiS( M )                 ! mode shapes at source
    COMPLEX,       INTENT( INOUT ) :: phiR( MaxM, Nrd )         ! mode shapes at receivers
    REAL (KIND=8), INTENT( INOUT ) :: rProf( NProf + 1 )        ! ranges of profiles
    COMPLEX,       INTENT( OUT   ) :: P( Nrd, Nr )              ! pressure field
    CHARACTER (LEN=50), INTENT( IN ) :: Option                  ! Cartesian or cylindrical coordinates
    INTEGER                  :: ird, M1
    REAL  (KIND=8)           :: w
    COMPLEX                  :: phiL( MaxM, Nrd ), SUM, kInt( MaxM ), phiINT( MaxM ), kL( MaxM ), kR( MaxM ), kMid( MaxM )
    COMPLEX      (KIND=8)    :: Hank( M ), const( M ), sumk( MaxM ), sumkinv( MaxM )
    CHARACTER    (LEN=80)    :: FileRoot

    ! Initialization
    rProf( NProf + 1 ) = 1e20   ! HUGE( rProf( NProf ) ), but huge produces an underflow later
    iProf              = 1

    ! Receiver depths at left  of segment
    CALL GetModes( FileRoot, iProf    , ifreq, MaxM, rd, Nrd, 'N', kL, phiL, M1, freqVec, Nfreq, Title )
    M = MIN( M, M1 )

    ! Receiver depths at right of segment
    CALL GetModes( FileRoot, iProf + 1, ifreq, MaxM, rd, Nrd, 'N', kR, phiR, M1, freqVec, Nfreq, Title )
    M = MIN( M, M1 )

    const(   1 : M ) = i * SQRT( 2.0 * pi ) * EXP( i * pi / 4.0 ) * phiS( 1 : M )
    sumk(    1 : M ) = 0.0
    sumkinv( 1 : M ) = 0.0
    IF ( Option( 1 : 1 ) == 'T' ) const( 1 : M ) = const( 1 : M ) / SQRT( kL( 1 : M ) )   ! translationally invariant

    ! March forward in range
    Ranges: DO ir = 1, Nr

       ! Crossing into new range segment?
       DO WHILE ( r( ir ) > 1000.0 * rProf( iProf + 1 ) )
          CALL NewSegment( FileRoot, r, rProf, ifreq, NProf, kL, kR, phiL, phiR, M, rd, Nrd, sumk, sumkinv )
       END DO

       ! Compute proportional distance, W, and interpolate
       IF ( ir > 1 ) THEN
          rLeft = MAX( r( ir - 1 ), 1000.0 * rProf( iProf ) )
       ELSE
          rLeft = 1000.0 * rProf( iProf )
       ENDIF

       rMid = 0.5 * ( r( ir ) + rLeft )
       w    = ( r( ir ) / 1000.0 - rProf( iProf ) ) / ( rProf( iProf + 1 ) - rProf( iProf ) )
       wMid = ( rMid    / 1000.0 - rProf( iProf ) ) / ( rProf( iProf + 1 ) - rProf( iProf ) )

       kInt(    1 : M ) = CMPLX( kL(      1 : M ) + W    * ( kR( 1 : M ) - kL( 1 : M ) ) )
       kMid(    1 : M ) = CMPLX( kL(      1 : M ) + wMid * ( kR( 1 : M ) - kL( 1 : M ) ) )
       sumk(    1 : M ) = sumk(    1 : M ) + kMid( 1 : M ) * ( r( ir ) - rLeft )
       sumkinv( 1 : M ) = sumkinv( 1 : M ) + ( r( ir ) - rLeft ) / kMid( 1 : M )

       IF ( Option( 4 : 4 ) /= 'I' ) THEN   ! coherent   case
          Hank( 1 : M ) = const( 1 : M ) * CMPLX( EXP(       -i * sumk( 1 : M ) ) )
       ELSE                                 ! incoherent case
          Hank( 1 : M ) = const( 1 : M ) * CMPLX( EXP( REAL( -i * sumk( 1 : M ) ) ), KIND=8 )
       ENDIF

       SELECT CASE( Option( 1 : 1 ) )
       CASE ( 'R' )                             ! Cylindrical coordinates
          IF ( r( ir ) == 0.0 ) THEN
             Hank( 1 : M ) = 0.0
          ELSE
             Hank( 1 : M ) = Hank( 1 : M ) / SQRT( kInt( 1 : M ) * r( ir ) )
          ENDIF
       CASE ( 'X' )                             ! Cartesian coordinates
          Hank( 1 : M ) = Hank( 1 : M ) / kInt ( 1 : M )
       CASE ( 'T' )                             ! Translationally invariant
          Hank( 1 : M ) = Hank( 1 : M ) / CMPLX( SQRT( kInt( 1 : M ) * sumkinv( 1 : M ) ) )
       CASE DEFAULT                             ! Scaled cylindrical coordinates
          Hank( 1 : M ) = Hank( 1 : M ) / SQRT( kInt( 1 : M ) )
       END SELECT

       ! For each receiver, add up modal contributions
       Depths: DO ird = 1, Nrd
          phiINT( 1 : M ) = CMPLX( phiL( 1 : M, ird ) + W * ( phiR( 1 : M, ird ) - phiL( 1 : M, ird ) ) )
          IF ( Option( 4 : 4 ) /= 'I' )  THEN   ! coherent   case
             P( ird, ir ) =       CMPLX( SUM(       phiINT( 1 : M ) * Hank( 1 : M ) ) )
          ELSE                                  ! incoherent case
             P( ird, ir ) = CMPLX( SQRT( SUM( ABS(  phiINT( 1 : M ) * Hank( 1 : M ) ) ** 2 ) ), KIND=4 )
          ENDIF
       END DO Depths
    END DO Ranges

  END SUBROUTINE EvaluateAD

  !**********************************************************************!

  SUBROUTINE NewSegment( FileRoot, r, rProf, ifreq, NProf, kL, kR, phiL, phiR, M, rd, Nrd, sumk, sumkinv )

    ! Treats the crossing into a new range segment

    USE ReadModes

    INTEGER,            INTENT( IN    ) :: NProf, Nrd                      ! max # modes, # receiver depths
    REAL,               INTENT( IN    ) :: rd( Nrd )                       ! receiver depths
    REAL      (KIND=8), INTENT( IN    ) :: rProf( NProf + 1 )              ! profile ranges
    REAL      (KIND=8), INTENT( IN )    :: r( * )                          ! receiver ranges
    CHARACTER (LEN=80), INTENT( IN    ) :: FileRoot
    INTEGER,            INTENT( IN    ) :: ifreq                                 ! frequency index
    INTEGER,            INTENT( INOUT ) :: M                                     ! # modes
    COMPLEX,            INTENT( INOUT ) :: phiL( MaxM, Nrd ), phiR( MaxM, Nrd ), kL( M ), kR( MaxM )  ! modes at left and right of segment
    COMPLEX   (KIND=8), INTENT( INOUT ) :: sumk( M ), sumkinv( M )               ! wavenumbers and their inverse, integrated over range
    INTEGER                             :: MR                                    ! # mode in new, right segment
    COMPLEX                             :: kMid( M )

    ! Do phase integral up to the new range
    IF ( ir > 1 ) THEN
       rLeft = MAX( r( ir - 1 ), 1000.0 * rProf( iProf ) )
    ELSE
       rLeft = 1000.0 * rProf( iProf )
    ENDIF

    rMid             = 0.5 * ( 1000.0 * rProf( iProf + 1 ) + rLeft )
    wMid             = ( rMid / 1000.0 - rProf( iProf ) ) / ( rProf( iProf + 1 ) - rProf( iProf ) )
    kMid(    1 : M ) = CMPLX( kL(      1 : M ) + wMid  * ( kR( 1 : M ) - kL( 1 : M ) ) )
    sumk(    1 : M ) =        sumk(    1 : M ) + kMid( 1 : M ) * ( 1000.0 * rProf( iProf + 1 ) - rLeft )
    sumkinv( 1 : M ) =        sumkinv( 1 : M ) + ( 1000.0 * rProf( iProf + 1 ) - rLeft ) / kMid( 1 : M )

    ! Copy right modes to left
    kL(   1 : M )          = kR(   1 : M )
    phiL( 1 : M, 1 : Nrd ) = phiR( 1 : M, 1 : Nrd ) 

    ! Read in the new right mode set
    iProf = iProf + 1

    IF ( iProf + 1 <= NProf ) THEN   ! is there a new mode set to read?
       CALL GetModes( FileRoot, iProf + 1, ifreq, MaxM, rd, Nrd, 'N', kR, phiR, MR, freqVec, Nfreq, Title )
       M = MIN( M, MR )
       WRITE( *, * ) 'New profile read', r( ir ), iProf + 1, rProf( iProf + 1 ), ' #modes=', M
    ENDIF

  END SUBROUTINE NewSegment
END MODULE EvaluateADMod
