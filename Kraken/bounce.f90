PROGRAM BOUNCE

  ! Program for computing the reflection coefficient vs angle, R( theta )
  ! Michael B. Porter

  USE KrakencMod
  USE ReadEnvironmentMod
  USE RefCoef
  USE sspMod
  USE AttenMod
  USE PekRoot
  USE FatalError

  IMPLICIT NONE
  INTEGER            :: IAllocStat, ii
  REAL      (KIND=8) :: kMin, kMax, freq0
  CHARACTER (LEN=80) :: FileRoot

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  Title = 'BOUNCE-'
  iProf = 1
  CALL ReadEnvironment( FileRoot, Title, freq, MaxMedium, TopOpt, NG, BotOpt, cLow, cHigh, RMax, ENVFile, PRTFile )

  IF ( cLow <= 0.0 ) CALL ERROUT( 'BOUNCE', 'cLow must be positive' )

  CLOSE( ENVFile )

  omega  = 2.0D0 * pi * freq
  omega2 = omega ** 2

  CALL ReadReflectionCoefficient( FileRoot, BotOpt( 1 : 1 ), TopOpt( 1 : 1 ), PRTFile )

  freq0 = freq   ! save the reference frequency for scaling the grid
  CALL UpdateSSPLoss( freq, freq0 )
  CALL UpdateHSLoss(  freq, freq0 )

  WRITE( PRTFile, * ) 'Writing internal reflection coeffecient table'

  ! Compute NkTab, the number of wavenumbers for which we will tabulate the reflection loss

  kMin = omega / cHigh
  kMax = omega / cLow
  IF ( cHigh > 1.0E6 ) kMin = 0.0

  NkTab = INT( 1000.0 * RMax * ( kMax - kMin ) / ( 2.0 * pi ) )
  WRITE( PRTFile, * ) 'NkTab = ', NkTab

  ALLOCATE( xTab( NkTab ), fTab( NkTab ), gTab( NkTab), ITab( NkTab ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( 'BOUNCE', 'Too many points in reflection coefficient' )

  CALL ComputeReflectionCoefficient
  CLOSE( PRTFile )

CONTAINS

  !**********************************************************************

  SUBROUTINE Initialize

    ! Initializes arrays defining difference equations

    LOGICAL           :: ElasticFlag = .FALSE.
    INTEGER           :: iAllocStat, j, Medium, N1, NPts
    REAL    (KIND=8)  :: Twoh
    COMPLEX (KIND=8)  :: cP2, cS2
    COMPLEX (KIND=8), ALLOCATABLE :: cP( : ), cS( : )
    CHARACTER (LEN=8) :: Task

    cmin          = 1.0E6
    FirstAcoustic = 0
    Loc( 1 )      = 0
    NPTS          = SUM( N( 1 : SSP%NMedia ) ) + SSP%NMedia

    ALLOCATE ( B1( NPTS ), B2( NPTS ), B3( NPTS ), B4( NPTS ), rho( NPTS ), CP( NPTS ), cS( NPTS ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( 'BOUNCE - Initialize', 'Insufficient memory: Reduce mesh.' )

    MediumLoop: DO Medium = 1, SSP%NMedia
       IF ( Medium /= 1 ) Loc( Medium ) = Loc( Medium - 1 ) + N( Medium - 1 ) + 1
       N1 = N(   Medium ) + 1
       ii = Loc( Medium ) + 1

       ! EvaluateSSP reads in the data for a medium

       Task = 'TAB'
       CALL EvaluateSSP( cP( ii ), cS( ii ), rho( ii ), Medium, N1, freq, Task )

       ! Load diagonals of the finite-difference equations

       IF ( cS( ii ) == ( 0.0, 0.0 ) ) THEN ! Case of an acoustic medium

          SSP%Material( Medium ) = 'ACOUSTIC'
          IF ( FirstAcoustic == 0 ) FirstAcoustic = Medium
          LastAcoustic = Medium

          cmin = MIN( MINVAL( DBLE( cP( ii : ii + N( Medium ) ) ) ), cmin )
          B1( ii : ii + N( Medium ) ) = -2.0 + h( Medium ) ** 2 * omega2 / cP( ii : ii + N( Medium ) ) ** 2

       ELSE                                 ! Case of an elastic medium

          IF ( SSP%sigma( Medium ) /= 0.0 ) &
             CALL ERROUT( 'BOUNCE - Initialize', 'Rough elastic interfaces are not allowed' )

          SSP%Material( Medium ) = 'ELASTIC'
          ElasticFlag = .TRUE.
          Twoh   = 2.0 * h( Medium )

          DO j = ii, ii + N( Medium )
             cmin = MIN( DBLE( cS( j ) ), cmin )

             cP2 = cP( j ) ** 2
             cS2 = cS( j ) ** 2

             B1(  j ) = Twoh / ( rho( j ) * cS2 )
             B2(  j ) = Twoh / ( rho( j ) * cP2 )
             B3(  j ) = 4.0 * Twoh * rho( j ) * cS2 * ( cP2 - cS2 ) / cP2
             B4(  j ) = Twoh * ( cP2 - 2.0 * cS2 ) / cP2
             rho( j ) = Twoh * omega2 * rho( j )
          END DO
       END IF
    END DO MediumLoop

    ! Bottom properties

    IF ( HSBot%BC( 1 : 1 ) == 'A' ) THEN
       IF ( HSBot%cS /= ( 0.0, 0.0 ) ) THEN   ! Elastic bottom:
          ElasticFlag = .TRUE.
          cmin = MIN( cmin, DBLE(  HSBot%cS ) )
       ELSE                                   ! Acoustic bottom:
          cmin = MIN( cmin, DBLE(  HSBot%cP ) )
       END IF
    ENDIF

    ! Top properties

    IF ( HSTop%BC( 1 : 1 ) == 'A' ) THEN
       IF (  HSTop%cS /= ( 0.0, 0.0 ) ) THEN   ! Elastic top:
          ElasticFlag = .TRUE.
          cmin = MIN( cmin, DBLE(  HSTop%cS ) )
       ELSE                                    ! Acoustic top:
          cmin = MIN( cmin, DBLE(  HSTop%cP ) )
       END IF
    END IF

    IF ( ElasticFlag ) cmin = 0.9 * cmin
    cLow = MAX( cLow, 0.99 * cmin )

  END SUBROUTINE Initialize

  !**********************************************************************

  SUBROUTINE ComputeReflectionCoefficient

    ! Computes the reflection coefficient for k in [ kmin, kmax ]

    USE BCImpedancecMod
    INTEGER                  :: ik, IPow, itheta, Loops
    REAL      (KIND=8)       :: k0, c0, Deltak, k1
    COMPLEX   (KIND=8)       :: x, f, g
    REAL      (KIND=8)       :: kx( NkTab ), kz( NkTab ), theta( NkTab ), R( NkTab ), phase( NkTab )
    COMPLEX   (KIND=8)       :: RCmplx( NkTab ), R1, R2

    N( 1 : SSP%NMedia ) = NG( 1 : SSP%NMedia )
    h( 1 : SSP%NMedia ) = ( SSP%Depth( 2 : SSP%NMedia + 1 ) - SSP%Depth( 1 : SSP%NMedia ) ) / N( 1 : SSP%NMedia )

    HV( 1 ) = h( 1 )
    CALL Initialize

    Deltak = ( kMax - kMin ) / ( NkTab - 1 )

    Wavenumber: DO ik = 1, NkTab
       k1 = kMin + ( ik - 1 ) * Deltak
       x  = k1 ** 2

       CALL BCImpedance( x, 'BOT', HSBot, f, g, IPow )  ! Bottom impedance
       CALL AcousticLayers( x, f, g, IPow  )            ! Shoot through acoustic layers
       xTab( ik ) = DBLE( x )
       fTab( ik ) = f
       gTab( ik ) = g
       ITab( ik ) = IPow
    END DO Wavenumber

    IF ( HSTop%BC( 1 : 1 ) == 'A' ) THEN
       c0 = DBLE( HSTop%cP )                    ! use upper halfspace speed for reference if a halfspace was specified
    ELSE
       c0 = 1500
    END IF

    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Reference sound speed = ', c0

    k0 = omega / c0   ! free-space wavenumber
    kx = SQRT( xTab ) ! horizontal wavenumber

    WHERE( k0 > kx )  ! propagating waves
       kz     = SQRT( k0 ** 2 - kx ** 2 )   ! vertical wavenumber
       theta  = RadDeg * ATAN2( kz, kx )    ! angle of incidence
       RCmplx =  - ( fTab - i * kz * gTab ) / ( fTab + i * kz * gTab )   ! complex reflection coef.
       R      = ABS( RCmplx )
       phase  = RadDeg * ATAN2( AIMAG( RCmplx ), REAL( RCmplx ) )
    ELSEWHERE         ! evanescent waves
       kz     = 0
       theta  = 0
       RCmplx = 1.0   ! complex reflection coef.
       R      = 1.0
       phase  = 180.0
    END WHERE

    ! unwrap the phase by counting loops in the complex plane
    Loops = 0

    Angle: DO itheta = NkTab - 1, 1, -1
       R1 = RCmplx( itheta )
       R2 = RCmplx( itheta + 1 )

       IF ( AIMAG( R1 ) > 0 .AND. AIMAG( R1 ) < 0 .AND. REAL( R2 ) < 0 ) Loops = Loops + 1
       IF ( AIMAG( R1 ) < 0 .AND. AIMAG( R2 ) > 0 .AND. REAL( R2 ) < 0 ) Loops = Loops - 1

       phase( itheta ) = phase( itheta ) - Loops * 360
    END DO Angle

    OPEN( FILE = TRIM( FileRoot ) // '.irc', UNIT = IRCFile, STATUS = 'UNKNOWN' )   ! Internal Reflection Coef. format
    WRITE( IRCFile, * ) '''', Title, '''', freq
    WRITE( IRCFile, * ) NkTab
    WRITE( IRCFile, FMT = "( 5G15.7, I5 )" ) ( xTab( ik ), fTab( ik ), gTab( ik ), ITab( ik ), ik = 1, NkTab )

    OPEN( FILE = TRIM( FileRoot ) // '.brc', UNIT = BRCFile, STATUS = 'UNKNOWN' )   ! Bottom   Reflection Coef. format
    !WRITE( BRCFile, * ) Title
    !WRITE( BRCFile, * ) freq
    WRITE( BRCFile, * ) NkTab
    DO ik = NkTab, 1, -1
       WRITE( BRCFile, * ) theta( ik ), R( ik ), phase( ik )
    END DO

    CLOSE( IRCFile )
    CLOSE( BRCFile )

  END SUBROUTINE ComputeReflectionCoefficient

  !**********************************************************************

  SUBROUTINE AcousticLayers( x, F, G, IPow )

    ! Shoot through acoustic layers

    INTEGER,       PARAMETER :: IPowF = -5
    REAL (KIND=8), PARAMETER :: Roof = 1.0E5, Floor = 1.0E-5
    INTEGER,          INTENT( INOUT ) :: IPow
    COMPLEX (KIND=8), INTENT( IN    ) :: x
    COMPLEX (KIND=8), INTENT( INOUT ) :: f, g
    INTEGER                  :: Medium
    REAL (KIND=8)            :: rhoM = 1.0
    COMPLEX (KIND=8)         :: P0, P1, P2, h2k2

    IF ( FirstAcoustic == 0 ) RETURN

    MediumLoop: DO Medium = LastAcoustic, FirstAcoustic, -1
       h2k2 = h( Medium ) ** 2 * x
       ii   = Loc( Medium ) + N( Medium ) + 1
       rhoM = rho( ii )
       P1   = -2.0 * g
       P2   = ( B1( ii ) - h2k2 ) * g - 2.0 * h( Medium ) * f * rhoM

       ! Shoot through a single medium

       DO ii = Loc( Medium ) + N( Medium ), Loc( Medium ) + 1, -1

          P0 = P1
          P1 = P2
          P2 = ( h2k2 - B1( ii ) ) * P1 - P0

          DO WHILE ( ABS( DBLE( P2 ) ) > Roof )
             P0   = Floor * P0
             P1   = Floor * P1
             P2   = Floor * P2
             IPow = IPow - IPowF
          END DO
       END DO

       ! g = P'/rho and g = -P since fP + gP'/rho = 0
       f = -( P2 - P0 ) / ( 2.0 * h( Medium ) ) / rhoM
       g = -P1
    END DO MediumLoop

  END SUBROUTINE AcousticLayers

END PROGRAM BOUNCE
