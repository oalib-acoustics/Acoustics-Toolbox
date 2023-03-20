MODULE BCImpedancecMod

  USE KrakencMod
  USE sspMod
  USE PekRoot

  IMPLICIT NONE

  INTEGER,          PARAMETER, PRIVATE :: iPowerR = 50, iPowerF = -50
  REAL    (KIND=8), PARAMETER, PRIVATE :: Roof = 1.0D+50, Floor = 1.0D-50
  INTEGER, PRIVATE          :: ii, j
  REAL    (KIND=8) :: two_h
  COMPLEX (KIND=8) :: two_x, xB3, four_h_x

  CONTAINS
  SUBROUTINE BCImpedance( x, BotTop, HS, f, g, iPower )

    ! Compute Boundary Condition Impedance

    USE RefCoef
    COMPLEX  (KIND=8), PARAMETER     :: zero = ( 0.0D0, 0.0D0 )
    COMPLEX  (KIND=8), INTENT( IN  ) :: x           ! wavenumber squared
    CHARACTER (LEN=3), INTENT( IN  ) :: BotTop      ! flag indicating bottom or top boundary
    TYPE( HSInfo ),    INTENT( IN  ) :: HS          ! structure containing halfspace parameters
    INTEGER,           INTENT( OUT ) :: iPower      ! powers of 10 for large dynamic range of the impedance
    COMPLEX  (KIND=8), INTENT( OUT ) :: f, g        ! functions defining the impedance for that value of x
    INTEGER           :: iTop = 0, iBot = 0, Medium
    REAL     (KIND=8) :: rhoInside = 1.
    COMPLEX  (KIND=8) :: gammaS, gammaP, gammaS2, gammaP2, mu, yV( 5 )
    COMPLEX  (KIND=8) :: kx, kz, RCmplx, cInside = 1500.
    TYPE(ReflectionCoef) :: RInt

    iPower = 0

    ! Get rho, c just inside the boundary
    ! There is at least one acoustic layer in the problem, except
    ! in the case where BOUNCE is used to get the refl. coef. for a purely elastic stack of layers.
    ! These are initialized above just to avoid a compiler warning

    SELECT CASE ( BotTop )
    CASE ( 'TOP' )
       IF ( FirstAcoustic > 0 ) THEN
          iTop      = Loc( FirstAcoustic ) + N( FirstAcoustic ) + 1
          rhoInside = rho( iTop )
          cInside   = PekerisRoot( omega2 * h( FirstAcoustic ) **2 / ( 2.0D0 + B1( iTop ) ) )
       END IF
    CASE ( 'BOT' )
       IF ( LastAcoustic > 0 ) THEN
          iBot      = Loc( LastAcoustic ) + N( LastAcoustic ) + 1
          rhoInside = rho( iBot )
          cInside   = PekerisRoot( omega2 * h( LastAcoustic ) **2 / ( 2.0D0 + B1( iBot ) ) )
       END IF
    END SELECT

    ! Return the impedance, depending on the type of boundary condition

    SELECT CASE ( HS%BC )
    CASE ( 'V' )                      ! Vacuum
       f  = 1.0D0
       !g = -i * PekerisRoot( omega2 / cInside ** 2 - x ) * SIGMA( 1 ) ** 2
       g  = 0.0D0
       yV = CMPLX( [ f, g, zero, zero, zero ] )
    CASE ( 'R' )                      ! Rigid
       f  = 0.0D0
       g  = 1.0D0
       yV = CMPLX( [ f, g, zero, zero, zero ] )
    CASE ( 'A' )                      ! Acousto-elastic half-space
       IF ( REAL( HS%cS ) > 0.0 ) THEN
          gammaS2 = x - omega2 / HS%cS ** 2
          gammaP2 = x - omega2 / HS%cP ** 2
          gammaS  = PekerisRoot( gammaS2 )
          gammaP  = PekerisRoot( gammaP2 )
          mu      = HS%rho * HS%cS ** 2

          yV( 1 ) = ( gammaS * gammaP - x ) / mu
          yV( 2 ) = ( ( gammaS2 + x ) ** 2 - 4.0D0 * gammaS * gammaP * x ) * mu
          yV( 3 ) = 2.0D0 * gammaS * gammaP - gammaS2 - x
          yV( 4 ) = gammaP * ( x - gammaS2 )
          yV( 5 ) = gammaS * ( gammaS2 - x )

          f = omega2 * yV( 4 )
          g = yV( 2 )
       ELSE
          gammaP = PekerisRoot( x - omega2 / HS%cP ** 2 )
          f      = gammaP
          g      = HS%rho
       END IF
    CASE ( 'F' )                    ! Tabulated reflection coefficient
       ! Compute the grazing angle, theta
       kx         = SQRT( x )
       kz         = SQRT( omega2 / cInside ** 2 - x )
       RInt%theta = RadDeg * DATAN2( DBLE( kz ), DBLE( kx ) )

       ! Evaluate R( ThetaInt )
       IF ( BotTop == 'TOP' ) THEN
          CALL InterpolateReflectionCoefficient( RInt, RTop, NTopPts, PRTFile )
       ELSE
          CALL InterpolateReflectionCoefficient( RInt, RBot, NBotPts, PRTFile )
       END IF

       ! Convert R(theta) to (f,g) in Robin BC
       RCmplx = RInt%R * EXP( i * RInt%phi )
       f      = 1.0D0
       g      = ( 1.0D0 + RCmplx ) / ( i * kz * ( 1.0D0 - RCmplx ) )     
    CASE ( 'P' )                    ! Precalculated reflection coef
       CALL InterpolateIRC( x, f, g, iPower, xTab, fTab, gTab, iTab, NkTab )
    END SELECT

    IF ( BotTop == 'TOP' ) g = -g    ! A top BC has the sign flipped relative to a bottom BC

    ! Shoot through elastic layers
    SELECT CASE ( BotTop )
    CASE ( 'TOP' )
       IF ( FirstAcoustic > 1 ) THEN   ! Shoot down from top

          DO Medium = 1, FirstAcoustic - 1
             CALL ElasticDN( x, yV, iPower, Medium )
          END DO

          f = omega2 * yV( 4 )
          g = yV( 2 )
       END IF
    CASE ( 'BOT' )
       IF ( LastAcoustic < SSP%NMedia ) THEN   ! Shoot up from bottom

          DO Medium = SSP%NMedia, LastAcoustic + 1, -1
             CALL ElasticUP( x, yV, iPower, Medium )
          END DO

          f = omega2 * yV( 4 )
          g = yV( 2 )
       END IF
    END SELECT

  END SUBROUTINE BCImpedance
  !**********************************************************************!
  SUBROUTINE ElasticUP( x, yV, iPower, Medium )

    ! Propagates through an elastic layer using compound matrix formulation

    INTEGER,          INTENT( IN    ) :: Medium
    INTEGER,          INTENT( INOUT ) :: iPower
    COMPLEX (KIND=8), INTENT( IN    ) :: x                ! trial eigenvalue, k2
    COMPLEX (KIND=8) :: xV( 5 ), yV( 5 ), zV( 5 )   ! solution of differential equation at 3 successive steps

    ! Euler's method for first step

    two_x    = 2.0D0 * x
    two_h    = 2.0D0 * h( Medium )
    four_h_x = 4.0D0 * h( Medium ) * x
    j        = Loc( Medium ) + N( Medium ) + 1
    xB3      = x * B3( j ) - rho( j )

    zV( 1 ) = yV( 1 ) - 0.5D0 * (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
    zV( 2 ) = yV( 2 ) - 0.5D0 * ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
    zV( 3 ) = yV( 3 ) - 0.5D0 * (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
    zV( 4 ) = yV( 4 ) - 0.5D0 * (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
    zV( 5 ) = yV( 5 ) - 0.5D0 * (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

    ! Modified midpoint method

    DO ii = N( Medium ), 1, -1
       j = j - 1

       xV = yV
       yV = zV

       xB3 = x * B3( j ) - rho( j )

       zV( 1 ) = xV( 1 ) - (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
       zV( 2 ) = xV( 2 ) - ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
       zV( 3 ) = xV( 3 ) - (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
       zV( 4 ) = xV( 4 ) - (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
       zV( 5 ) = xV( 5 ) - (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

       ! Scale if necessary
       IF ( ii /= 1 ) THEN
          IF ( ABS( DBLE( zV( 2 ) ) ) < Floor ) THEN
             zV     = Roof * zV
             yV     = Roof * yV
             iPower = iPower - iPowerR
          END IF

          IF ( ABS( DBLE( zV( 2 ) ) ) > Roof  ) THEN
             zV     = Floor * zV
             yV     = Floor * yV
             iPower = iPower - iPowerF
          END IF
       END IF
    END DO

    yV = ( xV + 2.0D0 * yV + zV ) / 4.0D0   ! Apply the standard filter at the terminal point

  END SUBROUTINE ElasticUP
  !**********************************************************************!
  SUBROUTINE ElasticDN( x, yV, iPower, Medium )

    ! Propagates through an elastic layer using compound matrix formulation

    INTEGER,          INTENT( IN    ) :: Medium
    INTEGER,          INTENT( INOUT ) :: iPower
    COMPLEX (KIND=8), INTENT( IN    ) :: x                ! trial eigenvalue, k2
    COMPLEX (KIND=8) :: xV( 5 ), yV( 5 ), zV( 5 )   ! solution of differential equation at 3 successive steps

    ! Euler's method for first step

    two_x    = 2.0D0 * x
    two_h    = 2.0D0 * h( Medium )
    four_h_x = 4.0D0 * h( Medium ) * x
    j        = Loc( Medium ) + 1
    xB3      = x * B3( j ) - rho( j )

    zV( 1 ) = yV( 1 ) + 0.5D0 * (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
    zV( 2 ) = yV( 2 ) + 0.5D0 * ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
    zV( 3 ) = yV( 3 ) + 0.5D0 * (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
    zV( 4 ) = yV( 4 ) + 0.5D0 * (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
    zV( 5 ) = yV( 5 ) + 0.5D0 * (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

    ! Modified midpoint method

    DO ii = 1, N( Medium )
       j = j + 1

       xV = yV
       yV = zV

       xB3 = x * B3( j ) - rho( j )

       zV( 1 ) = xV( 1 ) + (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
       zV( 2 ) = xV( 2 ) + ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
       zV( 3 ) = xV( 3 ) + (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
       zV( 4 ) = xV( 4 ) + (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
       zV( 5 ) = xV( 5 ) + (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

       ! Scale if necessary
       IF ( ii /= N( Medium ) ) THEN
          IF ( ABS( DBLE( zV( 2 ) ) ) < Floor ) THEN
             zV     = Roof * zV
             yV     = Roof * yV
             iPower = iPower - iPowerR
          END IF

          IF ( ABS( DBLE( zV( 2 ) ) ) > Roof  ) THEN
             zV     = Floor * zV
             yV     = Floor * yV
             iPower = iPower - iPowerF
          END IF
       END IF
    END DO

    yV = ( xV + 2.0D0 * yV + zV ) / 4.0D0   ! Apply the standard filter at the terminal point

  END SUBROUTINE ElasticDN

END MODULE BCImpedancecMod
