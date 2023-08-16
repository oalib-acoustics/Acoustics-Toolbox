MODULE ReflectMod

  USE bellhopMod
  IMPLICIT NONE
CONTAINS

  SUBROUTINE Reflect2D( is, HS, BotTop, nBdry3d, z_xx, z_xy, z_yy, kappa_xx, kappa_xy, kappa_yy, RefC, Npts, tradial )

    !USE norms
    USE RefCoef
    USE sspMod
    USE Cone   ! if using analytic formulas for the cone (conical seamount)

    INTEGER,              INTENT( IN ) :: Npts
    REAL     (KIND=8),    INTENT( IN ) :: tradial( 2 )
    REAL     (KIND=8),    INTENT( IN ) :: nBdry3d( 3 )                    ! Normal to the boundary
    CHARACTER (LEN=3),    INTENT( IN ) :: BotTop                          ! Flag indicating bottom or top reflection
    TYPE( HSInfo ),       INTENT( IN ) :: HS                              ! half-space properties
    TYPE(ReflectionCoef), INTENT( IN ) :: RefC( NPts )                    ! reflection coefficient
    INTEGER,              INTENT( INOUT ) :: is
    INTEGER           :: is1
    REAL     (KIND=8) :: c, cimag, gradc( 2 ), crr, crz, czz, rho         ! derivatives of sound speed
    REAL     (KIND=8) :: RM, RN, Tg, Th, rayt( 2 ), rayn( 2 ), rayt_tilde( 2 ), rayn_tilde( 2 ), cnjump, csjump  ! for curvature change
    REAL     (KIND=8) :: ck, co, si, cco, ssi, pdelta, rddelta, sddelta   ! for beam shift
    COMPLEX  (KIND=8) :: gamma1, gamma2, gamma1Sq, gamma2Sq, GK, Refl, ch, a, b, d, sb, delta, ddelta
    REAL     (KIND=8) :: kappa                                            ! Boundary curvature
    REAL     (KIND=8) :: tBdry( 2 ), nBdry( 2 )                           ! tangent, normal to boundary
    TYPE(ReflectionCoef) :: RInt
    REAL     (KIND=8) :: z_xx, z_xy, z_yy, &
         t_rot( 2 ), n_rot( 2 ), RotMat( 2, 2 ), kappaMat( 2, 2 ), DMat( 2, 2 )   ! for cone reflection
    REAL     (KIND=8) :: kappa_xx, kappa_xy, kappa_yy

    is  = is + 1
    is1 = is + 1

    ! following should perhaps be pre-calculated in Bdry3DMod
    nBdry( 1 ) = DOT_PRODUCT( nBdry3d( 1 : 2 ), tradial )
    nBdry( 2 ) = nBdry3d( 3 )

    !CALL ConeFormulas( z_xx, z_xy, z_yy, nBdry, xs_3D, tradial, ray2D( is )%x, BotTop )   ! special case of a conical seamount

    !!!! use kappa_xx or z_xx?
    ! kappa = ( z_xx * tradial( 1 ) ** 2 + 2 * z_xy * tradial( 1 ) * tradial( 2 ) + z_yy * tradial( 2 ) ** 2 )
    kappa = ( kappa_xx * tradial( 1 ) ** 2 + 2 * kappa_xy * tradial( 1 ) * tradial( 2 ) + kappa_yy * tradial( 2 ) ** 2 )

    IF ( BotTop == 'TOP' ) kappa = -kappa

    nBdry = nBdry / NORM2( nBdry )

    Th = DOT_PRODUCT( ray2D( is )%t, nBdry )  ! component of ray tangent normal to boundary

    tBdry = ray2D( is )%t - Th * nBdry        ! component of ray tangent along the boundary
    tBdry = tBdry / NORM2( tBdry )
    ! could also calculate tbdry as +/- of [ nbdry( 2), -nbdry( 1 ) ], but need sign

    Tg = DOT_PRODUCT( ray2D( is )%t, tBdry )  ! component of ray tangent along boundary

    ray2D( is1 )%NumTopBnc = ray2D( is )%NumTopBnc
    ray2D( is1 )%NumBotBnc = ray2D( is )%NumBotBnc
    ray2D( is1 )%x         = ray2D( is )%x
    ray2D( is1 )%t         = ray2D( is )%t - 2.0 * Th * nBdry  ! changing the ray direction

    ! Calculate the change in curvature
    ! Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).

    CALL EvaluateSSP2D( ray2D( is1 )%x, c, cimag, gradc, crr, crz, czz, rho, xs_3D, tradial, freq )

    ! incident and reflected (tilde) unit ray tangent and normal
    rayt       = c * ray2D( is  )%t                       ! unit tangent to ray
    rayt_tilde = c * ray2D( is1 )%t                       ! unit tangent to ray
    rayn       = +[ -rayt(       2 ), rayt(       1 ) ]   ! unit normal  to ray
    rayn_tilde = -[ -rayt_tilde( 2 ), rayt_tilde( 1 ) ]   ! unit normal  to ray

    RN = 2 * kappa / c ** 2 / Th    ! boundary curvature correction

    ! get the jumps (this could be simplified, e.g. jump in rayt is roughly 2 * Th * nbdry
    cnjump = -DOT_PRODUCT( gradc, rayn_tilde - rayn  )
    csjump = -DOT_PRODUCT( gradc, rayt_tilde - rayt )

    IF ( BotTop == 'TOP' ) THEN
       cnjump = -cnjump    ! flip sign for top reflection
       RN     = -RN
    END IF

    RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
    RN = RN + RM * ( 2 * cnjump - RM * csjump ) / c ** 2

    SELECT CASE ( Beam%Type( 2 : 2 ) )
    CASE ( 'D' )
       RN = 2.0 * RN
    CASE ( 'Z' )
       RN = 0.0
    END SELECT

    ray2D( is1 )%c   = c
    ray2D( is1 )%tau = ray2D( is )%tau
    ray2D( is1 )%p   = ray2D( is )%p + ray2D( is )%q * RN
    ray2D( is1 )%q   = ray2D( is )%q

    ! amplitude and phase change

    SELECT CASE ( HS%BC )
    CASE ( 'R' )                 ! rigid
       ray2D( is1 )%Amp   = ray2D( is )%Amp
       ray2D( is1 )%Phase = ray2D( is )%Phase
    CASE ( 'V' )                 ! vacuum
       ray2D( is1 )%Amp   = ray2D( is )%Amp
       ray2D( is1 )%Phase = ray2D( is )%Phase + pi
    CASE ( 'F' )                 ! file
       RInt%theta = RadDeg * ABS( ATAN2( Th, Tg ) )           ! angle of incidence (relative to normal to bathymetry)
       IF ( RInt%theta > 90 ) RInt%theta = 180. - RInt%theta  ! reflection coefficient is symmetric about 90 degrees
       CALL InterpolateReflectionCoefficient( RInt, RefC, Npts, PRTFile )
       ray2D( is1 )%Amp   = ray2D( is )%Amp   * RInt%R
       ray2D( is1 )%Phase = ray2D( is )%Phase + RInt%phi
    CASE ( 'A', 'G' )            ! half-space
       GK       = omega * Tg     ! wavenumber in direction parallel to bathymetry
       gamma1Sq = ( omega / c     ) ** 2 - GK ** 2 - i * tiny( omega )   ! tiny prevents g95 giving -zero, and wrong branch cut
       gamma2Sq = ( omega / HS%cP ) ** 2 - GK ** 2 - i * tiny( omega )
       gamma1   = SQRT( -gamma1Sq )
       gamma2   = SQRT( -gamma2Sq )

       Refl = ( HS%rho * gamma1 - rho * gamma2 ) / ( HS%rho * gamma1 + rho * gamma2 )
       ! write( *, * ) abs( Refl ), c, HS%cp, rho, HS%rho
       IF ( ABS( Refl ) < 1.0E-5 ) THEN   ! kill a ray that has lost its energy in reflection
          ray2D( is1 )%Amp   = 0.0
          ray2D( is1 )%Phase = ray2D( is )%Phase
       ELSE
          ray2D( is1 )%Amp   = ABS( Refl ) * ray2D(  is )%Amp
          ray2D( is1 )%Phase = ray2D( is )%Phase + ATAN2( AIMAG( Refl ), REAL( Refl ) )

          ! compute beam-displacement Tindle, Eq. (14)
          ! needs a correction to beam-width as well ...
          !  IF ( REAL( gamma2Sq ) < 0.0 ) THEN
          !     rhoW   = 1.0   ! density of water
          !     rhoWSq  = rhoW  * rhoW
          !     rhoHSSq = rhoHS * rhoHS
          !     DELTA = 2 * GK * rhoW * rhoHS * ( gamma1Sq - gamma2Sq ) /
          ! &( gamma1 * i * gamma2 *
          ! &( -rhoWSq * gamma2Sq + rhoHSSq * gamma1Sq ) )
          !     RV( is + 1 ) = RV( is + 1 ) + DELTA
          !  END IF

          if ( Beam%Type( 4 : 4 ) == 'S' ) then   ! beam displacement & width change (Seongil's version)

             ch = ray2D( is )%c / conjg( HS%cP )
             co = ray2D( is )%t( 1 ) * ray2D( is )%c
             si = ray2D( is )%t( 2 ) * ray2D( is )%c
             ck = omega / ray2D( is )%c

             a   = 2 * HS%rho * ( 1 - ch * ch )
             b   = co * co - ch * ch
             d   = HS%rho * HS%rho * si * si + b
             sb  = sqrt( b )
             cco = co * co
             ssi = si * si

             IF ( si /= 0.0 ) THEN
                delta = a * co / si / ( ck * sb * d )   ! Do we need an abs() on this???
             ELSE
                delta = 0.0
             END IF

             pdelta  = real( delta ) / ( ray2D( is )%c / co)
             ddelta  = -a / ( ck*sb*d ) - a*cco / ssi / (ck*sb*d) + a*cco / (ck*b*sb*d) &
                  -a*co / si / (ck*sb*d*d) * (2* HS%rho * HS%rho *si*co-2*co*si)
             rddelta = -real( ddelta )
             sddelta = rddelta / abs( rddelta )        

             print *, 'beam displacing ...'
             ray2D( is1 )%x( 1 ) = ray2D( is1 )%x( 1 ) + real( delta )   ! displacement
             ray2D( is1 )%tau    = ray2D( is1 )%tau + pdelta             ! phase change
             ray2D( is1 )%q      = ray2D( is1 )%q + sddelta * rddelta * si * c * ray2D( is )%p   ! beam-width change
          endif

       ENDIF
    END SELECT
  END SUBROUTINE Reflect2D

END MODULE ReflectMod
