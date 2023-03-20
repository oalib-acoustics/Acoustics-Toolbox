MODULE Influence3D

  ! Compute the beam influence, i.e. the contribution of a single beam to the complex pressure
  ! mbp 12/2018, based on much older subroutines

  USE bellhopMod
  USE WriteRay
  USE RayNormals
  USE SourceReceiverPositions
  USE ArrMod
  USE sspMod   ! used to construct image beams in the Cerveny style beam routines
  USE cross_products

  IMPLICIT NONE
  INTEGER,       PRIVATE :: itheta, iz, ir, is
  REAL (KIND=8), PRIVATE :: W, s, m, n, Amp, phase, const, phaseInt, Ratio1, &
       L1, L2, rayt( 3 ), &
       SrcDeclAngle, RcvrDeclAngle, SrcAzimAngle, RcvrAzimAngle, &
       rA, rB, lambda
  COMPLEX ( KIND=8 ) :: delay

CONTAINS
  SUBROUTINE Influence3DGeoHatRayCen( xs, alpha, beta, Dalpha, Dbeta, P )

  ! Geometrically-spreading beams with a hat-shaped beam in ray-centered coordinates

  REAL ( KIND=8 ), INTENT( IN  ) :: alpha, beta, Dalpha, Dbeta         ! ray take-off angle
  REAL ( KIND=8 ), INTENT( IN  ) :: xs( 3 )                            ! source coordinate
  COMPLEX        , INTENT( OUT ) :: P( Pos%Ntheta, Pos%Nrz, Pos%NRr )  ! complex pressure
  INTEGER          :: irA, irB, II
  REAL    (KIND=8) :: nA, nB, mA, mB, DS, deltaA, deltaB, t_rcvr( 2, Pos%Ntheta ), &
       a, b, KMAHphase( Beam%Nsteps ), DetQint
  REAL    (KIND=8) :: e1G( 3, Beam%Nsteps ), e2G( 3, Beam%Nsteps ), e1xe2( 3, Beam%Nsteps ), &
       xt( 3, Beam%Nsteps ),  xtxe1( 3, Beam%Nsteps ),  xtxe2( 3, Beam%Nsteps )
  REAL    (KIND=8) :: q_tilde( 2 ), q_hat( 2 ), DetQ( Beam%Nsteps ), dq_tilde( 2 ), dq_hat( 2 )
  COMPLEX (KIND=8) :: dtau( Beam%Nsteps - 1 )

  SrcDeclAngle = RadDeg * alpha          ! take-off declination angle in degrees
  SrcAzimAngle = RadDeg * beta           ! take-off azimuthal   angle in degrees

  DS     = SQRT( 2.0 ) * SIN( omega * xs( 3 ) * ray3D( 1 )%t( 3 ) )   ! Lloyd mirror pattern (for semi-coherent source)
  Ratio1 = SQRT( ABS( COS( alpha ) ) ) * SQRT( Dalpha * Dbeta ) / ray3D( 1 )%c

  DetQ     = ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) * ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) - &
             ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) * ray3D( 1 : Beam%Nsteps )%q_hat(   1 )

  dtau     = ray3D(   2 : Beam%Nsteps )%tau - ray3D(   1 : Beam%Nsteps - 1 )%tau
  
  ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) =                       Dalpha * ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) / ray3D( 1 )%c
  ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) =                       Dalpha * ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) / ray3D( 1 )%c
  ray3D( 1 : Beam%Nsteps )%q_hat(   1 ) = ABS( COS( ALPHA ) ) * Dbeta  * ray3D( 1 : Beam%Nsteps )%q_hat(   1 ) / ray3D( 1 )%c
  ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) = ABS( COS( ALPHA ) ) * Dbeta  * ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) / ray3D( 1 )%c

  ! phase shifts at caustics (rotations of Det_Q)
  KMAHphase( 1 ) = 0
  DO is = 2, Beam%Nsteps
     KMAHphase( is ) = KMAHphase( is - 1 )
     IF ( DetQ( is ) <= 0.0d0 .AND. DetQ( is - 1 ) > 0.0d0 .OR. &
          DetQ( is ) >= 0.0d0 .AND. DetQ( is - 1 ) < 0.0d0 ) KMAHphase( is ) = KMAHphase( is - 1 ) + pi / 2.
  END DO

  ! pre-calculate tangents, normals (e1, e2), etc.
  DO is = 1, Beam%Nsteps
     xt( :, is ) = ray3D( is )%x - xs   ! vector from the origin of the receiver plane to this point on the ray
     CALL RayNormal( ray3D( is )%t, ray3D( is )%phi, ray3D( is )%c, e1G( :, is ), e2G( :, is ) )

     ! e1xe2( :, is ) = cross_product( e1G( :, is ), e2G( :, is ) )
     e1xe2( :, is ) = ray3D( is )%c * ray3D( is )%t   ! unit tangent to ray
  END DO

  ! tangent along receiver bearing line
  t_rcvr( 1, : ) = COS( DegRad * Pos%theta( 1 : Pos%Ntheta ) )
  t_rcvr( 2, : ) = SIN( DegRad * Pos%theta( 1 : Pos%Ntheta ) )

  ReceiverDepths: DO iz = 1, Pos%Nrz
     ! precalculate vector from receiver to each step of ray
     xt( 3, : ) = ray3D( 1 : Beam%Nsteps )%x( 3 ) - Pos%rz( iz )
     DO is = 1, Beam%Nsteps
        xtxe1( :, is ) = cross_product( xt( :, is ), e1G( :, is ) )
        xtxe2( :, is ) = cross_product( xt( :, is ), e2G( :, is ) )
     END DO

     Radials: DO itheta = 1, Pos%Ntheta
        ! step along the beam ...
        ! Most of the time the beam makes no contribution to a receiver
        ! Therefore we try to test that quickly and move on to the next receiver

        Stepping: DO is = 2, Beam%Nsteps

           ! *** Compute coordinates of intercept: nB, mB, rB ***
           deltaA = -DOT_PRODUCT( t_rcvr( :, itheta ), e1xe2( 1 : 2, is - 1 ) )
           deltaB = -DOT_PRODUCT( t_rcvr( :, itheta ), e1xe2( 1 : 2, is     ) )

           ! Check for ray normal || radial of rcvr line
           IF ( ABS( deltaA ) < 1e-3 .OR. ABS( deltaB ) < 1e-3 )  THEN
              CYCLE Stepping
           END IF

           mA =  DOT_PRODUCT( t_rcvr( :, itheta ), xtxe1( 1 : 2, is - 1 ) ) / deltaA
           mB =  DOT_PRODUCT( t_rcvr( :, itheta ), xtxe1( 1 : 2, is     ) ) / deltaB

           nA = -DOT_PRODUCT( t_rcvr( :, itheta ), xtxe2( 1 : 2, is - 1 ) ) / deltaA
           nB = -DOT_PRODUCT( t_rcvr( :, itheta ), xtxe2( 1 : 2, is     ) ) / deltaB

           rA  = -DOT_PRODUCT( xt( :, is - 1 ), e1xe2( :, is - 1 ) ) / deltaA
           rB  = -DOT_PRODUCT( xt( :, is     ), e1xe2( :, is     ) ) / deltaB

           irA = MAX( MIN( INT( ( rA - Pos%Rr( 1 ) ) / Pos%Delta_r ) + 1, Pos%NRr ), 1 )  ! index of nearest rcvr before normal
           irB = MAX( MIN( INT( ( rB - Pos%Rr( 1 ) ) / Pos%Delta_r ) + 1, Pos%NRr ), 1 )  ! index of nearest rcvr before normal
           ! detect and skip duplicate points (happens at boundary reflection)
           ! IF ( irA /= irB .AND. NORM2( ray3D( is )%x - ray3D( is - 1 )%x ) > 1.0D3 * SPACING( ray3D( is )%x( 1 ) ) ) THEN  ! too slow
           IF ( irA /= irB .AND. NORM2( ray3D( is )%x - ray3D( is - 1 )%x ) > 1.0e-4 ) THEN

              ! *** Compute contributions to bracketted receivers ***
              dq_tilde = ray3D( is )%q_tilde - ray3D( is - 1 )%q_tilde
              dq_hat   = ray3D( is )%q_hat   - ray3D( is - 1 )%q_hat

              II = 0
              IF ( irB <= irA ) II = 1   ! going backwards in range

              Ranges: DO ir = irA + 1 - II, irB + II, SIGN( 1, irB - irA )
                 s = ( Pos%Rr( ir ) - rA ) / ( rB - rA )   

                 ! linear interpolation of q's
                 q_tilde = ray3D( is - 1 )%q_tilde + s * dq_tilde
                 q_hat   = ray3D( is - 1 )%q_hat   + s * dq_hat
                 n       = ABS( nA                 + s * ( nB - nA ) )   ! normal distance to ray
                 m       = ABS( mA                 + s * ( mB - mA ) )   ! normal distance to ray

                 ! represent (m, n) as a linear combination a q + b q
                 DetQint = q_tilde( 1 ) * q_hat( 2 ) - q_hat( 1 ) * q_tilde( 2 )   ! area of parallelogram formed by ray tube
                 IF ( DetQint == 0    ) CYCLE Ranges  ! receiver is outside the beam

                 a = ABS( ( -q_hat(   1 ) * m + q_hat(   2 ) * n ) / DetQint ) !!! sign is flipable inside abs
                 b = ABS( (  q_tilde( 1 ) * m - q_tilde( 2 ) * n ) / DetQint )
                 IF ( MAX( a, b ) > 1 ) CYCLE Ranges  ! receiver is outside the beam

                 W = ( 1 - a ) * ( 1 - b )   ! beamshape is a hat function: 1 on center, 0 on edge
                 delay = ray3D( is - 1 )%tau + s * dtau( is - 1 )

!take this out of the loop!!!
                 
                 const = Ratio1 * ray3D( is )%c / SQRT( ABS( DetQint ) ) * ray3D( is )%Amp
                 If ( Beam%RunType( 1 : 1 ) == 'S' ) const = DS * const   ! semi-coherent TL

                 Amp   = const * W   ! hat function

                 ! phase shift at caustics
                 phaseInt = ray3D( is - 1 )%Phase + KMAHphase( is - 1 )
                 IF ( DetQint <= 0.0d0 .AND. DetQ( is - 1 ) > 0.0d0 .OR. &
                      DetQint >= 0.0d0 .AND. DetQ( is - 1 ) < 0.0d0 ) phaseInt = phaseInt + pi / 2.
                 CALL ApplyContribution( alpha, beta, xs, P( itheta, iz, ir ) )
              END DO Ranges
           END IF
           mA     = mB
           deltaA = deltaB
        END DO Stepping
     END DO Radials
  END DO ReceiverDepths

END SUBROUTINE Influence3DGeoHatRayCen

!**********************************************************************!

SUBROUTINE Influence3DGeoHatCart( xs, alpha, beta, Dalpha, Dbeta, P, x_rcvrMat, t_rcvr )

  ! Geometrically-spreading beams with a hat-shaped beam

  REAL ( KIND=8 ), INTENT( IN  ) :: alpha, beta, Dalpha, Dbeta         ! ray take-off angle
  REAL ( KIND=8 ), INTENT( IN  ) :: xs( 3 )                            ! source coordinate
  REAL ( KIND=8 ), INTENT( IN  ) :: x_rcvrMat( 2, Pos%Ntheta, Pos%NRr ), t_rcvr( 2, Pos%Ntheta ) ! rcvr coordinates and tangent
  COMPLEX        , INTENT( OUT ) :: P( Pos%Ntheta, Pos%Nrz, Pos%NRr )  ! complex pressure
  INTEGER            :: irT( 1 ), irTT
  REAL    ( KIND=8 ) :: s, rlen, DS, x_ray( 3 ), n_ray_z( 3 ), n_ray_theta(3 ), &
       e1( 3 ), e2( 3 ), x_rcvr( 3 ), x_rcvr_ray( 3 ), &
       L_z, L_diag, e_theta( 3 ), m_prime, a, b, zMin, zMax, &
       Det_Q, Det_Qold, DetQint, &
       q_tilde( 2 ), q_hat( 2 ), dq_tilde( 2 ), dq_hat( 2 )
  COMPLEX ( KIND=8 ) :: dtau

  SrcDeclAngle = RadDeg * alpha          ! take-off declination angle in degrees
  SrcAzimAngle = RadDeg * beta           ! take-off azimuthal   angle in degrees

  DS     = SQRT( 2.0 ) * SIN( omega * xs( 3 ) * ray3D( 1 )%t( 3 ) )   ! Lloyd mirror pattern (for semi-coherent source)
  Ratio1 = SQRT( ABS( COS( alpha ) ) ) * SQRT( Dalpha * Dbeta ) / ray3D( 1 )%c

  ! scaling for geometric beams

  ray3D( 1 : Beam%Nsteps )%DetQ     = ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) * ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) - &
                                      ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) * ray3D( 1 : Beam%Nsteps )%q_hat(   1 )

  ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) =                       Dalpha * ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) / ray3D( 1 )%c
  ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) =                       Dalpha * ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) / ray3D( 1 )%c
  ray3D( 1 : Beam%Nsteps )%q_hat(   1 ) = ABS( COS( ALPHA ) ) * Dbeta  * ray3D( 1 : Beam%Nsteps )%q_hat(   1 ) / ray3D( 1 )%c
  ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) = ABS( COS( ALPHA ) ) * Dbeta  * ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) / ray3D( 1 )%c

  phase    = 0.0
  Det_QOld = ray3D( 1 )%DetQ   ! used to track phase changes at caustics (rotations of Det_Q)

  ! Compute nearest rcvr before normal
  rA  = NORM2( ray3D( 1 )%x( 1 : 2 ) - xs( 1 : 2 ) )         ! range of ray point
  irT = MINLOC( Pos%Rr( 1 : Pos%NRr ), MASK = Pos%Rr( 1 : Pos%NRr ) .GT. rA )        ! index of receiver
  ir  = irT( 1 )

  Stepping: DO is = 2, Beam%Nsteps
     ! Compute nearest rcvr before normal
     rB  = NORM2( ray3D( is )%x( 1 : 2 ) - xs( 1 : 2 ) )         ! range of ray point

     IF ( ABS( rB - rA ) > 1.0D3 * SPACING( rA ) ) THEN   ! jump to next step if duplicate point
        ! initialize the index of the receiver range
        IF ( is == 2 ) THEN
           IF ( rB > rA ) THEN   ! ray is moving outward in range
              ir = 1             ! index all the way in
           ELSE                  ! ray is moving inward in range
              ir = Pos%NRr       ! index all the way out
           END IF
        END IF

        x_ray = ray3D( is - 1 )%x

        ! compute normalized tangent (we need it to measure the step length)
        rayt = ray3D( is )%x - ray3D( is - 1 )%x
        rlen = NORM2( rayt )

        IF ( rlen > 1.0D3 * SPACING( ray3D( is )%x( 1 ) ) ) THEN  ! Make sure this is not a duplicate point
           rayt = rayt / rlen                                     ! unit tangent to ray
           CALL RayNormal_unit( rayt, ray3D( is )%phi, e1, e2 )   ! Get ray normals e1 and e2

           ! phase shifts at caustics
           Det_Q  = ray3D( is - 1 )%DetQ
           IF ( Det_Q <= 0.0d0 .AND. Det_QOld > 0.0d0 .OR. Det_Q >= 0.0d0 .AND. Det_QOld < 0.0d0 ) phase = phase + pi / 2.
           Det_Qold = Det_Q

           L1 = MAX( NORM2( ray3D( is - 1 )%q_tilde ), NORM2( ray3D( is )%q_tilde ) ) ! beamwidths
           L2 = MAX( NORM2( ray3D( is - 1 )%q_hat   ), NORM2( ray3D( is )%q_hat   ) )

           L_diag = SQRT( L1 ** 2 + L2 ** 2 )   ! worst case is when rectangle is rotated to catch the hypotenuse

           ! n_ray_theta = CROSS_PRODUCT( rayt, e_z )     ! normal to the ray in the horizontal receiver plane
           n_ray_theta = [ -rayt( 2 ), rayt( 1 ), 0.D0 ]  ! normal to the ray in the horizontal receiver plane

           ! *** Compute contributions to bracketted receivers ***
           dq_tilde = ray3D( is )%q_tilde - ray3D( is - 1 )%q_tilde
           dq_hat   = ray3D( is )%q_hat   - ray3D( is - 1 )%q_hat
           dtau     = ray3D( is )%tau     - ray3D( is - 1 )%tau

           Ranges: DO
              ! is r( ir ) contained in [ rA, rB ]? Then compute beam influence
              IF ( Pos%Rr( ir ) >= MIN( rA, rB ) .AND. Pos%Rr( ir ) < MAX( rA, rB ) ) THEN

                 Radials: DO itheta = 1, Pos%Ntheta   ! Loop over radials of receiver line

                    x_rcvr( 1 : 2 ) = x_rcvrMat( 1 : 2, itheta, ir )
                    m_prime = ABS( DOT_PRODUCT( x_rcvr( 1 : 2 ) - x_ray( 1 : 2 ), n_ray_theta( 1 : 2 ) ) )  ! normal distance from rcvr to ray segment

                    IF ( m_prime > L_diag ) CYCLE Radials

                    ! The set of possible receivers is a ring
                    ! However, extrapolating the beam backwards produces contributions with s negative and large
                    ! We do not want to accept these contributions--- they have the proper range but are 180 degrees
                    ! away from this segement of the ray
                    ! s = DOT_PRODUCT( x_rcvr( 1 : 2 ) - x_ray( 1 : 2 ), rayt( 1 : 2 ) / NORM2( rayt( 1 : 2 ) ) )   ! a distance along ray    (in x-y plane)
                    ! s = DOT_PRODUCT( x_rcvr( 1 : 2 ) - xs( 1 : 2 ), t_rcvr( 1 : 2, itheta ) ) ! a distance along radial (in x-y plane)
                    s = DOT_PRODUCT( x_rcvr( 1 : 2 ) - xs( 1 : 2 ), x_ray( 1 : 2 ) - xs( 1 : 2 ) ) ! vector to rcvr dotted into vector to ray point

                    ! The real receivers have an s-value in [0, R_len]
                    ! IF ( s < 0D0 .OR. s > NORM2( ray3D( is )%x( 1 : 2 ) - ray3D( is - 1 )%x( 1 : 2 ) ) ) THEN
                    ! IF ( ABS( s ) > NORM2( ray3D( is )%x( 1 : 2 ) - ray3D( is - 1 )%x( 1 : 2 ) ) ) THEN

                    IF ( s < 0D0 ) CYCLE Radials

                    ! calculate z-limits for the beam (could be pre-cacluated for each itheta)
                    e_theta      = [ -t_rcvr( 2, itheta ), t_rcvr( 1, itheta ), 0.0D0 ]  ! normal to the vertical receiver plane
                    ! n_ray_z    = CROSS_PRODUCT( rayt, e_theta )                        ! normal to the ray in the vertical receiver plane
                    n_ray_z( 3 ) = rayt( 1 ) * e_theta( 2 ) - rayt( 2 ) * e_theta( 1 )   ! normal to the ray in the vertical receiver plane

                    IF ( ABS( n_ray_z( 3 ) ) < 1D-9 ) CYCLE Radials   ! avoid divide by zero
                    L_z          = L_diag / ABS( n_ray_z( 3 ) )

                    zmin = MIN( ray3D( is - 1 )%x( 3 ), ray3D( is )%x( 3 ) ) - L_z  ! min depth of ray segment
                    zmax = MAX( ray3D( is - 1 )%x( 3 ), ray3D( is )%x( 3 ) ) + L_z  ! max depth of ray segment

                    ReceiverDepths: DO iz = 1, Pos%Nrz
                       x_rcvr( 3 ) = DBLE( Pos%rz( iz ) )   ! z coordinate of the receiver
                       IF ( x_rcvr( 3 ) < zmin .OR. x_rcvr( 3 ) > zmax ) CYCLE ReceiverDepths

                       x_rcvr_ray = x_rcvr - x_ray

                       !!! rlen factor could be built into rayt
                       ! linear interpolation of q's
                       s       = DOT_PRODUCT( x_rcvr_ray, rayt ) / rlen  ! proportional distance along ray
                       q_tilde = ray3D( is - 1 )%q_tilde + s * dq_tilde
                       q_hat   = ray3D( is - 1 )%q_hat   + s * dq_hat

                       L1 = NORM2( q_tilde )  ! beamwidth is length of the vector
                       L2 = NORM2( q_hat   )  ! beamwidth

                       n  = ABS( DOT_PRODUCT( x_rcvr_ray, e1 ) )         ! normal distance to ray
                       m  = ABS( DOT_PRODUCT( x_rcvr_ray, e2 ) )         ! normal distance to ray

                       ! represent (m, n) as a linear combination a q + b q
                       DetQint = q_tilde( 1 ) * q_hat( 2 ) - q_hat( 1 ) * q_tilde( 2 )   ! area of parallelogram formed by ray tube

                       IF ( DetQint == 0.0 ) CYCLE ReceiverDepths
                       a = ABS( ( -q_hat(   1 ) * m + q_hat(   2 ) * n ) / DetQint )
                       b = ABS( (  q_tilde( 1 ) * m - q_tilde( 2 ) * n ) / DetQint )

                       ! suppose qtilde( 2 ) = qhat( 1 ) = 0
                       !DetQint = q_tilde( 1 ) * q_hat( 2 )   ! area of parallelogram formed by ray tube
                       !a = ABS( n / q_tilde( 1 )  )
                       !b = ABS( m / q_hat(   2 )  )

                       IF ( MAX( a, b ) > 1 .OR. L1 == 0 .OR. L2 == 0 ) CYCLE ReceiverDepths   ! receiver is outside the beam

                       delay = ray3D( is - 1 )%tau + s * dtau

                       ! phase shift at caustics
                       phaseInt = ray3D( is - 1 )%Phase + phase
                       IF ( DetQint <= 0.0d0 .AND. Det_QOld > 0.0d0 .OR. &
                            DetQint >= 0.0d0 .AND. Det_QOld < 0.0d0 ) phaseInt = phaseInt + pi / 2.

!take this out of the loop !!!

                       const = Ratio1 * ray3D( is )%c / SQRT( ABS( DetQint ) ) * ray3D( is )%Amp
                       If ( Beam%RunType( 1 : 1 ) == 'S' ) const = DS * const   ! semi-coherent TL
   
                       W     = ( 1 - a ) * ( 1 - b )   ! hat function: 1 on center, 0 on edge
                       Amp   = const * W   ! hat function

                       CALL ApplyContribution( alpha, beta, xs, P( itheta, iz, ir ) )
                    END DO ReceiverDepths
                 END DO Radials

              END IF
              ! bump receiver index, ir, towards rB
              IF ( Pos%Rr( ir ) < rB ) THEN
                 IF ( ir >= Pos%NRr ) EXIT   ! jump out of the search and go to next step on ray
                 irTT = ir + 1               ! bump right
                 IF ( Pos%Rr( irTT ) >= rB ) EXIT
              ELSE
                 IF ( ir <= 1  ) EXIT        ! jump out of the search and go to next step on ray
                 irTT = ir - 1               ! bump left
                 IF ( Pos%Rr( irTT ) <= rB ) EXIT
              END IF
              ir = irTT
           END DO Ranges
        END IF
     END IF
     rA = rB
  END DO Stepping

END SUBROUTINE Influence3DGeoHatCart

!**********************************************************************!

SUBROUTINE Influence3DGeoGaussianRayCen( xs, alpha, beta, Dalpha, Dbeta, P )

  ! Geometrically-spreading beams with a hat-shaped beam in ray-centered coordinates

  REAL, PARAMETER                :: BeamWindow = 4.0                   ! kills beams outside e**(-0.5 * BeamWindow**2 )
  REAL ( KIND=8 ), INTENT( IN  ) :: alpha, beta, Dalpha, Dbeta         ! ray take-off angle
  REAL ( KIND=8 ), INTENT( IN  ) :: xs( 3 )                            ! source coordinate
  COMPLEX        , INTENT( OUT ) :: P( Pos%Ntheta, Pos%Nrz, Pos%NRr )  ! complex pressure
  INTEGER          :: irA, irB, II
  REAL    (KIND=8) :: nA, nB, mA, mB, DS, deltaA, deltaB, t_rcvr( 2, Pos%Ntheta ), &
       L1_stent, L2_stent, KMAHphase( Beam%Nsteps ), DetQint
  REAL    (KIND=8) :: e1G( 3, Beam%Nsteps ), e2G( 3, Beam%Nsteps ), e1xe2( 3, Beam%Nsteps ), &
                      xt( 3, Beam%Nsteps ),  xtxe1( 3, Beam%Nsteps ),  xtxe2( 3, Beam%Nsteps )
  REAL    (KIND=8) :: q_tilde( Beam%Nsteps ), q_hat( Beam%Nsteps ),  q_tilde_hat( Beam%Nsteps ), DetQ( Beam%Nsteps ), &
       dq_tilde( Beam%Nsteps - 1 ), dq_hat( Beam%Nsteps - 1 ), &
       MaxRadius_m( Beam%Nsteps - 1 ), MaxRadius_n( Beam%Nsteps - 1 )
  COMPLEX (KIND=8) :: dtau( Beam%Nsteps - 1 )

  SrcDeclAngle = RadDeg * alpha          ! take-off declination angle in degrees
  SrcAzimAngle = RadDeg * beta           ! take-off azimuthal   angle in degrees

  DS     = SQRT( 2.0 ) * SIN( omega * xs( 3 ) * ray3D( 1 )%t( 3 ) )   ! Lloyd mirror pattern (for semi-coherent source)
  Ratio1 = SQRT( ABS( COS( alpha ) ) ) / ( 2. * pi )

  DetQ     = ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) * ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) - &
             ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) * ray3D( 1 : Beam%Nsteps )%q_hat(   1 )
  q_tilde  =                       Dalpha * ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) / ray3D( 1 )%c
  q_hat    = ABS( COS( ALPHA ) ) * Dbeta  * ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) / ray3D( 1 )%c

  q_tilde_hat  =                       Dalpha * ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) / ray3D( 1 )%c * &
                 ABS( COS( ALPHA ) ) * Dbeta  * ray3D( 1 : Beam%Nsteps )%q_hat(   1 ) / ray3D( 1 )%c

  dq_tilde = q_tilde( 2 : Beam%Nsteps )     - q_tilde( 1 : Beam%Nsteps - 1 )
  dq_hat   = q_hat(   2 : Beam%Nsteps )     - q_hat(   1 : Beam%Nsteps - 1 )
  dtau     = ray3D(   2 : Beam%Nsteps )%tau - ray3D(   1 : Beam%Nsteps - 1 )%tau

  ! phase shifts at caustics (rotations of Det_Q)
  KMAHphase( 1 ) = 0
  DO is = 2, Beam%Nsteps
     KMAHphase( is ) = KMAHphase( is - 1 )
     IF ( DetQ( is ) <= 0.0d0 .AND. DetQ( is - 1 ) > 0.0d0 .OR. &
          DetQ( is ) >= 0.0d0 .AND. DetQ( is - 1 ) < 0.0d0 ) KMAHphase( is ) = KMAHphase( is - 1 ) + pi / 2.
  END DO

  ! pre-calculate tangents, normals (e1, e2), etc.
  DO is = 1, Beam%Nsteps
     xt( :, is ) = ray3D( is )%x - xs   ! vector from the origin of the receiver plane to this point on the ray
     CALL RayNormal( ray3D( is )%t, ray3D( is )%phi, ray3D( is )%c, e1G( :, is ), e2G( :, is ) )

     ! e1xe2( :, is ) = cross_product( e1G( :, is ), e2G( :, is ) )
     e1xe2( :, is ) = ray3D( is )%c * ray3D( is )%t   ! unit tangent to ray
  END DO

  DO is = 1, Beam%Nsteps - 1
     MaxRadius_m( is ) = BeamWindow * MAX( ABS( q_hat(   is ) ), ABS( q_hat(   is + 1 ) ) )
     MaxRadius_n( is ) = BeamWindow * MAX( ABS( q_tilde( is ) ), ABS( q_tilde( is + 1 ) ) )
  END DO
  
  ! tangent along receiver bearing line
  t_rcvr( 1, : ) = COS( DegRad * Pos%theta( 1 : Pos%Ntheta ) )
  t_rcvr( 2, : ) = SIN( DegRad * Pos%theta( 1 : Pos%Ntheta ) )

  ReceiverDepths: DO iz = 1, Pos%Nrz
     ! precalculate vector from receiver to each step of ray
     xt( 3, : ) = ray3D( 1 : Beam%Nsteps )%x( 3 ) - Pos%rz( iz )
     DO is = 1, Beam%Nsteps
        xtxe1( :, is ) = cross_product( xt( :, is ), e1G( :, is ) )
        xtxe2( :, is ) = cross_product( xt( :, is ), e2G( :, is ) )
     END DO

     Radials: DO itheta = 1, Pos%Ntheta
        ! *** Compute coordinates of intercept: nA, mA, rA ***
        is = 1
        deltaA = -DOT_PRODUCT( t_rcvr( :, itheta ), e1xe2( 1 : 2, is ) )

        ! Check for ray normal || radial of rcvr line
        IF ( ABS( deltaA ) < 1D3 * SPACING( deltaA ) ) THEN
           irA = 0   ! serves as a flag that this normal can't be used
        ELSE
           mA  =  DOT_PRODUCT( t_rcvr( :, itheta ), xtxe1( 1 : 2, is ) ) / deltaA
        END IF

        ! step along the beam ...
        ! Most of the time the beam makes no contribution to a receiver
        ! Therefore we try to test that quickly and move on to the next receiver

        Stepping: DO is = 2, Beam%Nsteps
           lambda = ray3D( is - 1 )%c / freq   ! local wavelength

           ! *** Compute coordinates of intercept: nB, mB, rB ***
           deltaA = -DOT_PRODUCT( t_rcvr( :, itheta ), e1xe2( 1 : 2, is - 1 ) )
           deltaB = -DOT_PRODUCT( t_rcvr( :, itheta ), e1xe2( 1 : 2, is     ) )

           ! Check for ray normal || radial of rcvr line
           IF ( ABS( deltaA ) < 1e-3 .OR. ABS( deltaB ) < 1e-3 )  THEN
              CYCLE Stepping
           END IF

           mA  =  DOT_PRODUCT( t_rcvr( :, itheta ), xtxe1( 1 : 2, is - 1 ) ) / deltaA
           mB  =  DOT_PRODUCT( t_rcvr( :, itheta ), xtxe1( 1 : 2, is     ) ) / deltaB
!!! stent needs to be applied here already
           ! Possible contribution if max possible beamwidth > min possible distance to receiver
           IF ( MaxRadius_m( is - 1 ) > MIN( ABS( mA ), ABS( mB ) ) .OR. ( mA * mB < 0 ) ) THEN
              nA = -DOT_PRODUCT( t_rcvr( :, itheta ), xtxe2( 1 : 2, is - 1 ) ) / deltaA
              nB = -DOT_PRODUCT( t_rcvr( :, itheta ), xtxe2( 1 : 2, is     ) ) / deltaB
              
              ! Possible contribution if max possible beamwidth > min possible distance to receiver
              IF ( MaxRadius_n( is - 1 ) > MIN( ABS( nA ), ABS( nB ) ) .OR. ( nA * nB < 0 ) ) THEN
                 !!! can generate an inexact result if the resulting integer is too big
                 !!! assumes uniform spacing in Pos%Rr
                 rA  = -DOT_PRODUCT( xt( :, is - 1 ), e1xe2( :, is - 1 ) ) / deltaA
                 rB  = -DOT_PRODUCT( xt( :, is     ), e1xe2( :, is     ) ) / deltaB

                 irA = MAX( MIN( INT( ( rA - Pos%Rr( 1 ) ) / Pos%Delta_r ) + 1, Pos%NRr ), 1 )  ! index of nearest rcvr before normal
                 irB = MAX( MIN( INT( ( rB - Pos%Rr( 1 ) ) / Pos%Delta_r ) + 1, Pos%NRr ), 1 )  ! index of nearest rcvr before normal

                 ! detect and skip duplicate points (happens at boundary reflection)
                 IF ( irA /= irB .AND. NORM2( ray3D( is )%x - ray3D( is - 1 )%x ) > 1.0D3 * SPACING( ray3D( is )%x( 1 ) ) ) THEN

                    ! *** Compute contributions to bracketted receivers ***
                    
                    II = 0
                    IF ( irB <= irA ) II = 1   ! going backwards in range

                    Ranges: DO ir = irA + 1 - II, irB + II, SIGN( 1, irB - irA )
                       W = ( Pos%Rr( ir ) - rA ) / ( rB - rA )

                       ! Within beam window?
                       n  = ABS( nA + W * ( nB - nA ) )
                       L1 = ABS( q_tilde( is - 1 ) + W * dq_tilde( is - 1 ) )   ! beamwidth
                       IF ( n > BeamWindow * L1 ) CYCLE Ranges                  ! in beamwindow?

                       m  = ABS( mA + W * ( mB - mA ) )
                       L2 = ABS( q_hat(   is - 1 ) + W * dq_hat(   is - 1 ) )   ! beamwidth
                       IF ( m > BeamWindow * L2 ) CYCLE Ranges                  ! in beamwwindow?

                       ! calculate the beamwidth (must be at least lambda, except in the nearfield)
                       ! comment out the following 2 lines to turn the stent off
                       L1_stent = MAX( L1, MIN( 0.2 * freq * REAL( ray3D( is - 1 )%tau ), pi * lambda ) )
                       L2_stent = MAX( L2, MIN( 0.2 * freq * REAL( ray3D( is - 1 )%tau ), pi * lambda ) )

                       DetQint = DetQ(   is - 1 )     + W * ( DetQ( is ) - DetQ( is - 1 )    )
                       delay   = ray3D(  is - 1 )%tau + W * dtau( is - 1 )

                       ! phase shift at caustics
                       phaseInt = ray3D( is - 1 )%Phase + KMAHphase( is - 1 )
                       IF ( DetQint <= 0.0d0 .AND. DetQ( is - 1 ) > 0.0d0 .OR. &
                            DetQint >= 0.0d0 .AND. DetQ( is - 1 ) < 0.0d0 ) phaseInt = PhaseInt + pi / 2.

                       DetQint = L1 * L2 * ray3D( 1 )%c ** 2 / ( Dalpha * Dbeta )  ! based on actual beamwidth

!take this out of the loop
                       const   = Ratio1 * ray3D( is )%c / SQRT( ABS( DetQint ) ) * ray3D( is )%Amp
                       IF ( Beam%RunType( 1 : 1 ) == 'S' ) const = DS * const      ! semi-coherent TL

                       ! Beam shape function (Gaussian)
                       ! The factor involving L1, L2 compensates for the change in intensity when widening the beam
                       W   = EXP( -.5 * ( ( n / L1_stent ) ** 2 + ( m / L2_stent ) ** 2 ) ) * L1 / L1_stent * L2 / L2_stent
                       Amp = const * W

                       CALL ApplyContribution( alpha, beta, xs, P( itheta, iz, ir ) )
                    END DO Ranges
                 END IF
              END IF
           END IF
           mA     = mB
           deltaA = deltaB
        END DO Stepping
     END DO Radials
  END DO ReceiverDepths

END SUBROUTINE Influence3DGeoGaussianRayCen


!**********************************************************************!

SUBROUTINE Influence3DGeoGaussianCart( xs, alpha, beta, Dalpha, Dbeta, P, x_rcvrMat, t_rcvr )

  ! Geometrically-spreading beams with a Gaussian-shaped beam

  REAL, PARAMETER                :: BeamWindow = 4.0                   ! kills beams outside e**(-0.5 * BeamWindow**2 )
  REAL ( KIND=8 ), INTENT( IN  ) :: alpha, beta, Dalpha, Dbeta         ! ray take-off angle
  REAL ( KIND=8 ), INTENT( IN  ) :: xs( 3 )                            ! source coordinate
  REAL ( KIND=8 ), INTENT( IN  ) :: x_rcvrMat( 2, Pos%Ntheta, Pos%NRr ), t_rcvr( 2, Pos%Ntheta ) ! rcvr coordinates and tangent
  COMPLEX        , INTENT( OUT ) :: P( Pos%Ntheta, Pos%Nrz, Pos%NRr )  ! complex pressure
  INTEGER            :: irT( 1 ), irTT
  REAL    ( KIND=8 ) :: rlen, a, b, DS, x_ray( 3 ), n_ray_z( 3 ), n_ray_theta( 3 ), &
       e1( 3 ), e2( 3 ), x_rcvr( 3 ), x_rcvr_ray( 3 ), &
       L1_stent, L2_stent, L_z, L_diag, e_theta( 3 ), m_prime, zMin, zMax, &
       Det_Q, Det_Qold, DetQint, &
       q_tilde( 2 ), q_hat( 2 ), dq_tilde( 2 ), dq_hat( 2 )
  COMPLEX ( KIND=8 ) :: dtau

  SrcDeclAngle = RadDeg * alpha          ! take-off declination angle in degrees
  SrcAzimAngle = RadDeg * beta           ! take-off azimuthal   angle in degrees

  DS     = SQRT( 2.0 ) * SIN( omega * xs( 3 ) * ray3D( 1 )%t( 3 ) )   ! Lloyd mirror pattern (for semi-coherent source)
  Ratio1 = SQRT( ABS( COS( alpha ) ) ) * SQRT( Dalpha * Dbeta ) / ray3D( 1 )%c / ( 2.0 * pi )

  ! scaling for geometric beams

  ray3D( 1 : Beam%Nsteps )%DetQ     = ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) * ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) - &
                                      ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) * ray3D( 1 : Beam%Nsteps )%q_hat(   1 )

  !ray3D%q_tilde =                       Dalpha * ray3D%q_tilde / ray3D( 1 )%c
  !ray3D%q_hat   = ABS( COS( ALPHA ) ) * Dbeta  * ray3D%q_hat   / ray3D( 1 )%c

  ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) =                       Dalpha * ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) / ray3D( 1 )%c
  ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) =                       Dalpha * ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) / ray3D( 1 )%c
  ray3D( 1 : Beam%Nsteps )%q_hat(   1 ) = ABS( COS( ALPHA ) ) * Dbeta  * ray3D( 1 : Beam%Nsteps )%q_hat(   1 ) / ray3D( 1 )%c
  ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) = ABS( COS( ALPHA ) ) * Dbeta  * ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) / ray3D( 1 )%c

  phase    = 0.0
  Det_QOld = ray3D( 1 )%DetQ   ! used to track phase changes at caustics (rotations of Det_Q)

  ! Compute nearest rcvr before normal
  rA  = NORM2( ray3D( 1 )%x( 1 : 2 ) - xs( 1 : 2 ) )         ! range of ray point
  irT = MINLOC( Pos%Rr( 1 : Pos%NRr ), MASK = Pos%Rr( 1 : Pos%NRr ) .GT. rA )   ! index of receiver
  ir  = irT( 1 )

  Stepping: DO is = 2, Beam%Nsteps
     lambda = ray3D( is - 1 )%c / freq   ! local wavelength
     ! Compute nearest rcvr before normal
     rB  = NORM2( ray3D( is )%x( 1 : 2 ) - xs( 1 : 2 ) )         ! range of ray point

     IF ( ABS( rB - rA ) > 1.0D3 * SPACING( rA ) ) THEN   ! jump to next step if duplicate point
        ! initialize the index of the receiver range
        IF ( is == 2 ) THEN
           IF ( rB > rA ) THEN   ! ray is moving outward in range
              ir = 1             ! index all the way in
           ELSE                  ! ray is moving inward in range
              ir = Pos%NRr       ! index all the way out
           END IF
        END IF

        x_ray = ray3D( is - 1 )%x

        ! compute normalized tangent (we need it to measure the step length)
        rayt = ray3D( is )%x - ray3D( is - 1 )%x
        rlen = NORM2( rayt )

        IF ( rlen > 1.0D3 * SPACING( ray3D( is )%x( 1 ) ) ) THEN  ! Make sure this is not a duplicate point
           rayt = rayt / rlen                                     ! unit tangent to ray
           CALL RayNormal_unit( rayt, ray3D( is )%phi, e1, e2 )   ! Get ray normals e1 and e2

           ! phase shifts at caustics
           Det_Q  = ray3D( is - 1 )%DetQ
           IF ( Det_Q <= 0.0d0 .AND. Det_QOld > 0.0d0 .OR. Det_Q >= 0.0d0 .AND. Det_QOld < 0.0d0 ) phase = phase + pi / 2.
           Det_Qold = Det_Q

           L1 = MAX( NORM2( ray3D( is - 1 )%q_tilde ), NORM2( ray3D( is )%q_tilde ) ) ! beamwidths
           L2 = MAX( NORM2( ray3D( is - 1 )%q_hat   ), NORM2( ray3D( is )%q_hat   ) )

           L_diag = SQRT( L1 ** 2 + L2 ** 2 )   ! worst case is when rectangle is rotated to catch the hypotenuse

           ! n_ray_theta = CROSS_PRODUCT( rayt, e_z )     ! normal to the ray in the horizontal receiver plane
           n_ray_theta = [ -rayt( 2 ), rayt( 1 ), 0.D0 ]  ! normal to the ray in the horizontal receiver plane

           ! *** Compute contributions to bracketted receivers ***
           dq_tilde = ray3D( is )%q_tilde - ray3D( is - 1 )%q_tilde
           dq_hat   = ray3D( is )%q_hat   - ray3D( is - 1 )%q_hat
           dtau     = ray3D( is )%tau     - ray3D( is - 1 )%tau

           Ranges: DO
              ! is r( ir ) contained in [ rA, rB ]? Then compute beam influence
              IF ( Pos%Rr( ir ) >= MIN( rA, rB ) .AND. Pos%Rr( ir ) < MAX( rA, rB ) ) THEN

                 Radials: DO itheta = 1, Pos%Ntheta   ! Loop over radials of receiver line

                    x_rcvr( 1 : 2 ) = x_rcvrMat( 1 : 2, itheta, ir )
                    m_prime = ABS( DOT_PRODUCT( x_rcvr( 1 : 2 ) - x_ray( 1 : 2 ), n_ray_theta( 1 : 2 ) ) )  ! normal distance from rcvr to ray segment

                    !!!!!!!IF ( m_prime > BeamWindow * L_diag ) CYCLE Radials

                    ! The set of possible receivers is a ring
                    ! However, extrapolating the beam backwards produces contributions with s negative and large
                    ! We do not want to accept these contributions--- they have the proper range but are 180 degrees
                    ! away from this segement of the ray
                    !!! pre-calculate unit ray tangent
                    ! s = DOT_PRODUCT( x_rcvr( 1 : 2 ) - x_ray( 1 : 2 ), rayt( 1 : 2 ) / NORM2( rayt( 1 : 2 ) ) )   ! a distance along ray    (in x-y plane)
                    ! s = DOT_PRODUCT( x_rcvr( 1 : 2 ) - xs( 1 : 2 ), t_rcvr( 1 : 2, itheta ) ) ! a distance along radial (in x-y plane)
                    s = DOT_PRODUCT( x_rcvr( 1 : 2 ) - xs( 1 : 2 ), x_ray( 1 : 2 ) - xs( 1 : 2 ) ) ! vector to rcvr dotted into vector to ray point

                    ! The real receivers have an s-value in [0, R_len]
                    ! IF ( s < 0D0 .OR. s > NORM2( ray3D( is )%x( 1 : 2 ) - ray3D( is - 1 )%x( 1 : 2 ) ) ) THEN
                    ! IF ( ABS( s ) > NORM2( ray3D( is )%x( 1 : 2 ) - ray3D( is - 1 )%x( 1 : 2 ) ) ) THEN

                    IF ( s < 0D0 ) CYCLE Radials

                    ! calculate z-limits for the beam (could be pre-cacluated for each itheta)
                    e_theta      = [ -t_rcvr( 2, itheta ), t_rcvr( 1, itheta ), 0.0D0 ]  ! normal to the vertical receiver plane
                    ! n_ray_z    = CROSS_PRODUCT( rayt, e_theta )                        ! normal to the ray in the vertical receiver plane
                    n_ray_z( 3 ) = rayt( 1 ) * e_theta( 2 ) - rayt( 2 ) * e_theta( 1 )   ! normal to the ray in the vertical receiver plane

                    IF ( ABS( n_ray_z( 3 ) ) < 1D-9 ) CYCLE Radials   ! avoid divide by zero
                    L_z          = BeamWindow * L_diag / ABS( n_ray_z( 3 ) )

                    zmin = MIN( ray3D( is - 1 )%x( 3 ), ray3D( is )%x( 3 ) ) - L_z  ! min depth of ray segment
                    zmax = MAX( ray3D( is - 1 )%x( 3 ), ray3D( is )%x( 3 ) ) + L_z  ! max depth of ray segment

                    ReceiverDepths: DO iz = 1, Pos%Nrz
                       x_rcvr( 3 ) = DBLE( Pos%rz( iz ) )   ! z coordinate of the receiver
                       IF ( x_rcvr( 3 ) < zmin .OR. x_rcvr( 3 ) > zmax ) CYCLE ReceiverDepths

                       x_rcvr_ray = x_rcvr - x_ray

                       !!! rlen factor could be built into rayt
                       ! linear interpolation of q's
                       s       = DOT_PRODUCT( x_rcvr_ray, rayt ) / rlen  ! proportional distance along ray
                       q_tilde = ray3D( is - 1 )%q_tilde + s * dq_tilde
                       q_hat   = ray3D( is - 1 )%q_hat   + s * dq_hat

                       L1 = NORM2( q_tilde )  ! beamwidth is length of the vector
                       L2 = NORM2( q_hat   )  ! beamwidth

                       n  = ABS( DOT_PRODUCT( x_rcvr_ray, e1 ) )         ! normal distance to ray
                       m  = ABS( DOT_PRODUCT( x_rcvr_ray, e2 ) )         ! normal distance to ray

                       ! represent (m, n) as a linear combination a q + b q
                       DetQint = q_tilde( 1 ) * q_hat( 2 ) - q_hat( 1 ) * q_tilde( 2 )   ! area of parallelogram formed by ray tube
                       IF ( DetQint == 0.0 ) CYCLE ReceiverDepths

                       a = ABS( ( -q_hat(   1 ) * m + q_hat(   2 ) * n ) / DetQint )
                       b = ABS( (  q_tilde( 1 ) * m - q_tilde( 2 ) * n ) / DetQint )

                       IF ( a + b > BeamWindow .OR. L1 == 0 .OR. L2 == 0 ) CYCLE ReceiverDepths   ! receiver is outside the beam

                       delay = ray3D( is - 1 )%tau + s * dtau

                       ! phase shift at caustics
                       phaseInt = ray3D( is - 1 )%Phase + phase
                       IF ( DetQint <= 0.0d0 .AND. Det_QOld > 0.0d0 .OR. &

                            DetQint >= 0.0d0 .AND. Det_QOld < 0.0d0 ) phaseInt = phaseInt + pi / 2.
!take this out of the loop
                       const = Ratio1 * ray3D( is )%c / SQRT( ABS( DetQint ) ) * ray3D( is )%Amp
                       IF ( Beam%RunType( 1 : 1 ) == 'S' ) const = DS * const   ! semi-coherent TL

                       ! calculate the beamwidth (must be at least lambda, except in the nearfield)
                       ! comment out the following 2 lines to turn the stent off
                       !L1_stent  = MAX( L1, MIN( 0.2 * freq * REAL( ray3D( is )%tau ), pi * lambda ) )
                       !L2_stent  = MAX( L2, MIN( 0.2 * freq * REAL( ray3D( is )%tau ), pi * lambda ) )

                       ! Beam shape function (Gaussian)
                       ! The factor involving L1, L2 compensates for the change in intensity when widening the beam

                       !W   = EXP( -.5 * ( n ** 2 + m ** 2 ) ) * L1 / L1_stent * L2 / L2_stent
                       W   = EXP( -.5 * ( a ** 2 + b ** 2 ) )   ! Gaussian
                       Amp = const * W

                       CALL ApplyContribution( alpha, beta, xs, P( itheta, iz, ir ) )
                    END DO ReceiverDepths
                 END DO Radials

              END IF
              ! bump receiver index, ir, towards rB
              IF ( Pos%Rr( ir ) < rB ) THEN
                 IF ( ir >= Pos%NRr ) EXIT   ! jump out of the search and go to next step on ray
                 irTT = ir + 1               ! bump right
                 IF ( Pos%Rr( irTT ) >= rB ) EXIT
              ELSE
                 IF ( ir <= 1  ) EXIT        ! jump out of the search and go to next step on ray
                 irTT = ir - 1               ! bump left
                 IF ( Pos%Rr( irTT ) <= rB ) EXIT
              END IF
              ir = irTT
           END DO Ranges
        END IF
     END IF
     rA = rB
  END DO Stepping

END SUBROUTINE Influence3DGeoGaussianCart

! **********************************************************************!
  
SUBROUTINE ApplyContribution( alpha, beta, xs, U )

  REAL ( KIND=8 ), INTENT( IN    ) :: alpha, beta         ! ray take-off angle
  REAL ( KIND=8 ), INTENT( IN    ) :: xs( 3 )             ! source coordinate
  COMPLEX,         INTENT( INOUT ) :: U
  
  SELECT CASE( Beam%RunType( 1 : 1 ) )
  CASE ( 'E' )      ! eigenrays
     CALL WriteRay3D( alpha, beta, is, xs )   ! produces no output if NR=1
  CASE ( 'A', 'a' ) ! arrivals
     rayt = ray3D( is )%x - ray3D( is - 1 )%x ! ray tangent !!! does this always need to be done???
     RcvrDeclAngle = RadDeg * ATAN2( rayt( 3 ), NORM2( rayt( 1 : 2 ) ) )
     RcvrAzimAngle = RadDeg * ATAN2( rayt( 2 ), rayt( 1 ) )

     CALL AddArr3D( omega, itheta, iz, ir, Amp, phaseInt, delay, &
          SrcDeclAngle, SrcAzimAngle, RcvrDeclAngle, RcvrAzimAngle, &
          ray3D( is )%NumTopBnc, ray3D( is )%NumBotBnc )
  CASE ( 'C'  )     ! coherent TL
     U = U + CMPLX( Amp * EXP( -i * ( omega * delay - phaseInt ) ) )
  CASE DEFAULT      ! incoherent/semi-coherent TL
      IF ( Beam%Type( 1 : 1 ) == 'B' ) THEN   ! Gaussian beam
          U = U + SNGL( 2. * pi * ( const * EXP( AIMAG( omega * delay ) ) ) ** 2 * W )
       ELSE
          U = U + SNGL(           ( const * EXP( AIMAG( omega * delay ) ) ) ** 2 * W )
       END IF
  END SELECT

END SUBROUTINE ApplyContribution

! **********************************************************************!

SUBROUTINE ScalePressure3D( Dalpha, Dbeta, c, epsilon, P, Ntheta, Nrz, Nr, RunType, freq )

  ! Scale the pressure field

  IMPLICIT NONE
  INTEGER,            INTENT( IN    ) :: Ntheta, Nrz, Nr
  REAL    ( KIND=8 ), INTENT( IN    ) :: Dalpha, Dbeta        ! angular spacing between rays
  REAL    ( KIND=8 ), INTENT( IN    ) :: freq, c              ! source frequency, nominal sound speed
  COMPLEX,            INTENT( INOUT ) :: P( Ntheta, Nrz, Nr ) ! Pressure field
  COMPLEX ( KIND=8 ), INTENT( IN    ) :: epsilon( 2 )
  CHARACTER (LEN=5 ), INTENT( IN    ) :: RunType
  COMPLEX ( KIND=8 )                  :: const

  ! Compute scale factor for field
  SELECT CASE ( RunType( 2 : 2 ) )
  CASE ( 'C' )   ! Cerveny Gaussian beams in Cartesian coordinates
     ! epsilon is normally imaginary here, so const is complex
     const = SQRT( epsilon( 1 ) * epsilon( 2 ) ) * freq * Dbeta * Dalpha / ( SQRT( c ) ) **3
     P( :, :, : ) = CMPLX( const, KIND=4 ) * P( :, :, : )
  CASE DEFAULT
     const = 1.0
  END SELECT

  IF ( RunType( 1 : 1 ) /= 'C' ) P = SQRT( REAL( P ) ) ! For incoherent run, convert intensity to pressure

END SUBROUTINE ScalePressure3D
END MODULE Influence3D
