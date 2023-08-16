MODULE Influence

  ! Compute the beam influence, i.e. the contribution of a single beam to the complex pressure
  ! mbp 12/2018, based on much older subroutines

  USE bellhopMod
  USE SourceReceiverPositions
  USE ArrMod
  USE sspMod   ! used to construct image beams in the Cerveny style beam routines
  USE WriteRay

  IMPLICIT NONE
  INTEGER,          PRIVATE :: iz, ir, iS
  REAL    (KIND=8), PRIVATE :: Ratio1 = 1.0D0   ! scale factor for a line source
  REAL    (KIND=8), PRIVATE :: W, s, n, Amp, phase, const, phaseInt, q0, q, qold, RcvrDeclAngle, rA, rB
  COMPLEX (KIND=8), PRIVATE :: delay

CONTAINS
  SUBROUTINE InfluenceCervenyRayCen( U, eps, alpha, iBeamWindow2, RadiusMax )

    ! Paraxial (Cerveny-style) beams in ray-centered coordinates

    INTEGER,          INTENT( IN    ) :: IBeamWindow2
    REAL    (KIND=8), INTENT( IN    ) :: alpha, RadiusMax             ! take-off angle
    COMPLEX,          INTENT( INOUT ) :: U( NRz_per_range, Pos%NRr )  ! complex pressure field
    COMPLEX (KIND=8), INTENT( IN    ) :: eps
    INTEGER          :: ir1, ir2, KMAHV( MaxN ), KMAH, image
    REAL    (KIND=8) :: nA, nB, nSq, c, zr
    REAL    (KIND=8) :: znV( Beam%Nsteps ), rnV( Beam%Nsteps )   ! ray normal
    COMPLEX (KIND=8) :: pVB( MaxN ), qVB( MaxN ), q, epsV( MaxN ), contri, gammaV( MaxN ), gamma, P_n, P_s
    COMPLEX (KIND=8) :: tau

!!! need to add logic related to NRz_per_range

    ! During reflection imag( q ) is constant and adjacent normals cannot bracket a segment of the TL
    ! line, so no special treatment is necessary

    IF ( Beam%Type( 2 : 2 ) == 'C' ) THEN
       epsV( 1 : Beam%Nsteps ) = i * ABS( ray2D( 1 : Beam%Nsteps )%q( 1 ) / ray2D( 1 : Beam%Nsteps )%q( 2 ) )
    ELSE
       epsV( 1 : Beam%Nsteps ) = eps
    END IF

    pVB(    1 : Beam%Nsteps ) = ray2D( 1 : Beam%Nsteps )%p( 1 ) + epsV( 1 : Beam%Nsteps ) * ray2D( 1 : Beam%Nsteps )%p( 2 )
    qVB(    1 : Beam%Nsteps ) = ray2D( 1 : Beam%Nsteps )%q( 1 ) + epsV( 1 : Beam%Nsteps ) * ray2D( 1 : Beam%Nsteps )%q( 2 )
    gammaV( 1 : Beam%Nsteps ) = pVB(   1 : Beam%Nsteps ) / qVB( 1 : Beam%Nsteps )

    ! pre-calculate ray normal based on tangent with c(s) scaling
    znV = -ray2D( 1 : Beam%Nsteps )%t( 1 ) * ray2D( 1 : Beam%Nsteps )%c
    rnV =  ray2D( 1 : Beam%Nsteps )%t( 2 ) * ray2D( 1 : Beam%Nsteps )%c

    IF ( Beam%RunType( 4 : 4 ) == 'R' ) Ratio1 = SQRT( ABS( COS( alpha ) ) )  ! point source

    ! compute KMAH index
    ! Following is incorrect for 'Cerveny'-style beamwidth (narrow as possible)
    KMAHV(  1 ) = 1

    DO iS = 2, Beam%Nsteps
       KMAHV(  iS ) = KMAHV( iS - 1 )
       CALL BranchCut( qVB( iS - 1 ), qVB( iS ), Beam%Type, KMAHV( iS ) )
    END DO

    RcvrDepths: DO iz = 1, NRz_per_range
       zR = Pos%Rz( iz )

       Images: DO image = 1, Beam%Nimage

          !!! This logic means that the first step along the ray is skipped
          !!! which is a problem if deltas is very large, e.g. isospeed problems
          !!! I fixed this in InfluenceGeoHatRayCen
          ir1 = HUGE( ir1 )

          Stepping: DO iS = 2, Beam%Nsteps

             ! Compute ray-centered coordinates, (znV, rnV)

             IF ( ABS( znV( iS ) ) < EPSILON( znV( iS ) ) ) CYCLE Stepping ! If normal parallel to TL-line, skip to next step on ray

             SELECT CASE ( image )     ! Images of beams
             CASE ( 1 )                ! True beam
                nB  = ( zR -                             ray2D( iS )%x( 2 )   ) / znV( iS )
             CASE ( 2 )                ! Surface-reflected beam
                rnV = -rnV
                nB  = ( zR - ( 2.0 * Bdry%Top%HS%Depth - ray2D( iS )%x( 2 ) ) ) / znV( iS )
             CASE ( 3 )                ! Bottom-reflected beam
                rnV = -rnV
                nB  = ( zR - ( 2.0 * Bdry%Bot%HS%Depth - ray2D( iS )%x( 2 ) ) ) / znV( iS )
             END SELECT

             rB  = ray2D( iS )%x( 1 ) + nB * rnV( iS )
             !!! following assumes uniform space in Pos%r
             ir2 = MAX( MIN( INT( ( rB - Pos%Rr( 1 ) ) / Pos%Delta_r ) + 1, Pos%NRr ), 1 ) ! index of receiver

             ! detect and skip duplicate points (happens at boundary reflection)
             IF ( ir1 >= ir2 .OR. ABS( ray2D( iS )%x( 1 ) - ray2D( iS - 1 )%x( 1 ) ) < 1.0D3 * SPACING( ray2D( iS )%x( 1 ) ) ) THEN
                rA  = rB
                nA  = nB
                ir1 = ir2
                CYCLE Stepping
             END IF

             RcvrRanges: DO ir = ir1 + 1, ir2    ! Compute influence for each rcvr
                W     = ( Pos%Rr( ir ) - rA ) / ( rB - rA )
                q     =    qVB( iS - 1 ) + W * (    qVB( iS ) -    qVB( iS - 1 ) )
                gamma = gammaV( iS - 1 ) + W * ( gammaV( iS ) - gammaV( iS - 1 ) )
                n     = nA + W * ( nB - nA )
                nSq   = n * n
                IF ( AIMAG( gamma ) > 0 ) THEN
                   WRITE( PRTFile, * ) 'Unbounded beam'
                   CYCLE RcvrRanges
                END IF

                IF ( -0.5 * omega * AIMAG( gamma ) * nSq < iBeamWindow2 ) THEN   ! Within beam window?
                   c      = ray2D( iS - 1 )%c
                   tau    = ray2D( iS - 1 )%tau + W * ( ray2D( iS )%tau - ray2D( iS - 1 )%tau )
                   contri = ratio1 * ray2D( iS )%Amp * SQRT( c * ABS( epsV( iS ) ) / q ) * &
                        EXP( -i * ( omega * ( tau + 0.5 * gamma * nSq ) - ray2D( iS )%phase ) )

                   SELECT CASE ( Beam%Component )
                   CASE ( 'P' )   ! pressure
                   CASE ( 'V' )   ! vertical component
                      P_n    = -i * omega * gamma * n * contri
                      P_s    = -i * omega / c         * contri
                      contri = c * DOT_PRODUCT( [ P_n, P_s ], ray2D( iS )%t ) 
                   CASE ( 'H' )   ! horizontal component
                      P_n    = -i * omega * gamma * n * contri
                      P_s    = -i * omega / c         * contri
                      contri = c * ( -P_n * ray2D( iS )%t( 2 ) + P_s * ray2D( iS )%t( 1 ) ) 
                   END SELECT

                   KMAH = KMAHV( iS - 1 )
                   CALL BranchCut( qVB( iS - 1 ), q, Beam%Type, KMAH ) ! Get correct branch of SQRT

                   IF ( KMAH  < 0  ) contri = -contri
                   IF ( image == 2 ) contri = -contri

                   SELECT CASE ( Beam%RunType( 1 : 1 ) )
                   CASE ( 'I', 'S' )   ! Incoherent or Semi-coherent TL
                      contri = ABS( contri ) ** 2
                   END SELECT

                   U( iz, ir ) = U( iz, ir ) + CMPLX( Hermite( n, RadiusMax, 2 * RadiusMax ) * contri )
                END IF
             END DO RcvrRanges
             rA  = rB
             nA  = nB
             ir1 = ir2
          END DO Stepping
       END DO Images
    END DO RcvrDepths

  END SUBROUTINE InfluenceCervenyRayCen

  ! **********************************************************************!

  SUBROUTINE InfluenceCervenyCart( U, eps, alpha, iBeamWindow2, RadiusMax )

    ! Paraxial (Cerveny-style) beams in Cartesian coordinates

    INTEGER,          INTENT( IN    ) :: IBeamWindow2
    REAL    (KIND=8), INTENT( IN    ) :: alpha, RadiusMax             ! take-off angle
    COMPLEX,          INTENT( INOUT ) :: U( NRz_per_range, Pos%NRr )  ! complex pressure field
    COMPLEX (KIND=8), INTENT( IN    ) :: eps
    INTEGER          :: KMAHV( MaxN ), KMAH, irA, irB, Image
    REAL    (KIND=8) :: x( 2 ), rayt( 2 ), rayn( 2 ), Tr, Tz, zr, Polarity = 1, &
         c, cimag, cs, cn, csq, gradc( 2 ), crr, crz, czz, rho, deltaz
    COMPLEX (KIND=8) :: pVB( MaxN ), qVB( MaxN ), q, epsV( MaxN ), contri, gammaV( MaxN ), gamma, const
    COMPLEX (KIND=8) :: tau

    ! need to add logic related to NRz_per_range

    ! During reflection imag( q ) is constant and adjacent normals cannot bracket a segment of the TL
    ! line, so no special treatment is necessary

    IF ( Beam%Type( 2 : 2 ) == 'C' ) THEN
       epsV( 1 : Beam%Nsteps ) = i * ABS( ray2D( 1 : Beam%Nsteps )%q( 1 ) / ray2D( 1 : Beam%Nsteps )%q( 2 ) )
    ELSE
       epsV( 1 : Beam%Nsteps ) = eps
    END IF

    pVB( 1 : Beam%Nsteps ) = ray2D( 1 : Beam%Nsteps )%p( 1 ) + epsV( 1 : Beam%Nsteps ) * ray2D( 1 : Beam%Nsteps )%p( 2 )
    qVB( 1 : Beam%Nsteps ) = ray2D( 1 : Beam%Nsteps )%q( 1 ) + epsV( 1 : Beam%Nsteps ) * ray2D( 1 : Beam%Nsteps )%q( 2 )

    IF ( Beam%RunType( 4 : 4 ) == 'R' ) Ratio1 = SQRT( ABS( COS( alpha ) ) )  ! point source

    ! Form gamma and KMAH index
    ! Treatment of KMAH index is incorrect for 'Cerveny' style beam width BeamType

    Stepping0: DO iS = 1, Beam%Nsteps

       rayt = ray2D( iS )%c * ray2D( iS )%t   ! unit tangent
       rayn = [ rayt( 2 ), -rayt( 1 ) ]       ! unit normal

       CALL EvaluateSSP( ray2D( iS )%x, c, cimag, gradc, crr, crz, czz, rho, Freq, 'TAB' )

       csq = c * c
       cS  = DOT_PRODUCT( gradc, rayt )
       cN  = DOT_PRODUCT( gradc, rayn )

       Tr  = rayt(  1 )
       Tz  = rayt(  2 )

       gammaV( iS ) = 0.0
       IF ( qVB( iS ) /= 0.0 ) gammaV( iS ) = 0.5 * ( pVB( iS ) / qVB( iS ) * Tr**2 + 2.0 * cN / csq * Tz * Tr - cS / csq * Tz**2 )

       IF ( iS == 1 ) THEN
          KMAHV( 1 ) = 1
       ELSE
          KMAHV( iS ) = KMAHV( iS - 1 )
          CALL BranchCut( qVB( iS - 1 ), qVB( iS ), Beam%Type, KMAHV( iS ) )
       END IF
    END DO Stepping0

    Stepping: DO iS = 3, Beam%Nsteps
       IF ( ray2D( iS     )%x( 1 ) > Pos%Rr( Pos%NRr ) ) RETURN
       rA = ray2D( iS - 1 )%x( 1 )
       rB = ray2D( iS     )%x( 1 )
       IF ( ABS( rB - rA ) < 1.0D3 * SPACING( rB ) ) CYCLE Stepping   ! don't process duplicate points

       ! Compute upper index on rcvr line
       !!! Assumes r is a vector of equally spaced points
       irA = MAX( MIN( INT( ( rA - Pos%Rr( 1 ) ) / Pos%Delta_r ) + 1, Pos%NRr ), 1 ) ! should be ", 0 )" ?
       irB = MAX( MIN( INT( ( rB - Pos%Rr( 1 ) ) / Pos%Delta_r ) + 1, Pos%NRr ), 1 )

       IF ( irA >= irB ) CYCLE Stepping

       RcvrRanges: DO ir = irA + 1, irB

          W     = ( Pos%Rr( ir ) - rA ) / ( rB - rA )

          x     = ray2D(  iS - 1 )%x    + W * ( ray2D(  iS )%x   -  ray2D(  iS - 1 )%x   )
          rayt  = ray2D(  iS - 1 )%t    + W * ( ray2D(  iS )%t   -  ray2D(  iS - 1 )%t   )
          c     = ray2D(  iS - 1 )%c    + W * ( ray2D(  iS )%c   -  ray2D(  iS - 1 )%c   )
          q     = qVB(    iS - 1 )      + W * ( qVB(    iS )     -  qVB(    iS - 1 )     )
          tau   = ray2D(  iS - 1 )%tau  + W * ( ray2D(  iS )%tau -  ray2D(  iS - 1 )%tau )
          gamma = gammaV( iS - 1 )      + W * ( gammaV( iS )     -  gammaV( iS - 1 )     )

          IF ( AIMAG( gamma ) > 0 ) THEN
             WRITE( PRTFile, * ) 'Unbounded beam'
             WRITE( PRTFile, * ) gammaV( iS - 1 ), gammaV( iS ), gamma
             CYCLE RcvrRanges
          END IF

          const = Ratio1 * SQRT( c * ABS( epsV( iS - 1 ) ) / q )

          ! Get correct branch of SQRT
          KMAH = KMAHV( iS - 1 )
          CALL BranchCut( qVB( iS - 1 ), q, Beam%Type, KMAH )
          IF ( KMAH < 0 ) const = -const

          RcvrDepths: DO iz = 1, NRz_per_range
             zR = Pos%Rz( iz )

             contri = 0.0
             ImageLoop: DO Image = 1, Beam%Nimage
                SELECT CASE ( Image )
                CASE ( 1 )   ! True beam
                   deltaz = zR - x( 2 )
                   Polarity = +1.0D0
                CASE ( 2 )   ! Surface reflected beam
                   deltaz = -zR + 2.0 * Bdry%Top%HS%Depth - x( 2 )
                   Polarity = -1.0D0
                CASE ( 3 )   ! Bottom  reflected beam
                   deltaz = -zR + 2.0 * Bdry%Bot%HS%Depth - x( 2 )
                   Polarity = +1.0D0   ! assumes rigid bottom
                END SELECT

                IF ( omega * AIMAG( gamma ) * deltaz ** 2 < iBeamWindow2 ) &
                     contri =  contri + Polarity * ray2D( iS )%Amp * Hermite( deltaz, RadiusMax, 2.0 * RadiusMax ) * &
                     EXP( -i * ( omega * ( tau + rayt( 2 ) * deltaz + gamma * deltaz**2 ) - ray2D( iS )%Phase ) )
             END DO ImageLoop

             ! contribution to field
             SELECT CASE( Beam%RunType( 1 : 1 ) )
             CASE ( 'C' )        ! coherent
                contri = const * contri
             CASE ( 'I', 'S' )   ! incoherent or semi-coherent
                contri = ABS( const * contri ) ** 2
             END SELECT
             U( iz, ir ) = U( iz, ir ) + CMPLX( contri )
          END DO RcvrDepths
       END DO RcvrRanges
    END DO Stepping

  END SUBROUTINE InfluenceCervenyCart

  ! **********************************************************************!

  SUBROUTINE InfluenceGeoHatRayCen( U, alpha, dalpha )

    ! Geometrically-spreading beams with a hat-shaped beam in ray-centered coordinates

    REAL (KIND=8), INTENT( IN    ) :: alpha, dalpha                 ! take-off angle
    COMPLEX,       INTENT( INOUT ) :: U( NRz_per_range, Pos%NRr )   ! complex pressure field
    INTEGER          :: irA, irB, II
    REAL    (KIND=8) :: nA, nB, zr, L, dq( Beam%Nsteps - 1 )
    REAL    (KIND=8) :: znV( Beam%Nsteps ), rnV( Beam%Nsteps ),  RcvrDeclAngleV ( Beam%Nsteps )
    COMPLEX (KIND=8) :: dtau( Beam%Nsteps - 1 )

    !!! need to add logic related to NRz_per_range

    q0           = ray2D( 1 )%c / Dalpha   ! Reference for J = q0 / q
    SrcDeclAngle = RadDeg * alpha          ! take-off angle in degrees

    dq   = ray2D( 2 : Beam%Nsteps )%q( 1 ) - ray2D( 1 : Beam%Nsteps - 1 )%q( 1 )
    dtau = ray2D( 2 : Beam%Nsteps )%tau    - ray2D( 1 : Beam%Nsteps - 1 )%tau

    ! pre-calculate ray normal based on tangent with c(s) scaling
    znV = -ray2D( 1 : Beam%Nsteps )%t( 1 ) * ray2D( 1 : Beam%Nsteps )%c
    rnV =  ray2D( 1 : Beam%Nsteps )%t( 2 ) * ray2D( 1 : Beam%Nsteps )%c

    RcvrDeclAngleV( 1 : Beam%Nsteps ) = RadDeg * ATAN2( ray2D( 1 : Beam%Nsteps )%t( 2 ), ray2D( 1 : Beam%Nsteps )%t( 1 ) )

    ! During reflection imag(q) is constant and adjacent normals cannot bracket a segment of the TL
    ! line, so no special treatment is necessary

    IF ( Beam%RunType( 4 : 4 ) == 'R' ) Ratio1 = SQRT( ABS( COS( alpha ) ) )  ! point source

    ray2D( 1 : Beam%Nsteps )%Amp = Ratio1 * SQRT( ray2D( 1 : Beam%Nsteps )%c ) * ray2D( 1 : Beam%Nsteps )%Amp   ! pre-apply some scaling

    RcvrDepths: DO iz = 1, NRz_per_range
       zR = Pos%Rz( iz )

       phase = 0.0
       qOld  = ray2D( 1 )%q( 1 )       ! used to track KMAH index

       IF ( ABS( znV( 1 ) ) < 1D-6 ) THEN   ! normal parallel to horizontal receiver line
          nA  = 1D10
          rA  = 1D10
          irA = 1
       ELSE
          nA  = ( zR - ray2D( 1 )%x( 2 )   ) / znV( 1 )
          rA  = ray2D( 1 )%x( 1 ) + nA * rnV( 1 )
          !!! following assumes uniform spacing in Pos%r
          irA = MAX( MIN( INT( ( rA - Pos%Rr( 1 ) ) / Pos%Delta_r ) + 1, Pos%NRr ), 1 ) ! index of receiver
       END IF

       Stepping: DO iS = 2, Beam%Nsteps

          ! Compute ray-centered coordinates, (znV, rnV)

          IF ( ABS( znV( iS ) ) < 1D-10 ) CYCLE Stepping   ! If normal parallel to TL-line, skip to next step on ray
          nB  = ( zR - ray2D( iS )%x( 2 )   ) / znV( iS )
          rB  = ray2D( iS )%x( 1 ) + nB * rnV( iS )

          !!! following assumes uniform spacing in Pos%r
          irB = MAX( MIN( INT( ( rB - Pos%Rr( 1 ) ) / Pos%Delta_r ) + 1, Pos%NRr ), 1 ) ! index of receiver

          ! detect and skip duplicate points (happens at boundary reflection)
          IF ( ABS( ray2D( iS )%x( 1 ) - ray2D( iS - 1 )%x( 1 ) ) < 1.0D3 * SPACING( ray2D( iS )%x( 1 ) ) .OR. irA == irB ) THEN
             rA  = rB
             nA  = nB
             irA = irB
             CYCLE Stepping
          END IF

          !!! this should be pre-computed
          q  = ray2D( iS - 1 )%q( 1 )
          IF ( q <= 0.0d0 .AND. qOld > 0.0d0 .OR. q >= 0.0d0 .AND. qOld < 0.0d0 ) phase = phase + pi / 2.  ! phase shifts at caustics
          qold = q

          RcvrDeclAngle = RcvrDeclAngleV( iS )

          ! *** Compute contributions to bracketed receivers ***

          II = 0
          IF ( irB <= irA ) II = 1   ! going backwards in range

          RcvrRanges: DO ir = irA + 1 - II, irB + II, SIGN( 1, irB - irA )  ! Compute influence for each rcvr
             W = ( Pos%Rr( ir ) - rA ) / ( rB - rA )
             n = ABS( nA                + W * ( nB - nA ) )
             q = ray2D( iS - 1 )%q( 1 ) + W * dq( iS - 1 )     ! interpolated amplitude
             L = ABS( q ) / q0   ! beam radius

             IF ( n < L ) THEN   ! in beamwindow?
                delay    = ray2D( iS - 1 )%tau + W * dtau( iS - 1 )   ! interpolated delay
                const    = ray2D( iS )%Amp / SQRT( ABS( q ) ) 
                W        = ( L - n ) / L   ! hat function: 1 on center, 0 on edge
                Amp      = const * W
                phaseInt = ray2D( iS - 1 )%Phase + phase
                !!! this should be precomputed
                IF ( q <= 0.0d0 .AND. qOld > 0.0d0 .OR. q >= 0.0d0 .AND. qOld < 0.0d0 ) phaseInt = phaseInt + pi / 2.   ! phase shifts at caustics

                CALL ApplyContribution( U( iz, ir ) )
             END IF
          END DO RcvrRanges
          rA  = rB
          nA  = nB
          irA = irB
       END DO Stepping
    END DO RcvrDepths

  END SUBROUTINE InfluenceGeoHatRayCen

  ! **********************************************************************!

  SUBROUTINE InfluenceGeoHatCart( U, alpha, Dalpha )

    ! Geometric, hat-shaped beams in Cartesisan coordinates

    REAL (KIND=8), INTENT( IN    ) :: alpha, dalpha                 ! take-off angle, angular spacing
    COMPLEX,       INTENT( INOUT ) :: U( NRz_per_range, Pos%NRr )   ! complex pressure field
    INTEGER          :: irT( 1 ), irTT
    REAL    (KIND=8) :: x_ray( 2 ), rayt( 2 ), rayn( 2 ), x_rcvr( 2, NRz_per_range ), rLen, RadiusMax, zMin, zMax, dqds
    COMPLEX (KIND=8) :: dtauds

    q0           = ray2D( 1 )%c / Dalpha   ! Reference for J = q0 / q
    SrcDeclAngle = RadDeg * alpha          ! take-off angle in degrees
    phase        = 0.0
    qOld         = ray2D( 1 )%q( 1 )       ! used to track KMAH index
    rA           = ray2D( 1 )%x( 1 )       ! range at start of ray

    ! what if never satisfied?
    ! what if there is a single receiver (ir = 0 possible)
    irT = MINLOC( Pos%Rr( 1 : Pos%NRr ), MASK = Pos%Rr( 1 : Pos%NRr ) > rA )   ! find index of first receiver to the right of rA
    ir  = irT( 1 )
    IF ( ray2D( 1 )%t( 1 ) < 0.0d0 .AND. ir > 1 ) ir = ir - 1  ! if ray is left-traveling, get the first receiver to the left of rA

    IF ( Beam%RunType( 4 : 4 ) == 'R' ) Ratio1 = SQRT( ABS( COS( alpha ) ) )  ! point source

    Stepping: DO iS = 2, Beam%Nsteps
       rB     = ray2D( iS     )%x( 1 )
       x_ray  = ray2D( iS - 1 )%x

       ! compute normalized tangent (compute it because we need to measure the step length)
       rayt = ray2D( iS )%x - ray2D( iS - 1 )%x
       rlen = NORM2( rayt )
       IF ( rlen < 1.0D3 * SPACING( ray2D( iS )%x( 1 ) ) ) CYCLE Stepping  ! if duplicate point in ray, skip to next step along the ray
       rayt = rayt / rlen                    ! unit tangent to ray
       rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal  to ray
       RcvrDeclAngle = RadDeg * ATAN2( rayt( 2 ), rayt( 1 ) )

       dqds   = ray2D( iS )%q( 1 ) - ray2D( iS - 1 )%q( 1 )
       dtauds = ray2D( iS )%tau    - ray2D( iS - 1 )%tau

       q  = ray2D( iS - 1 )%q( 1 )
       IF ( q <= 0.0d0 .AND. qOld > 0.0d0 .OR. q >= 0.0d0 .AND. qOld < 0.0d0 ) phase = phase + pi / 2.   ! phase shifts at caustics
       qold = q

       RadiusMax = MAX( ABS( ray2D( iS - 1 )%q( 1 ) ), ABS( ray2D( iS )%q( 1 ) ) ) / q0 / ABS( rayt( 1 ) ) ! beam radius projected onto vertical line

       ! depth limits of beam
       IF ( ABS( rayt( 1 ) ) > 0.5 ) THEN   ! shallow angle ray
          zmin   = min( ray2D( iS - 1 )%x( 2 ), ray2D( iS )%x( 2 ) ) - RadiusMax
          zmax   = max( ray2D( iS - 1 )%x( 2 ), ray2D( iS )%x( 2 ) ) + RadiusMax
       ELSE                                 ! steep angle ray
          zmin = -HUGE( zmin )
          zmax = +HUGE( zmax )
       END IF

       ! compute beam influence for this segment of the ray
       RcvrRanges: DO
          ! is Rr( ir ) contained in [ rA, rB )? Then compute beam influence
          IF ( Pos%Rr( ir ) >= MIN( rA, rB ) .AND. Pos%Rr( ir ) < MAX( rA, rB ) ) THEN
             
             x_rcvr( 1, 1 : NRz_per_range ) = Pos%Rr( ir )
             IF ( Beam%RunType( 5 : 5 ) == 'I' ) THEN
                x_rcvr( 2, 1 ) = Pos%Rz( ir )                  ! irregular   grid
             ELSE
                x_rcvr( 2, 1 : NRz_per_range ) = Pos%Rz( 1 : NRz_per_range )   ! rectilinear grid
             END IF

             RcvrDepths: DO iz = 1, NRz_per_range
                IF ( x_rcvr( 2, iz ) < zmin .OR. x_rcvr( 2, iz ) > zmax ) CYCLE RcvrDepths

                s         =      DOT_PRODUCT( x_rcvr( :, iz ) - x_ray, rayt ) / rlen ! proportional distance along ray
                n         = ABS( DOT_PRODUCT( x_rcvr( :, iz ) - x_ray, rayn ) )      ! normal distance to ray
                q         = ray2D( iS - 1 )%q( 1 ) + s * dqds               ! interpolated amplitude
                RadiusMax = ABS( q / q0 )                                   ! beam radius

                IF ( n < RadiusMax ) THEN
                   delay    = ray2D( iS - 1 )%tau + s * dtauds              ! interpolated delay
                   const    = Ratio1 * SQRT( ray2D( iS )%c / ABS( q ) ) * ray2D( iS )%Amp
                   W        = ( RadiusMax - n ) / RadiusMax   ! hat function: 1 on center, 0 on edge
                   Amp      = const * W
                   phaseInt = ray2D( iS - 1 )%Phase + phase
                   IF ( q <= 0.0d0 .AND. qOld > 0.0d0 .OR. q >= 0.0d0 .AND. qOld < 0.0d0 ) phaseInt = phaseInt + pi / 2.   ! phase shifts at caustics

                   CALL ApplyContribution( U( iz, ir ) )
                END IF
             END DO RcvrDepths
          END IF

          ! bump receiver index, ir, towards rB
          IF ( Pos%Rr( ir ) < rB ) THEN
             IF ( ir >= Pos%NRr        ) EXIT  ! go to next step on ray
             irTT = ir + 1                     ! bump right
             IF ( Pos%Rr( irTT ) >= rB ) EXIT
          ELSE
             IF ( ir <= 1              ) EXIT  ! go to next step on ray
             irTT = ir - 1                     ! bump left
             IF ( Pos%Rr( irTT ) <= rB ) EXIT
          END IF
          ir = irTT
       END DO RcvrRanges

       rA = rB
    END DO Stepping

  END SUBROUTINE InfluenceGeoHatCart

  ! **********************************************************************!

  SUBROUTINE InfluenceGeoGaussianCart( U, alpha, Dalpha )

    ! Geometric, Gaussian beams in Cartesian coordintes

    INTEGER,       PARAMETER       :: BeamWindow = 4               ! beam window: kills beams outside e**(-0.5 * ibwin**2 )
    REAL (KIND=8), INTENT( IN    ) :: alpha, dalpha                ! take-off angle, angular spacing
    COMPLEX,       INTENT( INOUT ) :: U( NRz_per_range, Pos%NRr )  ! complex pressure field
    INTEGER          :: irT( 1 ), irTT
    REAL    (KIND=8) :: x_ray( 2 ), rayt( 2 ), rayn( 2 ), x_rcvr( 2 ), rLen, RadiusMax, zMin, zMax, sigma, lambda, A, dqds
    COMPLEX (KIND=8) :: dtauds

    q0           = ray2D( 1 )%c / Dalpha   ! Reference for J = q0 / q
    SrcDeclAngle = RadDeg * alpha          ! take-off angle in degrees
    phase        = 0
    qOld         = ray2D( 1 )%q( 1 )       ! used to track KMAH index
    rA           = ray2D( 1 )%x( 1 )       ! range at start of ray

    ! what if never satisfied?
    ! what if there is a single receiver (ir = 0 possible)

    irT = MINLOC( Pos%Rr( 1 : Pos%NRr ), MASK = Pos%Rr( 1 : Pos%NRr ) > rA )      ! find index of first receiver to the right of rA
    ir  = irT( 1 )

    IF ( ray2D( 1 )%t( 1 ) < 0.0d0 .AND. ir > 1 ) ir = ir - 1  ! if ray is left-traveling, get the first receiver to the left of rA

    ! sqrt( 2 * pi ) represents a sum of Gaussians in free space
    IF ( Beam%RunType( 4 : 4 ) == 'R' ) THEN
       Ratio1 = SQRT( ABS( COS( alpha ) ) ) / SQRT( 2. * pi )   ! point source
    ELSE
       Ratio1 = 1 / SQRT( 2. * pi )                             ! line  source
    END IF

    Stepping: DO iS = 2, Beam%Nsteps

       rB    = ray2D( iS     )%x( 1 )
       x_ray = ray2D( iS - 1 )%x

       ! compute normalized tangent (compute it because we need to measure the step length)
       rayt = ray2D( iS )%x - ray2D( iS - 1 )%x
       rlen = NORM2( rayt )
       IF ( rlen < 1.0D3 * SPACING( ray2D( iS )%x( 1 ) ) ) CYCLE Stepping  ! if duplicate point in ray, skip to next step along the ray
       rayt = rayt / rlen
       rayn = [ -rayt( 2 ), rayt( 1 ) ]      ! unit normal to ray
       RcvrDeclAngle = RadDeg * ATAN2( rayt( 2 ), rayt( 1 ) )

       dqds   = ray2D( iS )%q( 1 ) - ray2D( iS - 1 )%q( 1 )
       dtauds = ray2D( iS )%tau    - ray2D( iS - 1 )%tau

       q  = ray2D( iS - 1 )%q( 1 )
       IF ( q <= 0.0 .AND. qOld > 0.0 .OR. q >= 0.0 .AND. qOld < 0.0 ) phase = phase + pi / 2.   ! phase shifts at caustics
       qold = q

       ! calculate beam width
       lambda    = ray2D( iS - 1 )%c / freq
       sigma     = MAX( ABS( ray2D( iS - 1 )%q( 1 ) ), ABS( ray2D( iS )%q( 1 ) ) ) / q0 / ABS( rayt( 1 ) ) ! beam radius projected onto vertical line
       sigma     = MAX( sigma, MIN( 0.2 * freq * REAL( ray2D( iS )%tau ), pi * lambda ) )
       RadiusMax = BeamWindow * sigma

       ! depth limits of beam
       IF ( ABS( rayt( 1 ) ) > 0.5 ) THEN   ! shallow angle ray
          zmin   = min( ray2D( iS - 1 )%x( 2 ), ray2D( iS )%x( 2 ) ) - RadiusMax
          zmax   = max( ray2D( iS - 1 )%x( 2 ), ray2D( iS )%x( 2 ) ) + RadiusMax
       ELSE                                 ! steep angle ray
          zmin = -HUGE( zmin )
          zmax = +HUGE( zmax )
       END IF

       ! compute beam influence for this segment of the ray
       RcvrRanges: DO
          ! is Rr( ir ) contained in [ rA, rB )? Then compute beam influence
          IF ( Pos%Rr( ir ) >= MIN( rA, rB ) .AND. Pos%Rr( ir ) < MAX( rA, rB ) ) THEN

             RcvrDepths: DO iz = 1, NRz_per_range
                IF ( Beam%RunType( 5 : 5 ) == 'I' ) THEN
                   x_rcvr = [ Pos%Rr( ir ), DBLE( Pos%Rz( ir ) ) ]   ! irregular   grid
                ELSE
                   x_rcvr = [ Pos%Rr( ir ), DBLE( Pos%Rz( iz ) ) ]   ! rectilinear grid
                END IF
                IF ( x_rcvr( 2 ) < zmin .OR. x_rcvr( 2 ) > zmax ) CYCLE RcvrDepths

                s      =      DOT_PRODUCT( x_rcvr - x_ray, rayt ) / rlen  ! proportional distance along ray
                n      = ABS( DOT_PRODUCT( x_rcvr - x_ray, rayn ) )       ! normal distance to ray
                q      = ray2D( iS - 1 )%q( 1 ) + s * dqds                ! interpolated amplitude
                sigma  = ABS( q / q0 )                                    ! beam radius
                sigma  = MAX( sigma, MIN( 0.2 * freq * REAL( ray2D( iS )%tau ), pi * lambda ) )  ! min pi * lambda, unless near

                IF ( n < BeamWindow * sigma ) THEN   ! Within beam window?
                   A        = ABS( q0 / q )
                   delay    = ray2D( iS - 1 )%tau + s * dtauds     ! interpolated delay
                   const    = Ratio1 * SQRT( ray2D( iS )%c / ABS( q ) ) * ray2D( iS )%Amp
                   W        = EXP( -0.5 * ( n / sigma ) ** 2 ) / ( sigma * A )   ! Gaussian decay
                   Amp      = const * W
                   phaseInt = ray2D( iS - 1 )%Phase + phase
                   IF ( q <= 0.0d0 .AND. qOld > 0.0d0 .OR. q >= 0.0d0 .AND. qOld < 0.0d0 ) phaseInt = phaseInt + pi / 2.  ! phase shifts at caustics

                   CALL ApplyContribution( U( iz, ir ) )
                END IF
             END DO RcvrDepths
          END IF

          ! receiver not bracketed; bump receiver index, ir, towards rB
          IF ( rB > Pos%Rr( ir ) ) THEN
             IF ( ir >= Pos%NRr        ) EXIT   ! go to next step on ray
             irTT = ir + 1                      ! bump right
             IF ( Pos%Rr( irTT ) >= rB ) EXIT   ! go to next step on ray
          ELSE
             IF ( ir <= 1              ) EXIT   ! go to next step on ray
             irTT = ir - 1                      ! bump left
             IF ( Pos%Rr( irTT ) <= rB ) EXIT   ! go to next step on ray
          END IF
          ir = irTT

       END DO RcvrRanges

       rA = rB
    END DO Stepping

  END SUBROUTINE InfluenceGeoGaussianCart

  ! **********************************************************************!
  
  SUBROUTINE ApplyContribution( U )
    COMPLEX, INTENT( INOUT ) :: U

    SELECT CASE( Beam%RunType( 1 : 1 ) )
    CASE ( 'E' )                ! eigenrays
       IF ( Title( 1 :  9 ) == 'BELLHOP- ' ) THEN   ! BELLHOP run
          CALL WriteRay2D( SrcDeclAngle, iS )
       ELSE                                         ! BELLHOP3D run
          CALL WriteRay3D( DegRad * SrcDeclAngle, DegRad * SrcAzimAngle, is )   ! produces no output if NR=1
       END IF
    CASE ( 'A', 'a' )           ! arrivals
       CALL AddArr( omega, iz, ir, Amp, phaseInt, delay, SrcDeclAngle, RcvrDeclAngle, ray2D( iS )%NumTopBnc, ray2D( iS )%NumBotBnc )
    CASE ( 'C' )                ! coherent TL
       U = U + CMPLX( Amp * EXP( -i * ( omega * delay - phaseInt ) ) )
                     ! omega * n * n / ( 2 * ray2d( iS )%c**2 * delay ) ) ) )   ! curvature correction
    CASE DEFAULT                ! incoherent/semicoherent TL
       IF ( Beam%Type( 1 : 1 ) == 'B' ) THEN   ! Gaussian beam
          U = U + SNGL( SQRT( 2. * pi ) * ( const * EXP( AIMAG( omega * delay ) ) ) ** 2 * W )
       ELSE
          U = U + SNGL(                   ( const * EXP( AIMAG( omega * delay ) ) ) ** 2 * W )
       END IF
    END SELECT

  END SUBROUTINE ApplyContribution
                 
  ! **********************************************************************!

  SUBROUTINE InfluenceSGB( U, alpha, Dalpha )

    ! Bucker's Simple Gaussian Beams in Cartesian coordinates

    REAL (KIND=8), INTENT( IN    ) :: alpha, dalpha                 ! take-off angle, angular spacing
    COMPLEX,       INTENT( INOUT ) :: U( NRz_per_range, Pos%NRr )   ! complex pressure field
    REAL    (KIND=8)  :: x( 2 ), rayt( 2 ), A, beta, cn, CPA, deltaz, DS, sint, SX1, thet
    COMPLEX (KIND=8)  :: contri, tau

    Ratio1 = SQRT(  COS( alpha ) )
    phase  = 0
    qOld   = 1.0
    BETA   = 0.98  ! Beam Factor
    A      = -4.0 * LOG( BETA ) / Dalpha**2
    CN     = Dalpha * SQRT( A / pi )
    rA     = ray2D( 1 )%x( 1 )
    ir     = 1

    Stepping: DO iS = 2, Beam%Nsteps

       rB = ray2D( iS )%x( 1 )

       ! phase shifts at caustics
       q  = ray2D( iS - 1 )%q( 1 )
       IF ( q < 0.0d0 .AND. qOld >= 0.0d0 .OR. q > 0.0d0 .AND. qOld <= 0.0d0 ) phase = phase + pi / 2.
       qold = q

       RcvrRanges: DO WHILE ( ABS( rB - rA ) > 1.0D3 * SPACING( rA ) .AND. rB > Pos%Rr( ir ) )   ! Loop over bracketed receiver ranges

          W     = ( Pos%Rr( ir ) - rA ) / ( rB - rA )
          x     = ray2D( iS - 1 )%x      + W * ( ray2D( iS )%x      - ray2D( iS - 1 )%x )
          rayt  = ray2D( iS - 1 )%t      + W * ( ray2D( iS )%t      - ray2D( iS - 1 )%t )
          q     = ray2D( iS - 1 )%q( 1 ) + W * ( ray2D( iS )%q( 1 ) - ray2D( iS - 1 )%q( 1 ) )
          tau   = ray2D( iS - 1 )%tau    + W * ( ray2D( iS )%tau    - ray2D( iS - 1 )%tau )

          ! following is incorrect because ray doesn't always use a step of deltas
          SINT  = ( iS - 1 ) * Beam%deltas + W * Beam%deltas

          IF ( q < 0.0d0 .AND. qOld >= 0.0d0 .OR. q > 0.0d0 .AND. qOld <= 0.0d0 ) phase = phase + pi / 2. ! phase shifts at caustics

          RcvrDepths: DO iz = 1, NRz_per_range
             deltaz =  Pos%Rz( iz ) - x( 2 )   ! ray to rcvr distance
             ! Adeltaz    = ABS( deltaz )
             ! IF ( Adeltaz < RadiusMax ) THEN
             SELECT CASE( Beam%RunType( 1 : 1 ) )
             CASE ( 'E' )         ! eigenrays
                SrcDeclAngle = RadDeg * alpha   ! take-off angle in degrees
                CALL WriteRay2D( SrcDeclAngle, iS )
             CASE DEFAULT         ! coherent TL
                CPA    = ABS( deltaz * ( rB - rA ) ) / SQRT( ( rB - rA )**2 + ( ray2D( iS )%x( 2 ) - ray2D( iS - 1 )%x( 2 ) )**2  )
                DS     = SQRT( deltaz **2 - CPA **2 )
                SX1    = SINT + DS
                thet   = ATAN( CPA / SX1 )
                delay  = tau + rayt( 2 ) * deltaz
                contri = Ratio1 * CN * ray2D( iS )%Amp * EXP( -A * thet ** 2 - &
                     i * ( omega * delay - ray2D( iS )%Phase - phase ) ) / SQRT( SX1 )
                U( iz, ir ) = U( iz, ir ) + CMPLX( contri )
             END SELECT
             ! END IF
          END DO RcvrDepths

          qOld = q
          ir   = ir + 1
          IF ( ir > Pos%NRr ) RETURN
       END DO RcvrRanges

       rA = rB
    END DO Stepping

  END SUBROUTINE InfluenceSGB

  ! **********************************************************************!

  SUBROUTINE BranchCut( q1C, q2C, BeamType, KMAH )

    ! Checks for a branch cut crossing and updates KMAH accordingly

    COMPLEX  (KIND=8), INTENT( IN )    :: q1C, q2C
    CHARACTER (LEN=4), INTENT( IN )    :: BeamType
    INTEGER,           INTENT( INOUT ) :: KMAH
    REAL     (KIND=8)                  :: q1, q2

    SELECT CASE ( BeamType( 2 : 2 ) )
    CASE ( 'W' )   ! WKBeams
       q1 = REAL( q1C )
       q2 = REAL( q2C )
       IF ( ( q1 < 0.0 .AND. q2 >= 0.0 ) .OR. &
            ( q1 > 0.0 .AND. q2 <= 0.0 ) ) KMAH = -KMAH
    CASE DEFAULT
       IF ( REAL( q2C ) < 0.0 ) THEN
          q1 = AIMAG( q1C )
          q2 = AIMAG( q2C )
          IF ( ( q1 < 0.0 .AND. q2 >= 0.0 ) .OR. &
               ( q1 > 0.0 .AND. q2 <= 0.0 ) ) KMAH = -KMAH
       END IF
    END SELECT

  END SUBROUTINE BranchCut

  ! **********************************************************************!

  SUBROUTINE ScalePressure( Dalpha, c, r, U, NRz, Nr, RunType, freq )

    ! Scale the pressure field

    REAL,              PARAMETER       :: pi = 3.14159265
    INTEGER,           INTENT( IN    ) :: NRz, Nr
    REAL     (KIND=8), INTENT( IN    ) :: r( Nr )         ! ranges
    REAL     (KIND=8), INTENT( IN    ) :: Dalpha, freq, c ! angular spacing between rays, source frequency, nominal sound speed
    COMPLEX,           INTENT( INOUT ) :: U( NRz, Nr )    ! Pressure field
    CHARACTER (LEN=5), INTENT( IN    ) :: RunType
    REAL     (KIND=8)                  :: const, factor

    ! Compute scale factor for field
    SELECT CASE ( RunType( 2 : 2 ) )
    CASE ( 'C' )   ! Cerveny Gaussian beams in Cartesian coordinates
       const = -Dalpha * SQRT( freq ) / c
    CASE ( 'R' )   ! Cerveny Gaussian beams in Ray-centered coordinates
       const = -Dalpha * SQRT( freq ) / c
    CASE DEFAULT
       const = -1.0
    END SELECT

    IF ( RunType( 1 : 1 ) /= 'C' ) U = SQRT( REAL( U ) ) ! For incoherent run, convert intensity to pressure

    ! scale and/or incorporate cylindrical spreading
    Ranges: DO ir = 1, Nr
       IF ( RunType( 4 : 4 ) == 'X' ) THEN   ! line source
          factor = -4.0 * SQRT( pi ) * const
       ELSE                                  ! point source
          IF ( r ( ir ) == 0 ) THEN
             factor = 0.0D0                  ! avoid /0 at origin, return pressure = 0
          ELSE
             factor = const / SQRT( ABS( r( ir ) ) )
          END IF
       END IF
       U( :, ir ) = SNGL( factor ) * U( :, ir )
    END DO Ranges

  END SUBROUTINE ScalePressure

  ! **********************************************************************!

  REAL (KIND=8 ) FUNCTION Hermite( x, x1, x2 )

    ! Calculates a smoothing function based on the h0 hermite cubic
    ! x is the point where the function is to be evaluated
    ! returns:
    ! [  0, x1  ] = 1
    ! [ x1, x2  ] = cubic taper from 1 to 0
    ! [ x2, inf ] = 0

    REAL (KIND=8 ), INTENT( IN  ) :: x, x1, x2
    REAL (KIND=8 )                :: Ax, u

    Ax  = ABS( x  )

    IF ( Ax <= x1 ) THEN
       Hermite = 1.0d0
    ELSE IF ( Ax >= x2 ) THEN
       Hermite = 0.0d0
    ELSE
       u       = ( Ax - x1 ) / ( x2 - x1 )
       Hermite = ( 1.0d0 + 2.0d0 * u ) * ( 1.0d0 - u ) ** 2
    END IF

    !hermit = hermit / ( 0.5 * ( x1 + x2 ) )

  END FUNCTION Hermite

END MODULE Influence
