MODULE Reflect3DMod

  USE bellhopMod

  IMPLICIT NONE
CONTAINS

  SUBROUTINE Reflect3D( is, HS, BotTop, nBdry, z_xx, z_xy, z_yy, kappa_xx, kappa_xy, kappa_yy, RefC, Npts )

    ! *** If reflecting off surface need jump conditions ***

    USE RefCoef
    USE sspMod

    INTEGER,           INTENT( IN    ) :: Npts                        ! Number of points in the reflection coefficient
    REAL    (KIND= 8), INTENT( INOUT ) :: nBdry( 3 )                  ! normal to the boundary (changes if cone reflection)
    CHARACTER (LEN=3), INTENT( IN    ) :: BotTop                      ! bottom or top flag
    TYPE( HSInfo ),    INTENT( IN    ) :: HS                          ! halfspace parameters
    TYPE(ReflectionCoef), INTENT( IN ) :: RefC( NPts )                ! reflection coefficient
    INTEGER,           INTENT( INOUT ) :: is                          ! index of the ray step
    INTEGER          :: is1
    REAL    (KIND=8) :: rayt( 3 ), rayn1( 3 ), rayn2( 3 )             ! unit ray tangent and normals
    REAL    (KIND=8) :: rayt_tilde( 3 ), rayn1_tilde( 3 ), rayn2_tilde( 3 ), cn1jump, cn2jump, csjump 
    REAL    (KIND=8) :: c, cimag, gradc( 3 ), cxx, cyy, czz, cxy, cxz, cyz, rho   ! derivatives of sound speed in cartesian coordinates
    REAL    (KIND=8) :: RM, R1, R2, R3, Tg, Th                        ! curvature corrections on reflection
    COMPLEX (KIND=8) :: gamma1, gamma2, gamma1Sq, gamma2Sq, GK, Refl
    TYPE(ReflectionCoef) :: RInt
    REAL    (KIND=8) :: tBdry( 3 )                                    ! tangent to the boundary
    REAL    (KIND=8) :: e1( 3 ), e2( 3 )                              ! ray normals for ray-centered coordinates
    REAL    (KIND=8) :: p_tilde_in(  2 ), p_hat_in(  2 ), q_tilde_in(  2 ), q_hat_in(  2 ), p_tilde_out( 2 ), p_hat_out( 2 )
    REAL    (KIND=8) :: z_xx, z_xy, z_yy, kappa_xx, kappa_xy, kappa_yy, t_rot( 2 ), n_rot( 2 )
    REAL    (KIND=8) :: RotMat( 2, 2 ), kappaMat( 2, 2 ), DMat( 2, 2 ), DMatTemp( 2, 2 )

    is  = is + 1
    is1 = is + 1

    !CALL ConeFormulas(    z_xx, z_xy, z_yy, nBdry, xs, ray3D( is )%x, BotTop ) ! analytic formulas for the curvature of the seamount
    !CALL ParabotFormulas( z_xx, z_xy, z_yy, nBdry )                    ! analytic formulas for the curvature of the parabolic bottom
    !write( *, * ) 'z_xx, z_xy, z_yy', z_xx, z_xy, z_yy, 'nBdry', nBdry
    kappaMat( 1, 1 ) = z_xx / 2
    kappaMat( 1, 2 ) = z_xy / 2
    kappaMat( 2, 1 ) = z_xy / 2
    kappaMat( 2, 2 ) = z_yy / 2

    ! get normal and tangential components of ray in the reflecting plane

    Th = DOT_PRODUCT( ray3D( is )%t, nBdry )  ! component of ray tangent normal to boundary

    tBdry = ray3D( is )%t - Th * nBdry        ! component of ray tangent along the boundary, in the reflection plane
    tBdry = tBdry / NORM2( tBdry )

    Tg = DOT_PRODUCT( ray3D( is )%t, tBdry )  ! component of ray tangent along the boundary

    ray3D( is1 )%NumTopBnc = ray3D( is )%NumTopBnc
    ray3D( is1 )%NumBotBnc = ray3D( is )%NumBotBnc
    ray3D( is1 )%x         = ray3D( is )%x
    ray3D( is1 )%t         = ray3D( is )%t - 2.0 * Th * nBdry   ! reflect the ray

    CALL EvaluateSSP3D( ray3D( is1 )%x, c, cimag, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, 'TAB' )
    ray3D( is1 )%c         = c
    ray3D( is1 )%tau       = ray3D( is )%tau

    ! Calculate the ray normals, rayn1, rayn2, and a unit tangent

    CALL CalcTangent_Normals( ray3D( is  )%t, nBdry, rayt,       rayn1,       rayn2       ) ! incident
    CALL CalcTangent_Normals( ray3D( is1 )%t, nBdry, rayt_tilde, rayn1_tilde, rayn2_tilde ) ! reflected

    ! rotation matrix to get surface curvature in and perpendicular to the reflection plane
    ! we use only the first two elements of the vectors because we want the projection in the x-y plane
    t_rot = rayt(  1 : 2 ) / NORM2( rayt(  1 : 2 ) )
    n_rot = rayn2( 1 : 2 ) / NORM2( rayn2( 1 : 2 ) )

    RotMat( 1 : 2, 1 ) = t_rot
    RotMat( 1 : 2, 2 ) = n_rot

    ! apply the rotation to get the matrix D of curvatures (see Popov 1977 for definition of DMat)
    ! DMat = RotMat^T * kappaMat * RotMat, with RotMat anti-symmetric
    DMatTemp = MATMUL( 2 * kappaMat, RotMat )
    DMat     = MATMUL( TRANSPOSE( RotMat ), DMatTemp )

    ! normal and tangential derivatives of the sound speed
    cn1jump =  DOT_PRODUCT( gradc, -rayn1_tilde - rayn1 )
    cn2jump =  DOT_PRODUCT( gradc, -rayn2_tilde - rayn2 )
    csjump  = -DOT_PRODUCT( gradc,  rayt_tilde  - rayt  )

!!! not sure if cn2 needs a sign flip also
!!$  IF ( BotTop == 'TOP' ) THEN
!!$     cn1jump = -cn1jump    ! flip sign for top reflection
!!$     cn2jump = -cn2jump    ! flip sign for top reflection
!!$  END IF

    RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
    CALL CurvatureCorrection2( ray3D( is ), ray3D( is1 ) )

    ! account for phase change

    SELECT CASE ( HS%BC )
    CASE ( 'R' )                 ! rigid
       ray3D( is1 )%Amp   = ray3D( is )%Amp
       ray3D( is1 )%Phase = ray3D( is )%Phase
    CASE ( 'V' )                 ! vacuum
       ray3D( is1 )%Amp   = ray3D( is )%Amp
       ray3D( is1 )%Phase = ray3D( is )%Phase + pi
    CASE ( 'F' )                 ! file
       RInt%theta = RadDeg * ABS( ATAN2( Th, Tg ) )           ! angle of incidence (relative to normal to bathymetry)
       IF ( RInt%theta > 90 ) RInt%theta = 180. - RInt%theta  ! reflection coefficient is symmetric about 90 degrees
       CALL InterpolateReflectionCoefficient( RInt, RefC, Npts, PRTFile )
       ray3D( is1 )%Amp   = ray3D( is )%Amp * RInt%R
       ray3D( is1 )%Phase = ray3D( is )%Phase + RInt%phi
    CASE ( 'A', 'G' )            ! half-space
       GK       = omega * Tg     ! wavenumber in direction parallel to bathymetry
       gamma1Sq = ( omega / c     ) ** 2 - GK ** 2 - i * tiny( omega )   ! tiny prevents g95 giving -zero, and wrong branch cut
       gamma2Sq = ( omega / HS%cP ) ** 2 - GK ** 2 - i * tiny( omega )
       gamma1   = SQRT( -gamma1Sq )
       gamma2   = SQRT( -gamma2Sq )

       Refl = ( HS%rho * gamma1 - gamma2 ) / ( HS%rho * gamma1 + gamma2 )
       !write( *, * ) abs( Refl ), c, HS%cp, rho, HS%rho       
       ! Hack to make a wall (where the bottom slope is more than 80 degrees) be a perfect reflector
       !!!!!!!!!!!
       IF ( ABS( RadDeg * ATAN2( nBdry( 3 ), NORM2( nBdry( 1 : 2 ) ) ) ) < 0 ) THEN   ! was 60 degrees
          Refl = 1
       END IF

       IF ( ABS( Refl ) < 1.0E-5 ) THEN   ! kill a ray that has lost its energy in reflection
          ray3D( is1 )%Amp   = 0.0
          ray3D( is1 )%Phase = ray3D( is )%Phase
       ELSE
          ray3D( is1 )%Amp   = ABS( Refl ) * ray3D(  is )%Amp
          ray3D( is1 )%Phase = ray3D( is )%Phase + ATAN2( AIMAG( Refl ), REAL( Refl ) )
       ENDIF
    END SELECT

  CONTAINS
    SUBROUTINE CurvatureCorrection2( ray, rayOut )

      USE RayNormals

      TYPE( ray3DPt )    :: ray, rayOut

      ! Note that Tg, Th need to be multiplied by c to normalize tangent; hence, c^2 below
      ! added the SIGN in R2 to make ati and bty have a symmetric effect on the beam
      ! not clear why that's needed

      R1 = 2 / c ** 2 * DMat( 1, 1 ) / Th + RM * ( 2 * cn1jump - RM * csjump ) / c ** 2
      R2 = 2 / c *      DMat( 1, 2 ) * SIGN( 1.0D0, -Th )      + RM * cn2jump  / c ** 2
      R3 = 2 *          DMat( 2, 2 ) * Th

      ! z-component of unit tangent is sin( theta ); we want cos( theta )
      !R1 = R1 * ( 1 - ( ray%c * ray%t( 3 ) ) ** 2 )
      ! write( *, * )  1 - ( ray%c * ray%t( 3 ) ) ** 2

      ! *** curvature correction ***

      CALL RayNormal( ray%t, ray%phi, ray%c, e1, e2 )  ! Compute ray normals e1 and e2

      ! rotate p-q from e1, e2 system, onto rayn1, rayn2 system

      RotMat( 1, 1 ) = DOT_PRODUCT( rayn1, e1 )
      RotMat( 1, 2 ) = DOT_PRODUCT( rayn1, e2 )
      RotMat( 2, 1 ) = -RotMat( 1, 2 )             ! same as DOT_PRODUCT( rayn2, e1 )
      RotMat( 2, 2 ) = DOT_PRODUCT( rayn2, e2 )

      p_tilde_in = RotMat( 1, 1 ) * ray%p_tilde + RotMat( 1, 2 ) * ray%p_hat
      p_hat_in   = RotMat( 2, 1 ) * ray%p_tilde + RotMat( 2, 2 ) * ray%p_hat

      q_tilde_in = RotMat( 1, 1 ) * ray%q_tilde + RotMat( 1, 2 ) * ray%q_hat
      q_hat_in   = RotMat( 2, 1 ) * ray%q_tilde + RotMat( 2, 2 ) * ray%q_hat

      ! here's the actual curvature change

      p_tilde_out = p_tilde_in + q_tilde_in * R1 - q_hat_in * R2
      p_hat_out   = p_hat_in   + q_tilde_in * R2 + q_hat_in * R3
      !p_hat_out   = p_hat_in   + q_tilde_in * R2 - q_hat_in * R3 ! this one good

      ! rotate p-q back to e1, e2 system (RotMat^(-1) = RotMat^T)

      rayOut%p_tilde = RotMat( 1, 1 ) * p_tilde_out + RotMat( 2, 1 ) * p_hat_out
      rayOut%p_hat   = RotMat( 1, 2 ) * p_tilde_out + RotMat( 2, 2 ) * p_hat_out

      rayOut%q_tilde = ray%q_tilde
      rayOut%q_hat   = ray%q_hat

      !write( *, * ) 'p', ray%p_tilde, ray%p_hat, rayOut%p_tilde, rayOut%p_hat
      !write( *, * ) 'q', ray%q_tilde, ray%q_hat, rayOut%q_tilde, rayOut%q_hat

      ! Logic below fixes a bug when the |dot product| is infinitesimally greater than 1 (then ACos is complex)
      rayOut%phi = ray%phi + 2 * ACOS( MAX( MIN( DOT_PRODUCT( rayn1, e1 ), 1.0D0 ), -1.0D0 ) ) !!!What happens to torsion?
      !write( *, * ) rayn1, e1
      ! write( *, * ) DOT_PRODUCT( rayn1, e1 ), ACOS( MAX( MIN( DOT_PRODUCT( rayn1, e1 ), 1.0D0 ), -1.0D0 ) )
      ! f, g, h continuation; needs curvature corrections
      ! ray3D( is1 )%f    = ray3D( is )%f ! + ray3D( is )%DetQ * R1
      ! ray3D( is1 )%g    = ray3D( is )%g ! + ray3D( is )%DetQ * R1
      ! ray3D( is1 )%h    = ray3D( is )%h ! + ray3D( is )%DetQ * R2
      ! ray3D( is1 )%DetP = ray3D( is )%DetP
      ! ray3D( is1 )%DetQ = ray3D( is )%DetQ

    END SUBROUTINE CurvatureCorrection2

    ! ********************************************************************** !

    SUBROUTINE CalcTangent_Normals( ray3Dt, nBdry, rayt, rayn1, rayn2 )

      ! given the ray tangent vector (not normalized) and the unit normal to the bdry
      ! calculate a unit tangent and the two ray normals

      USE cross_products

      REAL (KIND=8), INTENT( IN )  :: ray3Dt( 3 ), nBdry( 3 )
      REAL (KIND=8), INTENT( OUT ) :: rayt( 3 ), rayn1( 3 ), rayn2( 3 )

      ! incident unit ray tangent and normal
      rayt  = c * ray3Dt                      ! unit tangent to ray
      rayn2 = -cross_product( rayt, nBdry )   ! ray tangent x boundary normal gives refl. plane normal
      rayn2 = rayn2 / NORM2( rayn2 )          ! unit normal to refl. plane
      rayn1 = -cross_product( rayt, rayn2 )   ! ray tangent x refl. plane normal is first ray normal

    END SUBROUTINE CalcTangent_Normals

    ! ********************************************************************** !

    SUBROUTINE ParabotFormulas( z_xx, z_xy, z_yy, nBdry )

      ! analytic formula for the parabolic bottom

      REAL (KIND=8), INTENT( OUT ) :: z_xx, z_xy, z_yy, nBdry( 3 )
      REAL (KIND=8) :: Radius, x, y, z, z_r, z_rr, z_x, z_y, a, b, c, RLen

      IF ( BotTop == 'BOT' ) THEN
         a = 1 / 1000.
         b = 0
         c = 250

         x = ray3D( is )%x( 1 )
         y = ray3D( is )%x( 2 )

         Radius = SQRT( x ** 2 + y ** 2 )   ! radius of seamount at the bounce point

         z = 500.0 * SQRT( 1 + Radius / c )

         z_r = 1 / ( 2 * a * z + b )
         z_x = z_r * x / Radius
         z_y = z_r * y / Radius

         nBdry = [ -z_x, -z_y, 1.D0 ]       !     normal to bdry (outward pointing)
         RLen = NORM2( nBdry )
         nBdry = nBdry / Rlen   ! unit normal to bdry
         ! curvatures

         ! z = z_xx / 2 * x^2 + z_xy * xy + z_yy * y^2 ! coefs are the kappa matrix
         ! RLen is an extra factor based on Eq. 4.4.18 in Cerveny's book
         ! It represents a length along the tangent for a delta x step

         z_rr = -1 / ( 4 * a * a * z * z * z )
         z_xx = ( z_rr * x * x / Radius ** 2 + z_r * y * y / Radius ** 3 ) / RLen
         z_xy = ( z_rr * x * y / Radius ** 2 - z_r * x * y / Radius ** 3 ) / RLen
         z_yy = ( z_rr * y * y / Radius ** 2 + z_r * x * x / Radius ** 3 ) / RLen

         ! write( *, * ) 'x=', x, 'y=', y, 'z_x', z_x, 'z_y', z_y, 'z_r', z_r, 'zrr', z_rr, 'Radius', Radius, 'RLen', RLen 
      END IF
    END SUBROUTINE ParabotFormulas
  END SUBROUTINE Reflect3D
END MODULE Reflect3DMod
