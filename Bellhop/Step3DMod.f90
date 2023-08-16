MODULE Step3DMod

  USE bellhopMod
  USE Bdry3DMod
  USE sspMod
  IMPLICIT NONE
CONTAINS

  SUBROUTINE Step3D( ray0, ray2 )

    ! Does a single step along the ray
    ! x denotes the ray coordinate, ( x, y, z )
    ! t denotes the scaled tangent to the ray (previously (xi, eta, zeta) )
    ! c * t would be the unit tangent

    USE RayNormals

    ! rays
    TYPE( ray3DPt ) :: ray0, ray1, ray2
    INTEGER         :: iSegx0, iSegy0, iSegz0
    REAL  (KIND=8 ) :: gradc0( 3 ), gradc1( 3 ), gradc2( 3 ), &
         c0, cimag0, csq0, cxx0, cyy0, czz0, cxy0, cxz0, cyz0, cnn0, cmn0, cmm0, &
         c1, cimag1, csq1, cxx1, cyy1, czz1, cxy1, cxz1, cyz1, cnn1, cmn1, cmm1, &
         c2, cimag2,       cxx2, cyy2, czz2, cxy2, cxz2, cyz2, c_mat0( 2, 2 ), c_mat1( 2, 2 ), &
         urayt0( 3 ), urayt1( 3 ), h, halfh, hw0, hw1, w0, w1, rho
    REAL  (KIND=8 ) :: gradcjump( 3 ), csjump, cn1jump, cn2jump, tBdry( 3 ), nBdry( 3 ), RM, R1, R2, Tg, Th, &
         rayt( 3 ), rayn1( 3 ), rayn2( 3 )
    REAL  (KIND=8 ) :: e1( 3 ), e2( 3 )                              ! ray normals for ray-centered coordinates
    REAL  (KIND=8 ) :: p_tilde_in(  2 ), p_hat_in(  2 ), q_tilde_in(  2 ), q_hat_in(  2 ), &
         p_tilde_out( 2 ), p_hat_out( 2 ), RotMat( 2, 2 )

    ! The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
    ! to the Heun (second order Runge-Kutta method).
    ! However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

    ! *** Phase 1 (an Euler step)

    !write( *, * ) '______________________________________________________________'
    !write( *, * ) 'in segment ', ISegBoty
    !write( *, * ) 'current coord ', ray0%x

    CALL EvaluateSSP3D( ray0%x, c0, cimag0, gradc0, cxx0, cyy0, czz0, cxy0, cxz0, cyz0, rho, freq, 'TAB' )
    CALL RayNormal( ray0%t, ray0%phi, c0, e1, e2 ) ! Compute ray normals e1 and e2
    CALL Get_c_partials( cxx0, cxy0, cxz0, cyy0, cyz0, czz0, e1, e2, cnn0, cmn0, cmm0 ) ! Compute second partials of c along ray normals

    csq0   = c0 * c0

    iSegx0 = iSegx     ! make note of current layer
    iSegy0 = iSegy
    iSegz0 = iSegz

    h = Beam%deltas       ! initially set the step h, to the basic one, deltas
    urayt0 = c0 * ray0%t  ! unit tangent

    CALL ReduceStep3D( ray0%x, urayt0, iSegx0, iSegy0, iSegz0, h ) ! reduce h to land on boundary

    halfh = 0.5 * h       ! first step of the modified polygon method is a half step

    ray1%x = ray0%x    + halfh * urayt0
    ray1%t = ray0%t    - halfh * gradc0 / csq0

    ! ray1%f    = ray0%f    + halfh * ( c0 * ray0%DetP - cnn0 / csq0 * ray0%DetQ )
    ! ray1%g    = ray0%g    + halfh * ( c0 * ray0%DetP - cmm0 / csq0 * ray0%DetQ )
    ! ray1%h    = ray0%h    + halfh * (                - cmn0 / csq0 * ray0%DetQ )
    ! ray1%DetP = ray0%DetP + halfh / csq0 * ( -cmm0 * ray0%f - cnn0 * ray0%g + 2.0 * cmn0 * ray0%h )
    ! ray1%DetQ = ray0%DetQ + halfh * c0 * ( ray0%f + ray0%g )

    ray1%phi  = ray0%phi  + halfh * ( 1.0 / c0 ) * ray0%t( 3 ) * &
         ( ray0%t( 2 ) * gradc0( 1 ) - ray0%t( 1 ) * gradc0( 2 ) ) / &
         ( ray0%t( 1 ) ** 2 + ray0%t( 2 ) ** 2 )

    c_mat0( 1, : ) = -[ cnn0, cmn0 ] / csq0
    c_mat0( 2, : ) = -[ cmn0, cmm0 ] / csq0

    ray1%p_tilde = ray0%p_tilde + halfh * MATMUL( c_mat0, ray0%q_tilde )
    ray1%q_tilde = ray0%q_tilde + halfh *         c0    * ray0%p_tilde 

    ray1%p_hat   = ray0%p_hat   + halfh * MATMUL( c_mat0, ray0%q_hat )
    ray1%q_hat   = ray0%q_hat   + halfh *         c0    * ray0%p_hat 

    ! *** Phase 2
    
    CALL EvaluateSSP3D( ray1%x, c1, cimag1, gradc1, cxx1, cyy1, czz1, cxy1, cxz1, cyz1, rho, freq, 'TAB' )
    rayt = c1 * ray1%t           ! unit tangent to ray
    CALL RayNormal_unit( rayt, ray1%phi, e1, e2 )
    CALL Get_c_partials( cxx1, cxy1, cxz1, cyy1, cyz1, czz1, e1, e2, cnn1, cmn1, cmm1 ) ! Compute second partials of c along ray normals

    csq1   = c1 * c1
    urayt1 = c1 * ray1%t   ! unit tangent

    CALL ReduceStep3D( ray0%x, urayt1, iSegx0, iSegy0, iSegz0, h ) ! reduce h to land on boundary

    ! use blend of f' based on proportion of a full step used.
    w1  = h / ( 2.0d0 * halfh )
    w0  = 1.0d0 - w1
    hw0 = h * w0
    hw1 = h * w1

    ray2%x = ray0%x    + hw0 * urayt0        + hw1 * urayt1
    ray2%t = ray0%t    - hw0 * gradc0 / csq0 - hw1 * gradc1 / csq1
    !write( *, * ) 'final coord ', ray2%x, w0, w1

    ! ERROR: need to do hw0 and hw1 blend here as well !!!
    ! ray2%f    = ray0%f    + h * ( c1 * ray1%DetP - cnn1 / csq1 * ray1%DetQ )
    ! ray2%g    = ray0%g    + h * ( c1 * ray1%DetP - cmm1 / csq1 * ray1%DetQ )
    ! ray2%h    = ray0%h    + h * (                - cmn1 / csq1 * ray1%DetQ )
    ! ray2%DetP = ray0%DetP + h / csq1 * ( -cmm1 * ray1%f - cnn1 * ray1%g + 2.0 * cmn1 * ray1%h )
    ! ray2%DetQ = ray0%DetQ + h * c1 * ( ray1%f + ray1%g )

    ray2%phi  = ray0%phi  + h * ( 1.0 / c1 ) * ray1%t( 3 ) * &
         ( ray1%t( 2 ) * gradc1( 1 ) - ray1%t( 1 ) * gradc1( 2 ) ) / &
         ( ray1%t( 1 ) ** 2 + ray1%t( 2 ) ** 2 )
    ray2%tau  = ray0%tau + hw0 / CMPLX( c0, cimag0, KIND=8 ) + hw1 / CMPLX( c1, cimag1, KIND=8 )

    ray2%Amp       = ray0%Amp
    ray2%Phase     = ray0%Phase
    ray2%NumTopBnc = ray0%NumTopBnc
    ray2%NumBotBnc = ray0%NumBotBnc

    c_mat1( 1, : ) = -[ cnn1, cmn1 ] / csq1
    c_mat1( 2, : ) = -[ cmn1, cmm1 ] / csq1

    ray2%p_tilde = ray0%p_tilde + hw0 * MATMUL( c_mat0, ray0%q_tilde ) + hw1 * MATMUL( c_mat1, ray1%q_tilde )
    ray2%q_tilde = ray0%q_tilde + hw0 *         c0    * ray0%p_tilde   + hw1 *         c1    * ray1%p_tilde
    ray2%p_hat   = ray0%p_hat   + hw0 * MATMUL( c_mat0, ray0%q_hat   ) + hw1 * MATMUL( c_mat1, ray1%q_hat )
    ray2%q_hat   = ray0%q_hat   + hw0 *         c0    * ray0%p_hat     + hw1 *         c1    * ray1%p_hat 

    ! *** If we crossed an interface, apply jump condition ***

    CALL EvaluateSSP3D( ray2%x, c2, cimag2, gradc2, cxx2, cyy2, czz2, cxy2, cxz2, cyz2, rho, freq, 'TAB' )

    ray2%c = c2

    IF ( iSegx /= iSegx0 .OR. &
         iSegy /= iSegy0 .OR. &
         iSegz /= iSegz0 ) THEN

       gradcjump =  gradc2 - gradc0  ! this is precise only for c-linear layers

!!! what if we cross isegx, isegy, or isegz at the same time?
       IF ( iSegz /= iSegz0 ) THEN
          nBdry = [ 0D0, 0D0, -SIGN( 1.0D0, ray2%t( 3 ) ) ]   ! inward normal to layer
       ELSE IF ( iSegx /= iSegx0 ) THEN
          nBdry = [ -SIGN( 1.0D0, ray2%t( 1 ) ), 0D0, 0D0 ]   ! inward normal to x-segment
       ELSE
          nBdry = [ 0D0, -SIGN( 1.0D0, ray2%t( 2 ) ), 0D0 ]   ! inward normal to y-segment
       END IF
       CALL CurvatureCorrection

    END IF

  CONTAINS
    SUBROUTINE CurvatureCorrection

      ! correct p-q due to jumps in the gradient of the sound speed

      USE cross_products

      Th    = DOT_PRODUCT( ray2%t, nBdry )   ! component of ray tangent, normal to boundary
      tBdry = ray2%t - Th * nBdry            ! tangent, along the boundary, in the reflection plane
      tBdry = tBdry / NORM2( tBdry )         ! unit boundary tangent
      Tg    = DOT_PRODUCT( ray2%t, tBdry )   ! component of ray tangent, along the boundary

      rayt = c2 * ray2%t                     ! unit tangent to ray

      rayn2 = cross_product( rayt, nBdry )   ! ray tangent x boundary normal gives refl. plane normal
      rayn2 = rayn2 / NORM2( rayn2 )         ! unit normal
      rayn1 = -cross_product( rayt, rayn2 )  ! ray tangent x refl. plane normal is first ray normal

      ! normal and tangential derivatives of the sound speed
      cn1jump = DOT_PRODUCT( gradcjump, rayn1 )
      cn2jump = DOT_PRODUCT( gradcjump, rayn2 )
      csjump  = DOT_PRODUCT( gradcjump, rayt  )

      RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
      R1 = RM * ( 2 * cn1jump - RM * csjump ) / c2 ** 2
      R2 = RM * cn2jump / c2 ** 2

      ! *** curvature correction ***

      CALL RayNormal_unit( rayt, ray2%phi, e1, e2 )   ! Compute ray normals e1 and e2

      RotMat( 1, 1 ) = DOT_PRODUCT( rayn1, e1 )
      RotMat( 1, 2 ) = DOT_PRODUCT( rayn1, e2 )
      RotMat( 2, 1 ) = -RotMat( 1, 2 )             ! DOT_PRODUCT( rayn2, e1 )
      RotMat( 2, 2 ) = DOT_PRODUCT( rayn2, e2 )

      ! rotate p-q values in e1, e2 system, onto rayn1, rayn2 system

      p_tilde_in = RotMat( 1, 1 ) * ray2%p_tilde + RotMat( 1, 2 ) * ray2%p_hat
      p_hat_in   = RotMat( 2, 1 ) * ray2%p_tilde + RotMat( 2, 2 ) * ray2%p_hat

      q_tilde_in = RotMat( 1, 1 ) * ray2%q_tilde + RotMat( 1, 2 ) * ray2%q_hat
      q_hat_in   = RotMat( 2, 1 ) * ray2%q_tilde + RotMat( 2, 2 ) * ray2%q_hat

      ! here's the actual curvature change

      p_tilde_out = p_tilde_in - q_tilde_in * R1 - q_hat_in * R2  
      p_hat_out   = p_hat_in   - q_tilde_in * R2

      ! rotate p back to e1, e2 system, q does not change
      ! Note RotMat^(-1) = RotMat^T

      ray2%p_tilde = RotMat( 1, 1 ) * p_tilde_out + RotMat( 2, 1 ) * p_hat_out
      ray2%p_hat   = RotMat( 1, 2 ) * p_tilde_out + RotMat( 2, 2 ) * p_hat_out

    END SUBROUTINE CurvatureCorrection

    !**********************************************************************!

    SUBROUTINE Get_c_partials( cxx, cxy, cxz, cyy, cyz, czz, e1, e2, cnn, cmn, cmm )

      ! Computes the second partials of c along ray normals

      REAL (KIND=8), INTENT( IN  ) :: cxx, cxy, cxz, cyy, cyz, czz  ! curvature of sound speed (cartesian)
      REAL (KIND=8), INTENT( IN  ) :: e1( 3 ), e2( 3 )              ! principal normals
      REAL (KIND=8), INTENT( OUT ) :: cnn, cmn, cmm                 ! curvature of sound speed (ray-centered)

      cnn = cxx * e1( 1 )**2 + cyy * e1( 2 )**2 + czz * e1( 3 )**2 + 2.0 * cxy * e1( 1 ) * e1( 2 ) + &
           2.0 * cxz * e1( 1 ) * e1( 3 ) + 2.0 * cyz * e1( 2 ) * e1( 3 )

      cmn = cxx * e1( 1 ) * e2( 1 ) + cyy * e1( 2 ) * e2( 2 ) + czz * e1( 3 ) * e2( 3 ) + &
           cxy * ( e1( 1 ) * e2( 2 ) + e2( 1 ) * e1( 2 ) ) + cxz * ( e1( 1 ) * e2( 3 ) + e2( 1 ) * e1( 3 ) ) +  &
           cyz * ( e1( 2 ) * e2( 3 ) + e2( 2 ) * e1( 3 ) )

      cmm = cxx * e2( 1 )**2 + cyy * e2( 2 )**2 + czz * e2( 3 )**2 + 2.0 * cxy * e2( 1 ) * e2( 2 ) + &
           2.0 * cxz * e2( 1 ) * e2( 3 ) + 2.0 * cyz * e2( 2 ) * e2( 3 )

      RETURN
    END SUBROUTINE Get_c_partials
  END SUBROUTINE Step3D

  ! **********************************************************************!

  SUBROUTINE ReduceStep3D( x0, urayt, iSegx0, iSegy0, iSegz0, h )

    INTEGER,       INTENT( IN    ) :: iSegx0, iSegy0, iSegz0
    REAL (KIND=8), INTENT( IN    ) :: x0( 3 ), urayt( 3 )  ! ray coordinate and tangent
    REAL (KIND=8), INTENT( INOUT ) :: h                    ! reduced step size
    REAL (KIND=8) :: hInt, hTop, hBot, hxSeg, hySeg, hTopDiag, hBotDiag, hBoxx, hBoxy, hBoxz        ! step sizes
    REAL (KIND=8) :: d( 3 ), d0( 3 ), tri_n( 3 )
    REAL (KIND=8) :: x( 3 )                                ! ray coordinate after full trial step
    REAL (KIND=8) :: xSeg( 2 ), ySeg( 2 )

    ! Detect SSP interface or boundary crossing and reduce step, if necessary, to land on that crossing.
    ! Keep in mind possibility that user put source right on an interface
    ! and that multiple events can occur (crossing interface, top, and bottom in a single step).

    x = x0 + h * urayt ! make a trial step

    ! interface crossing in depth
    ! Step reduction is not done for the top or bottom layer
    ! Instead the SSP is extrapolated
    ! This prevents problems when the boundaries are outside the domain of the SSP
    hInt = huge( hInt )
    IF ( ABS( urayt( 3 ) ) > EPSILON( hInt ) ) THEN
       IF        ( SSP%z( iSegz0     ) > x(  3 ) .AND. iSegz0     > 1  ) THEN
          hInt = ( SSP%z( iSegz0     ) - x0( 3 ) ) / urayt( 3 )
          ! write( *, * ) 'layer crossing', iSegz0, h1
       ELSE IF ( SSP%z( iSegz0 + 1 ) < x(  3 ) .AND. iSegz0 + 1 < SSP%Nz ) THEN
          hInt = ( SSP%z( iSegz0 + 1 ) - x0( 3 ) ) / urayt( 3 )
          ! write( *, * ) 'layer crossing', iSegz0, h1
       END IF
    END IF

    ! top crossing
    hTop = huge( hTop )
    d  = x - Topx              ! vector from top to ray
    IF ( DOT_PRODUCT( Topn, d )  > EPSILON( hTop ) ) THEN
       d0 = x0 - Topx   ! vector from top    node to ray origin
       hTop = -DOT_PRODUCT( d0, Topn ) / DOT_PRODUCT( urayt, Topn )
       ! write( *, * ) 'top crossing'
    END IF

    ! bottom crossing
    hBot = huge( hBot )
    d  = x - Botx              ! vector from bottom to ray
    IF ( DOT_PRODUCT( Botn, d ) > EPSILON( hBot ) ) THEN
       d0 = x0 - Botx   ! vector from bottom node to ray origin
       hBot = -DOT_PRODUCT( d0, Botn ) / DOT_PRODUCT( urayt, Botn )
       ! write( *, * ) 'bottom crossing'
    END IF

    ! top/bottom/ocean segment crossing in x
    hxSeg = huge( hxSeg )
    xSeg( 1 ) = MAX( xTopSeg( 1 ), xBotSeg( 1 ) )
    xSeg( 2 ) = MIN( xTopSeg( 2 ), xBotSeg( 2 ) )
    
    IF ( SSP%Type == 'H' ) THEN   ! ocean segment
       xSeg( 1 ) = MAX( xSeg( 1 ), SSP%Seg%x( iSegx0     ) )
       xSeg( 2 ) = MIN( xSeg( 2 ), SSP%Seg%x( iSegx0 + 1 ) )
    END IF

    IF ( ABS( urayt( 1 ) ) > EPSILON( hxSeg ) ) THEN
       IF          ( x(  1 ) < xSeg( 1 ) ) THEN
          hxSeg = -( x0( 1 ) - xSeg( 1 ) ) / urayt( 1 )
          ! write( *, * ) 'segment crossing in x'
       ELSE IF     ( x(  1 ) > xSeg( 2 ) ) THEN
          hxSeg = -( x0( 1 ) - xSeg( 2 ) ) / urayt( 1 )
          ! write( *, * ) 'segment crossing in x'
       END IF
    END IF

    ! top/bottom/ocean segment crossing in y
    hySeg = huge( hySeg )
    ySeg( 1 ) = MAX( yTopSeg( 1 ), yBotSeg( 1 ) )
    ySeg( 2 ) = MIN( yTopSeg( 2 ), yBotSeg( 2 ) )

    IF ( SSP%Type == 'H' ) THEN   ! ocean segment
       ySeg( 1 ) = MAX( ySeg( 1 ), SSP%Seg%y( iSegy0     ) )
       ySeg( 2 ) = MIN( ySeg( 2 ), SSP%Seg%y( iSegy0 + 1 ) )
    END IF

    IF ( ABS( urayt( 2 ) ) > EPSILON( hySeg ) ) THEN
       IF          ( x(  2 ) < ySeg( 1 ) ) THEN
          hySeg = -( x0( 2 ) - ySeg( 1 ) ) / urayt( 2 )
          !write( *, * ) 'segment crossing in y'
          !write( *, * ) x( 2 ), ySeg( 1 ), ySeg( 2 )
       ELSE IF     ( x(  2 ) > ySeg( 2 ) ) THEN
          hySeg = -( x0( 2 ) - ySeg( 2 ) ) / urayt( 2 )
          !write( *, * ) 'segment crossing in y'
          !write( *, * ) x( 2 ), ySeg( 1 ), ySeg( 2 )
       END IF
    END IF

    ! triangle crossing within a top segment
    hTopDiag = huge( hTopDiag )
    d     = x  - Topx   ! vector from bottom node to ray end
    d0    = x0 - Topx   ! vector from bottom node to ray origin
    tri_n = [ -Top_deltay, Top_deltax, 0.0d0 ]

    IF ( ( DOT_PRODUCT( tri_n, d0 ) > 0.0d0 .AND. DOT_PRODUCT( tri_n, d ) <= 0.0d0 ) .OR. &
         ( DOT_PRODUCT( tri_n, d0 ) < 0.0d0 .AND. DOT_PRODUCT( tri_n, d ) >= 0.0d0 )  ) THEN
       hTopDiag = -DOT_PRODUCT( d0, tri_n ) / DOT_PRODUCT( urayt, tri_n )
       ! write( *, * ) 'diagonal crossing'
    END IF

    ! triangle crossing within a bottom segment
    hBotDiag = huge( hBotDiag )
    d     = x  - Botx   ! vector from bottom node to ray end
    d0    = x0 - Botx   ! vector from bottom node to ray origin
    tri_n = [ -Bot_deltay, Bot_deltax, 0.0d0 ]

    IF ( ( DOT_PRODUCT( tri_n, d0 ) > 0.0d0 .AND. DOT_PRODUCT( tri_n, d ) <= 0.0d0 ) .OR. &
         ( DOT_PRODUCT( tri_n, d0 ) < 0.0d0 .AND. DOT_PRODUCT( tri_n, d ) >= 0.0d0 )  ) THEN
       hBotDiag = -DOT_PRODUCT( d0, tri_n ) / DOT_PRODUCT( urayt, tri_n )
       ! write( *, * ) 'diagonal crossing'
    END IF

    ! ray mask using a box centered at ( 0, xs_3D, ys_3D )
    hBoxx    = huge( hBoxx )
    hBoxy    = huge( hBoxy )
    hBoxz    = huge( hBoxz )
    
    IF ( ABS( x( 1 ) - xs_3D( 1 ) ) > Beam%Box%x ) hBoxx = ( Beam%Box%x - ABS( ( x0( 1 ) - xs_3D( 1 ) ) ) ) / ABS( urayt( 1 ) )
    IF ( ABS( x( 2 ) - xs_3D( 2 ) ) > Beam%Box%y ) hBoxy = ( Beam%Box%y - ABS( ( x0( 2 ) - xs_3D( 2 ) ) ) ) / ABS( urayt( 2 ) )
    IF ( ABS( x( 3 )              ) > Beam%Box%z ) hBoxz = ( Beam%Box%z - ABS(   x0( 3 )                ) ) / ABS( urayt( 3 ) )
    
    h = MIN( h, hInt, hTop, hBot, hxSeg, hySeg, hTopDiag, hBotDiag, hBoxx, hBoxy, hBoxz )  ! take limit set by shortest distance to a crossing

    IF ( h < 1.0d-4 * Beam%deltas ) THEN   ! is it taking an infinitesimal step?
       h = 1.0d-4 * Beam%deltas            ! make sure we make some motion
       iSmallStepCtr = iSmallStepCtr + 1   ! keep a count of the number of sequential small steps
    ELSE
       iSmallStepCtr = 0                   ! didn't do a small step so reset the counter
    END IF

  END SUBROUTINE ReduceStep3D

END MODULE Step3DMod

