
MODULE Cone
  USE bellhopMod
  USE MathConstants
  IMPLICIT NONE

CONTAINS

  SUBROUTINE ConeFormulas2D( z_xx, z_xy, z_yy, nBdry, xs, tradial, xray, BotTop )

    ! analytic formula for the conical seamount
    REAL (KIND=8),     INTENT( IN  ) :: xs( 3 ), tradial( 2 ), xray( 2 )
    CHARACTER (LEN=3), INTENT( IN  ) :: BotTop                  ! Flag indicating bottom or top reflection
    REAL (KIND=8),     INTENT( OUT ) :: z_xx, z_xy, z_yy, nBdry( 2 )
    REAL (KIND=8) :: nBdry3dCone( 3 ), theta, phi, Rray, Radius, x, y   ! for cone reflection

  IF ( BotTop == 'BOT' ) THEN
     phi   = 15 * DegRad   ! 15 degree angle of seamount
     Rray  = xray( 1 )   ! cylindrical range of ray from source

     ! get bearing from cone to ray using (x,y) coordinate of the ray
     x = xs( 1 ) + Rray * tradial( 1 )
     y = xs( 2 ) + Rray * tradial( 2 )
     theta = atan2( y, x )

     nBdry3dCone =  [ -cos( theta ) * sin( phi ), -sin( theta ) * sin( phi ), cos( phi ) ]

     ! nBdry3dCone = -[ -cos( theta ) * tan( phi ), -sin( theta ) * tan( phi ), 1.0D0 ]
     ! nBdry3dCone = -[ ray2D( is )%x( 1 ), ray2D( is )%x( 2 ), -1.0D0  ]
     ! nBdry3dCone = NBdry3dCone / NORM2( NBdry3dCone )

     nBdry( 1 ) = DOT_PRODUCT( nBdry3dCone( 1 : 2 ), tradial )
     nBdry( 2 ) = nBdry3dCone( 3 )

     ! curvatures
     Radius = SQRT( x ** 2 + y ** 2 )   ! radius of seamount at the bounce point

     ! z = z_xx / 2 * x^2 + z_xy * xy + z_yy / 2 * y^2; coefs are the kappa matrix

     z_xx =   y * y / Radius**3 * tan( phi )
     z_xy = - x * y / Radius**3 * tan( phi )
     z_yy =   x * x / Radius**3 * tan( phi )
  END IF
END SUBROUTINE ConeFormulas2D

  SUBROUTINE ConeFormulas3D( z_xx, z_xy, z_yy, nBdry, xs, xray, BotTop )

    ! analytic formula for the conical seamount

    REAL (KIND=8),     INTENT( IN  ) :: xs( 3 ), xray( 3 )
    CHARACTER (LEN=3), INTENT( IN  ) :: BotTop                  ! Flag indicating bottom or top reflection
    REAL (KIND=8),     INTENT( OUT ) :: z_xx, z_xy, z_yy, nBdry( 3 )
    REAL (KIND=8) :: theta, phi, Radius, x, y, RLen

    IF ( BotTop == 'BOT' ) THEN
       phi = 15 * DegRad     ! 15 degree angle of seamount

       x = xray( 1 )
       y = xray( 2 )

       theta = atan2( y, x )   ! bearing from origin (center of cone) to ray
       ! could write the folowing in terms of x, y, and z ...
       nBdry = [ -cos( theta ) * sin( phi ), -sin( theta ) * sin( phi ), cos( phi ) ]   ! unit normal to bdry

       ! curvatures

       Radius = SQRT( x ** 2 + y ** 2 )   ! radius of seamount at the bounce point
       ! z = z_xx / 2 * x^2 + z_xy * xy + z_yy * y^2 ! coefs are the kappa matrix

       RLen = SQRT( 1 + TAN( phi ) ** 2 )   ! 1 + fx^2 + fy^2
       z_xx =   y * y / Radius**3 * tan( phi ) / RLen
       z_xy = - x * y / Radius**3 * tan( phi ) / RLen
       z_yy =   x * x / Radius**3 * tan( phi ) / RLen
    END IF
  END SUBROUTINE ConeFormulas3D

END MODULE CONE

