MODULE sspmod

  ! Holds SSP input by user and associated variables

  ! This module is very similar to the one used by the other programs in the Acoustics Toolbox
  ! However, it returns the SSP *and* its derivatives

  ! Also, a greater premium has been placed on returning this info quickly, since BELLHOP calls it at every step
  ! Therefore more information is pre-computed

  USE FatalError
  IMPLICIT NONE

  SAVE
  INTEGER, PARAMETER, PRIVATE :: ENVFile = 5, PRTFile = 6
  INTEGER, PARAMETER     :: MaxSSP = 100001
  INTEGER                :: iSegr = 1, iSegx = 1, iSegy = 1, iSegz = 1
  INTEGER,       PRIVATE :: iz
  REAL (KIND=8), PRIVATE :: Depth, W
  REAL (KIND=8)          :: zTemp, betaPowerLaw = 1, fT = 1D20
  REAL (KIND=8)          :: alphaR = 1500, betaR = 0, alphaI = 0, betaI = 0, rhoR = 1

  TYPE rxyz_vector
    REAL(KIND=8), ALLOCATABLE :: r(:), x(:), y(:), z(:)
  END TYPE rxyz_vector

  ! SSP
  TYPE SSPStructure
    INTEGER              :: NPts, Nr, Nx, Ny, Nz
    REAL    (KIND=8)     :: z( MaxSSP ), rho( MaxSSP )
    COMPLEX (KIND=8)     :: c( MaxSSP ), cz( MaxSSP ), n2( MaxSSP ), n2z( MaxSSP ), cSpline( 4, MaxSSP )
    COMPLEX (KIND=8)     :: cCoef( 4, MaxSSP ), CSWork( 4, MaxSSP )   ! for PCHIP coefs.
    REAL (KIND=8), ALLOCATABLE :: cMat( :, : ), czMat( :, : ), cMat3( :, :, : ), czMat3( :, :, : )
    TYPE ( rxyz_vector ) :: Seg
    CHARACTER (LEN=1)    :: Type
    CHARACTER (LEN=2)    :: AttenUnit
  END TYPE SSPStructure

  TYPE( SSPStructure ) :: SSP

  ! *** Halfspace properties structure ***

  TYPE HSInfo
     REAL     (KIND=8) :: alphaR, alphaI, betaR, betaI  ! compressional and shear wave speeds/attenuations in user units
     COMPLEX  (KIND=8) :: cP, cS                 ! P-wave, S-wave speeds
     REAL     (KIND=8) :: rho, Depth             ! density, depth
     CHARACTER (LEN=1) :: BC                     ! Boundary condition type
     CHARACTER (LEN=6) :: Opt
  END TYPE HSInfo

  TYPE BdryPt
     TYPE( HSInfo )   :: HS
  END TYPE

  TYPE BdryType
     TYPE( BdryPt )   :: Top, Bot
  END TYPE BdryType

  TYPE(BdryType) :: Bdry

CONTAINS

  SUBROUTINE EvaluateSSP( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )

    ! Call the particular profil routine indicated by the SSP%Type and perform Task
    !   Task = 'TAB'  then tabulate cp, cs, rhoT 
    !   Task = 'INI' then initialize

    REAL      (KIND=8), INTENT( IN  ) :: freq
    REAL      (KIND=8), INTENT( IN  ) :: x( 2 )      ! r-z coordinate where SSP is to be evaluated
    CHARACTER ( LEN=3), INTENT( IN  ) :: Task
    REAL      (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
    REAL      (KIND=8)                :: gradc_3d( 3 ), cxx, cyy, cxy, cxz, cyz
    REAL      (KIND=8)                :: x3( 3 )

    SELECT CASE ( SSP%Type )
    CASE ( 'N' )  !  N2-linear profile option
       CALL n2Linear( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )
    CASE ( 'C' )  !  C-linear profile option
       CALL cLinear(  x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )
    CASE ( 'P' )  !  monotone PCHIP ACS profile option
       CALL cPCHIP(   x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )
    CASE ( 'S' )  !  Cubic spline profile option
       CALL cCubic(   x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )
    CASE ( 'Q' )
       CALL Quad(     x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )
    CASE ( 'H' )
       ! this is called by BELLHOP3D only once, during READIN
       ! possibly the logic should be changed to call EvaluateSSP2D or 3D
       x3 = [ 0.0D0, 0.0D0, x( 2 ) ]
       CALL Hexahedral( x3, c, cimag, gradc_3d, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, Task )
    CASE ( 'A' )  !  Analytic profile option
       CALL Analytic( x, c, cimag, gradc, crr, crz, czz, rho )
    CASE DEFAULT
       WRITE( PRTFile, * ) 'Profile option: ', SSP%Type
       CALL ERROUT( 'BELLHOP: EvaluateSSP', 'Invalid profile option' )
    END SELECT

  END SUBROUTINE EvaluateSSP
  
!**********************************************************************!

SUBROUTINE EvaluateSSP2D( x2D, c, cimag, gradc, crr, crz, czz, rho, xs, tradial, freq )

  ! Converts cartesian gradients to polar
  ! Called from BELLHOP3D to get a 2D slice out of the 3D SSP

  REAL (KIND=8), INTENT( IN  ) :: x2D( 2 ), xs( 3 ), tradial( 2 ), freq
  REAL (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 2 ), czz, crz, crr, rho
  REAL (KIND=8)                :: x( 3 ), gradc3D(3 ), cxx, cyy, cxy, cxz, cyz

  ! convert polar coordinate to cartesian
  x = [ xs( 1 ) + x2D( 1 ) * tradial( 1 ), xs( 2 ) + x2D( 1 ) * tradial( 2 ), x2D( 2 ) ]

  CALL EvaluateSSP3D( x, c, cimag, gradc3D, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, 'TAB' )

  gradc( 1 )  = DOT_PRODUCT( tradial, gradc3D( 1 : 2 ) )  ! r derivative
  gradc( 2 )  = gradc3D( 3 )                              ! z derivative

  crz = tradial( 1 ) * cxz + tradial( 2 ) * cyz
  crr = cxx * ( tradial( 1 ) )**2 + 2.0 * cxy * tradial( 1 ) * tradial( 2 ) + cyy * ( tradial( 2 ) )**2

  RETURN
END SUBROUTINE EvaluateSSP2D

  !**********************************************************************!

  SUBROUTINE EvaluateSSP3D( x, c, cimag, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, Task )

    ! Call the particular profil routine indicated by the SSP%Type and perform Task
    !   Task = 'TAB'  then tabulate cp, cs, rhoT 
    !   Task = 'INI' then initialize

    REAL      (KIND=8), INTENT( IN  ) :: freq
    REAL      (KIND=8), INTENT( IN  ) :: x( 3 )      ! x-y-z coordinate where SSP is to be evaluated
    CHARACTER ( LEN=3), INTENT( IN  ) :: Task
    REAL      (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 3 ), cxx, cyy, czz, cxy, cxz, cyz, rho
    REAL      (KIND=8)                :: x_rz( 2 ), gradc_rz( 2 ), crr, crz

    x_rz = [ 0.0D0, x( 3 ) ]   ! convert x-y-z coordinate to cylindrical coordinate

    SELECT CASE ( SSP%Type )
    CASE ( 'N' )
       CALL n2Linear( x_rz, c, cimag, gradc_rz, crr, crz, czz, rho, freq, Task )
    CASE ( 'C' )
       CALL cLinear(  x_rz, c, cimag, gradc_rz, crr, crz, czz, rho, freq, Task )
    CASE ( 'S' )
       CALL cCubic(   x_rz, c, cimag, gradc_rz, crr, crz, czz, rho, freq, Task )
    CASE ( 'H' )
       CALL Hexahedral( x, c, cimag, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, Task )
    CASE ( 'A' )
       CALL Analytic3D( x, c, cimag, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho )
    CASE DEFAULT
       WRITE( PRTFile, * ) 'Profile option: ', SSP%Type
       CALL ERROUT( 'BELLHOP3D: EvaluateSSP3D', 'Invalid profile option' )
    END SELECT

    SELECT CASE ( SSP%Type )
    CASE ( 'N', 'C', 'S' )
       gradc = [ 0.0D0, 0.0D0, gradc_rz( 2 ) ]

       cxx   = 0.0D0
       cyy   = 0.0D0
       cxy   = 0.0D0
       cxz   = 0.0D0
       cyz   = 0.0D0
    END SELECT

  END SUBROUTINE EvaluateSSP3D

  !**********************************************************************!

  SUBROUTINE n2Linear( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )

    ! N2-linear interpolation of SSP data

    REAL     (KIND=8), INTENT( IN  ) :: freq
    REAL     (KIND=8), INTENT( IN  ) :: x( 2 )   ! r-z coordinate where sound speed is evaluated
    CHARACTER (LEN=3), INTENT( IN  ) :: Task
    REAL     (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho ! sound speed and its derivatives
    
    IF ( Task == 'INI' ) THEN   ! read in SSP data
       Depth     = x( 2 )
       CALL ReadSSP( Depth, freq )
              
       SSP%n2(  1 : SSP%NPts ) = 1.0 / SSP%c( 1 : SSP%NPts ) ** 2

       ! compute gradient, n2z
       DO iz = 2, SSP%Npts
          SSP%n2z( iz - 1 ) = ( SSP%n2(   iz ) - SSP%n2(   iz - 1 ) ) / &
                              ( SSP%z(    iz ) - SSP%z(    iz - 1 ) )
       END DO
    ELSE                         ! return SSP info

       IF ( x( 2 ) < SSP%z( iSegz ) .OR. x( 2 ) > SSP%z( iSegz + 1 ) ) THEN
          DO iz = 2, SSP%NPts   ! Search for bracketing Depths
             IF ( x( 2 ) < SSP%z( iz ) ) THEN
                iSegz = iz - 1
                EXIT
             END IF
          END DO
       END IF

       W = ( x( 2 ) - SSP%z( iSegz ) ) / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )

       c     = REAL(  1.0D0 / SQRT( ( 1.0D0 - W ) * SSP%n2( iSegz ) + W * SSP%n2( iSegz + 1 ) ) )
       cimag = AIMAG( 1.0D0 / SQRT( ( 1.0D0 - W ) * SSP%n2( iSegz ) + W * SSP%n2( iSegz + 1 ) ) )

       gradc = [ 0.0D0, -0.5D0 * c * c * c * REAL( SSP%n2z( iSegz ) ) ]
       crr   = 0.0d0
       crz   = 0.0d0
       czz   = 3.0d0 * gradc( 2 ) * gradc( 2 ) / C

       rho   = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )
    END IF

  END SUBROUTINE n2Linear

  !**********************************************************************!

  SUBROUTINE cLinear( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )

    ! c-linear interpolation of SSP data

    REAL     (KIND=8), INTENT( IN  ) :: freq
    REAL     (KIND=8), INTENT( IN  ) :: x( 2 )   ! r-z coordinate where sound speed is evaluated
    CHARACTER (LEN=3), INTENT( IN  ) :: Task
    REAL     (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho ! sound speed and its derivatives
    
    IF ( Task == 'INI' ) THEN   ! read in SSP data
       Depth     = x( 2 )
       CALL ReadSSP( Depth, freq )
    ELSE                        ! return SSP info

       IF ( x( 2 ) < SSP%z( iSegz ) .OR. x( 2 ) > SSP%z( iSegz + 1 ) ) THEN
          DO iz = 2, SSP%NPts   ! Search for bracketing Depths
             IF ( x( 2 ) < SSP%z( iz ) ) THEN
                iSegz = iz - 1
                EXIT
             END IF
          END DO
       END IF

       c     = REAL(  SSP%c( iSegz ) + ( x( 2 ) - SSP%z( iSegz ) ) * SSP%cz( iSegz ) )
       cimag = AIMAG( SSP%c( iSegz ) + ( x( 2 ) - SSP%z( iSegz ) ) * SSP%cz( iSegz ) )
       gradc = [ 0.0D0, REAL( SSP%cz( iSegz ) ) ]
       crr   = 0.0d0
       crz   = 0.0d0
       czz   = 0.0d0

       W     = ( x( 2 ) - SSP%z( iSegz ) ) / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )
       rho   = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )
    END IF

  END SUBROUTINE cLinear

  !**********************************************************************!

  SUBROUTINE cPCHIP( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )

    ! This implements the new monotone piecewise cubic Hermite interpolating
    ! polynomial (PCHIP) algorithm for the interpolation of the sound speed c.

    USE pchipMod
    REAL     (KIND=8), INTENT( IN  ) :: freq
    REAL     (KIND=8), INTENT( IN  ) :: x( 2 )   ! r-z coordinate where sound speed is evaluated
    CHARACTER (LEN=3), INTENT( IN  ) :: Task
    REAL     (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho ! sound speed and its derivatives
    REAL     (KIND=8) :: xt
    COMPLEX  (KIND=8) :: c_cmplx

    IF ( Task == 'INI' ) THEN   ! read in SSP data

       Depth     = x( 2 )
       CALL ReadSSP( Depth, freq )

       !                                                               2      3
       ! compute coefficients of std cubic polynomial: c0 + c1*x + c2*x + c3*x
       !

       CALL PCHIP( SSP%z, SSP%c, SSP%NPts, SSP%cCoef, SSP%CSWork )

    ELSE                        ! return SSP info

       IF ( x( 2 ) < SSP%z( iSegz ) .OR. x( 2 ) > SSP%z( iSegz + 1 ) ) THEN
          DO iz = 2, SSP%NPts   ! Search for bracketing Depths
             IF ( x( 2 ) < SSP%z( iz ) ) THEN
                iSegz = iz - 1
                EXIT
             END IF
          END DO
       END IF

       xt = x( 2 ) - SSP%z( iSegz )
       c_cmplx = SSP%cCoef( 1, iSegz ) &
             + ( SSP%cCoef( 2, iSegz ) &
             + ( SSP%cCoef( 3, iSegz ) &
             +   SSP%cCoef( 4, iSegz ) * xt ) * xt ) * xt

       c     = REAL(  c_cmplx )
       cimag = AIMAG( c_cmplx )

       gradc = [ 0.0D0, REAL( SSP%cCoef( 2, iSegz ) &
                  + ( 2.0D0 * SSP%cCoef( 3, iSegz ) &
                    + 3.0D0 * SSP%cCoef( 4, iSegz ) * xt ) * xt ) ]

       crr   = 0.0D0
       crz   = 0.0D0
       czz   = REAL( 2.0D0 * SSP%cCoef( 3, iSegz ) + 6.0D0 * SSP%cCoef( 4, iSegz ) * xt )

       W     = ( x( 2 ) - SSP%z( iSegz ) ) / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )
       rho   = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )   ! linear interp of density

    END IF

  END SUBROUTINE cPCHIP

  !**********************************************************************!

  SUBROUTINE cCubic( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )

    ! Cubic spline interpolation

    REAL     (KIND=8), INTENT( IN )  :: freq
    REAL     (KIND=8), INTENT( IN  ) :: x( 2 )   ! r-z coordinate where sound speed is evaluated
    CHARACTER (LEN=3), INTENT( IN  ) :: Task
    REAL     (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho ! sound speed and its derivatives
    INTEGER                          :: iBCBeg, iBCEnd
    REAL     (KIND=8)                :: hSpline
    COMPLEX  (KIND=8)                :: c_cmplx, cz_cmplx, czz_cmplx
    
    IF ( Task == 'INI' ) THEN

       ! *** Task 'INIT' for initialization ***

       Depth     = x( 2 )
       CALL ReadSSP( Depth, freq )

       SSP%cSpline( 1, 1 : SSP%NPts ) = SSP%c( 1 : SSP%NPts )
       
       ! Compute spline coefs
       iBCBeg = 0
       iBCEnd = 0
       CALL CSpline( SSP%z, SSP%cSpline(   1, 1 ), SSP%NPts, iBCBeg, iBCEnd, SSP%NPts )
    ELSE

       ! *** Section to return SSP info ***

       IF ( x( 2 ) < SSP%z( iSegz ) .OR. x( 2 ) > SSP%z( iSegz + 1 ) ) THEN
          DO iz = 2, SSP%NPts   ! Search for bracketing Depths
             IF ( x( 2 ) < SSP%z( iz ) ) THEN
                iSegz = iz - 1
                EXIT
             END IF
          END DO

       END IF

       hSpline = x( 2 ) - SSP%z( iSegz )

       ! c   = Spline(   SSP%cSpline( 1, iSegz ), hSpline )
       ! cz  = SplineX(  SSP%cSpline( 1, iSegz ), hSpline )
       ! czz = SplineXX( SSP%cSpline( 1, iSegz ), hSpline )

       CALL SplineALL( SSP%cSpline( 1, iSegz ), hSpline, c_cmplx, cz_cmplx, czz_cmplx )

       c     = DBLE(  c_cmplx )
       cimag = AIMAG( c_cmplx )
       gradc = [ 0.0D0, DBLE( cz_cmplx ) ]
       czz   = DBLE( czz_cmplx )
       crr   = 0.0d0
       crz   = 0.0d0

       ! linear interpolation for density
       W   = ( x( 2 ) - SSP%z( iSegz ) ) / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )
       rho = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )
    END IF

  END SUBROUTINE cCubic

  !**********************************************************************!

  SUBROUTINE Quad( x, c, cimag, gradc, crr, crz, czz, rho, freq, Task )

    ! Bilinear quadrilatteral interpolation of SSP data in 2D

    INTEGER,           PARAMETER      :: SSPFile = 40
    REAL      (KIND=8), INTENT( IN  ) :: freq
    REAL      (KIND=8), INTENT( IN  ) :: x( 2 )   ! r-z coordinate where sound speed is evaluated
    CHARACTER (LEN=3),  INTENT( IN  ) :: Task
    REAL      (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho ! sound speed and its derivatives
    INTEGER                           :: AllocateStatus, iSegT, iz2
    REAL      (KIND=8)                :: c1, c2, cz1, cz2, cr, cz, s1, s2, delta_r, delta_z
    
    IF ( Task == 'INI' ) THEN

       !  *** read in SSP data ***

       Depth = x( 2 )
       CALL ReadSSP( Depth, freq )

      ! Read the 2D SSP matrix
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Using range-dependent sound speed'

       READ( SSPFile,  * ) SSP%Nr
       WRITE( PRTFile, * ) 'Number of SSP ranges = ', SSP%Nr

       IF ( SSP%Nr < 2 ) THEN
          CALL ERROUT( 'sspMod: Quad', 'You must have a least two profiles in your 2D SSP field'  )
       END IF

       ALLOCATE( SSP%cMat( SSP%NPts, SSP%Nr ), SSP%czMat( SSP%NPts - 1, SSP%Nr ), SSP%Seg%r( SSP%Nr ), STAT = AllocateStatus )
       IF ( AllocateStatus /= 0 ) CALL ERROUT( 'READIN: Quad', 'Insufficient memory to store SSP'  )

       READ( SSPFile,  * ) SSP%Seg%r( 1 : SSP%Nr )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Profile ranges (km):'
       WRITE( PRTFile, FMT="( F10.2 )"  ) SSP%Seg%r( 1 : SSP%Nr )

       SSP%Seg%r = 1000.0 * SSP%Seg%r   ! convert km to m

       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Sound speed matrix:'
       WRITE( PRTFile, * ) ' Depth (m )     Soundspeed (m/s)'
       DO iz2 = 1, SSP%NPts
          READ(  SSPFile, * ) SSP%cMat( iz2, : )
          ! WRITE( PRTFile, FMT="( 'iSegz depth = ', F10.2, ' m' )"  ) SSP%z( iz2 )
          ! WRITE( PRTFile, FMT="( 12F10.2 )"  ) SSP%cMat( iz2, : )
          WRITE( PRTFile, FMT="( 12F10.2 )"  ) SSP%z( iz2 ), SSP%cMat( iz2, : )
       END DO

       CLOSE( SSPFile )

       ! calculate cz
       DO iSegt = 1, SSP%Nr
          DO iz2 = 2, SSP%NPts
             delta_z = ( SSP%z( iz2 ) - SSP%z( iz2 - 1 ) )
             SSP%czMat( iz2 - 1, iSegt ) = ( SSP%cMat( iz2, iSegt ) - SSP%cMat( iz2 - 1, iSegt ) ) / delta_z
          END DO
       END DO

       SSP%Nz = SSP%NPts
       RETURN

    ELSE

       ! *** Section to return SSP info ***

       ! check depth-layer contains x( 2 ) in [ SSP%z( iSegz ), SSP%z( iSegz + 1 ) ]
       IF ( x( 2 ) < SSP%z( iSegz ) .OR. x( 2 ) > SSP%z( iSegz + 1 ) ) THEN
          DO iz = 2, SSP%NPts   ! Search for bracketing Depths
             IF ( x( 2 ) < SSP%z( iz ) ) THEN
                iSegz = iz - 1
                EXIT
             END IF
          END DO
       END IF

       ! The following tries to be more efficient than the code above by searching away from the current layer
       ! rather than searching through all the layers
       ! However, seems to be no faster
       ! Also, this code caused a problem on at/tests/Gulf for the range-dep. test cases
!!$     IF ( x( 2 ) < SSP%z( iSegz ) .AND. iSegz > 1 ) THEN
!!$        DO iz = iSegz - 1, 1, -1   ! Search for bracketing Depths
!!$           IF ( x( 2 ) > SSP%z( iz ) ) THEN
!!$              iSegz = iz
!!$              EXIT
!!$           END IF
!!$        END DO
!!$     END IF
!!$
!!$     IF ( x( 2 ) > SSP%z( iSegz + 1 ) .AND. iSegz < SSP%NPts - 2 ) THEN
!!$        DO iz = iSegz + 2, SSP%NPts   ! Search for bracketing Depths
!!$           IF ( x( 2 ) < SSP%z( iz ) ) THEN
!!$              iSegz = iz - 1
!!$              EXIT
!!$           END IF
!!$        END DO
!!$     END IF

       ! Check that x is inside the box where the sound speed is defined
       IF ( x( 1 ) < SSP%Seg%r( 1 ) .OR. x( 1 ) > SSP%Seg%r( SSP%Nr ) ) THEN ! .OR. &
          WRITE( PRTFile, * ) 'ray is outside the box where the ocean soundspeed is defined'
          WRITE( PRTFile, * ) ' x = ( r, z ) = ', x
          CALL ERROUT( 'sspMod: Quad', 'ray is outside the box where the soundspeed is defined' )
       END IF

       ! check range-segment contains x( 1 ) in [ SSP%Seg%r( iSSP%Seg ), SSP%Seg%r( iSegr + 1 ) )
       IF ( x( 1 ) < SSP%Seg%r( iSegr ) .OR. x( 1 ) >= SSP%Seg%r( iSegr + 1 ) ) THEN
          DO iSegT = 2, SSP%Nr   ! Search for bracketing segment ranges
             IF ( x( 1 ) < SSP%Seg%r( iSegT ) ) THEN
                iSegr = iSegT - 1
                EXIT
             END IF
          END DO
       END IF

       ! for this depth, x( 2 ), get the sound speed at both ends of the segment
       cz1 = SSP%czMat( iSegz, iSegr )
       cz2 = SSP%czMat( iSegz, iSegr + 1 )

       s2      = x( 2 ) - SSP%z( iSegz )
       delta_z = SSP%z( iSegz + 1 ) - SSP%z( iSegz )
       
       c1 = SSP%cMat( iSegz, iSegr     ) + s2 * cz1
       c2 = SSP%cMat( iSegz, iSegr + 1 ) + s2 * cz2

       ! s1 = proportional distance of x( 1 ) in range
       delta_r = ( SSP%Seg%r( iSegr + 1 ) - SSP%Seg%r( iSegr ) )
       s1 = ( x( 1 ) - SSP%Seg%r( iSegr ) ) / delta_r
       s1 = MIN( s1, 1.0D0 )   ! force piecewise constant extrapolation for points outside the box
       s1 = MAX( s1, 0.0D0 )   ! "

       c     = ( 1.0D0 - s1 ) * c1  + s1 * c2

       ! interpolate the attenuation !!!! This will use the wrong segment if the ssp in the envil is sampled at different depths
       s2    = s2 / delta_z   ! convert to a proportional depth in the layer
       cimag = AIMAG( ( 1.0D0 - s2 ) * SSP%c( Isegz )  + s2 * SSP%c( Isegz + 1 ) )   ! volume attenuation is taken from the single c(z) profile

       cz  = ( 1.0D0 - s1 ) * cz1 + s1 * cz2

       cr  = ( c2  - c1  ) / delta_r
       crz = ( cz2 - cz1 ) / delta_r

       gradc = [ cr, cz ]
       crr   = 0.0D0
       czz   = 0.0D0

       ! linear interpolation for density
       W   = ( x( 2 ) - SSP%z( iSegz ) ) / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )
       rho = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )
    END IF

  END SUBROUTINE Quad

  !**********************************************************************!

  SUBROUTINE Hexahedral( x, c, cimag, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, Task )

    ! Trilinear hexahedral interpolation of SSP data in 3D
    ! assumes a rectilinear case (not the most general hexahedral)

    INTEGER,            PARAMETER     :: SSPFile = 40
    REAL      (KIND=8), INTENT( IN  ) :: freq
    REAL      (KIND=8), INTENT( IN  ) :: x( 3 )   ! x-y-z coordinate where sound speed is evaluated
    CHARACTER (LEN =3), INTENT( IN  ) :: Task
    REAL      (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 3 ), cxx, cyy, czz, cxy, cxz, cyz, rho ! sound speed and its derivatives
    INTEGER                           :: AllocateStatus, iSegxt, iSegyt, iy2, iz2, iSegxTT( 1 ), iSegyTT( 1 ), iSegzTT( 1 )
    REAL      (KIND=8)                :: c1, c2, c11, c12, c21, c22, cz11, cz12, cz21, cz22, cz1, cz2, &
                                         cx, cy, cz, s1, s2, s3

    IF ( Task == 'INI' ) THEN

       !  *** Section to read in SSP data ***

       ! Read dummy SSP information from the environmental file
       ! This is over-ridden by the info in the SSP file
       ! However, cz info is used in beam selection
       Depth = x( 3 )
       CALL ReadSSP( Depth, freq )

       ! Read the 3D SSP matrix

       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Reading sound speed profile from file'

       ! x coordinates
       READ( SSPFile,  * ) SSP%Nx
       WRITE( PRTFile, * ) 'Number of points in x = ', SSP%Nx
       ALLOCATE( SSP%Seg%x( SSP%Nx ), STAT = AllocateStatus )
       IF ( AllocateStatus /= 0 ) CALL ERROUT( 'sspMod: Hexahedral', 'Insufficient memory to store SSP'  )
       READ( SSPFile,  * ) SSP%Seg%x
       !WRITE( PRTFile, * )
       !WRITE( PRTFile, * ) 'x-coordinates of SSP (km):'
       !WRITE( PRTFile, FMT="( F10.2 )"  ) SSP%Seg%x( 1 : SSP%Nx )

       ! y coordinates
       READ( SSPFile,  * ) SSP%Ny
       WRITE( PRTFile, * ) 'Number of points in y = ', SSP%Ny
       ALLOCATE( SSP%Seg%y( SSP%Ny ), STAT = AllocateStatus )
       IF ( AllocateStatus /= 0 ) CALL ERROUT( 'sspMod: Hexahedral', 'Insufficient memory to store SSP'  )
       READ( SSPFile,  * ) SSP%Seg%y
       !WRITE( PRTFile, * )
       !WRITE( PRTFile, * ) 'y-coordinates of SSP (km):'
       !WRITE( PRTFile, FMT="( F10.2 )"  ) SSP%Seg%y( 1 : SSP%Ny )

       ! z coordinates
       READ( SSPFile,  * ) SSP%Nz
       WRITE( PRTFile, * ) 'Number of points in z = ', SSP%Nz
       ALLOCATE( SSP%Seg%z( SSP%Nz ), STAT = AllocateStatus )
       IF ( AllocateStatus /= 0 ) CALL ERROUT( 'sspMod: Hexahedral', 'Insufficient memory to store SSP'  )
       READ( SSPFile,  * ) SSP%Seg%z
       !WRITE( PRTFile, * )
       !WRITE( PRTFile, * ) 'z-coordinates of SSP (km):'
       !WRITE( PRTFile, FMT="( F10.2 )"  ) SSP%Seg%z( 1 : SSP%Nz )

       ! SSP matrix should be bigger than 2x2x2
       IF ( SSP%Nx < 2 .OR. SSP%Ny < 2 .OR. SSP%Nz < 2 ) THEN
          CALL ERROUT( 'READIN: Hexahedral', &
               'You must have a least two points in x, y, z directions in your 3D SSP field'  )
       END IF

       ALLOCATE( SSP%cMat3( SSP%Nx, SSP%Ny, SSP%Nz ), SSP%czMat3( SSP%Nx, SSP%Ny, SSP%Nz - 1 ), STAT = AllocateStatus )
       IF ( AllocateStatus /= 0 ) CALL ERROUT( 'sspMod: Hexahedral', 'Insufficient memory to store SSP'  )

       WRITE( PRTFile, * )
       ! WRITE( PRTFile, * ) 'Sound speed matrix:'
       DO iz2 = 1, SSP%Nz
          DO iy2 = 1, SSP%Ny
             READ(  SSPFile, * ) SSP%cMat3( :, iy2, iz2 )
             ! WRITE( PRTFile, FMT="( 'z-section = ', F10.2, '  m' )"  ) SSP%Seg%z( iz2 )
             ! WRITE( PRTFile, FMT="( 'y-section = ', F10.2, ' km' )"  ) SSP%Seg%y( iy2 )
             ! WRITE( PRTFile, FMT="( 12F10.2 )"  ) SSP%cMat3( :, iy2, iz2 )
          END DO
       END DO

       CLOSE( SSPFile )

       SSP%Seg%x = 1000.0 * SSP%Seg%x   ! convert km to m
       SSP%Seg%y = 1000.0 * SSP%Seg%y   ! convert km to m

       ! calculate cz
       DO iSegxt = 1, SSP%Nx
          DO iSegyt = 1, SSP%Ny
             DO iz2 = 2, SSP%Nz
                SSP%czMat3( iSegxt, iSegyt, iz2 - 1 ) = &
                     ( SSP%cMat3( iSegxt, iSegyt, iz2 ) - SSP%cMat3( iSegxt, iSegyt, iz2 - 1 ) ) / &
                     ( SSP%Seg%z(                 iz2 ) - SSP%Seg%z(                 iz2 - 1 ) )
             END DO
          END DO
       END DO

       ! over-ride the SSP%z values read in from the environmental file with these new values
       SSP%NPts = SSP%Nz
       SSP%z( 1 : SSP%Nz ) = SSP%Seg%z

       RETURN
    ELSE

       ! *** Section to return SSP info ***

       ! Check that x is inside the box where the sound speed is defined
       IF ( x( 1 ) < SSP%Seg%x( 1 ) .OR. x( 1 ) > SSP%Seg%x( SSP%Nx ) .OR. &
            x( 2 ) < SSP%Seg%y( 1 ) .OR. x( 2 ) > SSP%Seg%y( SSP%Ny ) ) THEN ! .OR. &
!!$          x( 3 ) < SSP%Seg%z( 1 ) .OR. x( 3 ) > SSP%Seg%z( SSP%Nz ) ) THEN
          WRITE( PRTFile, * ) 'ray is outside the box where the ocean soundspeed is defined'
          WRITE( PRTFile, * ) ' x = ( x, y, z ) = ', x
          CALL ERROUT( 'sspMod: Hexahedral', 'ray is outside the box where the soundspeed is defined' )
       END IF

       ! check x-segment contains x( 1 ) in [ SSP%Seg%x( iSegx ), SSP%Seg%x( iSegx + 1 ) )
       !IF ( x( 1 ) <  SSP%Seg%x( iSegx     ) ) iSegx = MAX( 1,          iSegx - 1 )   ! bump left
       !IF ( x( 1 ) >= SSP%Seg%x( iSegx + 1 ) ) iSegx = MIN( SSP%Nx - 1, iSegx + 1 )   ! bump right

       IF ( x( 1 ) < SSP%Seg%x( iSegx ) .OR. x( 1 ) >= SSP%Seg%x( iSegx + 1 ) ) THEN
!!$          DO iSegxT = 2, SSP%Nx   ! Search for bracketing segment ranges
!!$             IF ( x( 1 ) < SSP%Seg%x( iSegxT ) ) THEN
!!$                iSegx = iSegxT - 1
!!$                EXIT
!!$             END IF
!!$          END DO

          iSegxTT = MAXLOC( SSP%Seg%x( : ), SSP%Seg%x( : ) < x( 1 ) )

          IF ( iSegxTT( 1 ) > 0 .AND. iSegxTT( 1 ) < SSP%Nx ) THEN  ! iSegx MUST LIE IN [ 1, SSP%Nx - 1 ]
             iSegx = iSegxTT( 1 )   
          END IF

       END IF

       ! check y-segment contains x( 2 ) in [ SSP%Seg%y( iSegy ), SSP%Seg%y( iSegy + 1 ) )
       IF ( x( 2 ) < SSP%Seg%y( iSegy ) .OR. x( 2 ) >= SSP%Seg%y( iSegy + 1 ) ) THEN
!!$          DO iSegyT = 2, SSP%Ny   ! Search for bracketing segment ranges
!!$             IF ( x( 2 ) < SSP%Seg%y( iSegyT ) ) THEN
!!$                iSegy = iSegyT - 1
!!$                EXIT
!!$             END IF
!!$          END DO

          iSegyTT = MAXLOC( SSP%Seg%y( : ), SSP%Seg%y( : ) < x( 2 ) )

          IF ( iSegyTT( 1 ) > 0 .AND. iSegyTT( 1 ) < SSP%Ny ) THEN  ! iSegx MUST LIE IN [ 1, SSP%Ny - 1 ]
             iSegy = iSegyTT( 1 )   
          END IF

       END IF

       ! check depth-layer contains x( 3 ) in [ SSP%Seg%z( iSegz ), SSP%Seg%z( iSegz + 1 ) ]
       IF ( x( 3 ) < SSP%Seg%z( iSegz ) .OR. x( 3 ) > SSP%Seg%z( iSegz + 1 ) ) THEN
!!$          DO iz = 2, SSP%Nz   ! Search for bracketing Depths
!!$             IF ( x( 3 ) < SSP%Seg%z( iz ) ) THEN
!!$                iSegz = iz - 1
!!$                EXIT
!!$             END IF
!!$          END DO

          iSegzTT = MAXLOC( SSP%Seg%z( : ), SSP%Seg%z( : ) < x( 3 ) )

          IF ( iSegzTT( 1 ) > 0 .AND. iSegzTT( 1 ) < SSP%Nz ) THEN  ! iSegx MUST LIE IN [ 1, SSP%Nz - 1 ]
             iSegz = iSegzTT( 1 )   
          END IF

       END IF

       ! cz at the corners of the current rectangle
       cz11 = SSP%czMat3( iSegx,     iSegy    , iSegz )
       cz12 = SSP%czMat3( iSegx + 1, iSegy    , iSegz )
       cz21 = SSP%czMat3( iSegx,     iSegy + 1, iSegz )
       cz22 = SSP%czMat3( iSegx + 1, iSegy + 1, iSegz )

       ! for this depth, x( 3 ) get the sound speed at the corners of the current rectangle
       s3 = x( 3 ) - SSP%Seg%z( iSegz )
       c11 = SSP%cMat3( iSegx,     iSegy    , iSegz ) + s3 * cz11
       c21 = SSP%cMat3( iSegx + 1, iSegy    , iSegz ) + s3 * cz12
       c12 = SSP%cMat3( iSegx,     iSegy + 1, iSegz ) + s3 * cz21
       c22 = SSP%cMat3( iSegx + 1, iSegy + 1, iSegz ) + s3 * cz22

       ! s1 = proportional distance of x( 1 ) in x
       s1 = ( x( 1 ) - SSP%Seg%x( iSegx ) ) / ( SSP%Seg%x( iSegx + 1 ) - SSP%Seg%x( iSegx ) )
       s1 = MIN( s1, 1.0D0 )   ! force piecewise constant extrapolation for points outside the box
       s1 = MAX( s1, 0.0D0 )   ! "

       ! s2 = proportional distance of x( 2 ) in y
       s2 = ( x( 2 ) - SSP%Seg%y( iSegy ) ) / ( SSP%Seg%y( iSegy + 1 ) - SSP%Seg%y( iSegy ) )
       s2 = MIN( s2, 1.0D0 )   ! force piecewise constant extrapolation for points outside the box
       s2 = MAX( s2, 0.0D0 )   ! "

       ! interpolate the soundspeed in the x direction, at the two endpoints in y (top and bottom sides of rectangle)
       c1 = c11 + s1 * ( c21 - c11 )
       c2 = c12 + s1 * ( c22 - c12 )
       !c  = ( 1.0D0 - s2 ) * c1  + s2 * c2   ! interpolation in y
       cy = ( c2 - c1 ) / ( SSP%Seg%y( iSegy + 1 ) - SSP%Seg%y( iSegy ) )

       ! interpolate the soundspeed in the y direction, at the two endpoints in x (left and right sides of rectangle)
       c1 = c11 + s2 * ( c12 - c11 )
       c2 = c21 + s2 * ( c22 - c21 )

       c  = c1  + s1 * ( c2 - c1   )   ! interpolation in x

       ! interpolate the attenuation !!!! This will use the wrong segment if the ssp in the envfil is sampled at different depths
       s3 = s3 / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )   ! convert s3 to a proportional distance in the layer
       cimag = AIMAG( ( 1.0D0 - s3 ) * SSP%c( Isegz )  + s3 * SSP%c( Isegz + 1 ) )   ! volume attenuation is taken from the single c(z) profile

       cx = ( c2 - c1 ) / ( SSP%Seg%x( iSegx + 1 ) - SSP%Seg%x( iSegx ) )

       ! same thing on cz
       cz1 = cz11 + s2 * ( cz21 - cz11 )
       cz2 = cz12 + s2 * ( cz22 - cz12 )
       cz  = cz1  + s1 * ( cz2  - cz1  )   ! interpolation in z

       !gradc = [ cx, cy, cz ]
       gradc( 1 ) = cx
       gradc( 2 ) = cy
       gradc( 3 ) = cz

       cxx   = 0.0D0
       cyy   = 0.0D0
       czz   = 0.0D0
       cxy   = 0.0D0
       cxz   = 0.0D0
       cyz   = 0.0D0
       ! write( *, FMT="( 'SSP', I3, 2X, 4F10.2, F10.4 )" ) layer, x, c, cz

       ! linear interpolation for density
       W   = ( x( 3 ) - SSP%z( iSegz ) ) / ( SSP%z( iSegz + 1 ) - SSP%z( iSegz ) )
       rho = ( 1.0D0 - W ) * SSP%rho( iSegz ) + W * SSP%rho( iSegz + 1 )
    END IF

  END SUBROUTINE Hexahedral

!**********************************************************************!

  SUBROUTINE Analytic( x, c, cimag, gradc, crr, crz, czz, rho )

    REAL (KIND=8), INTENT( IN  ) :: x( 2 )
    REAL (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
    REAL (KIND=8)                :: c0, cr, cz, DxtDz, xt

    iSegz = 1
    c0    = 1500.0
    rho   = 1.0

    ! homogeneous halfspace was removed since BELLHOP needs to get gradc just a little below the boundaries, on ray reflection

!!$  IF ( x( 2 ) < 5000.0 ) THEN
    xt    = 2.0 * ( x( 2 ) - 1300.0 ) / 1300.0
    DxtDz = 2.0 / 1300.0
    c     = C0 * ( 1.0 + 0.00737*( xt - 1.0 + EXP( -xt ) ) )
    cimag = 0.
    cz    = C0 * 0.00737 * ( 1.0 - EXP( -xt ) ) * DxtDz
    czz   = C0 * 0.00737 * EXP( -xt ) * DxtDz ** 2
!!$  ELSE
!!$     ! Homogeneous half-space
!!$     xt   = 2.0 * ( 5000.0 - 1300.0 ) / 1300.0
!!$     c    = C0 * ( 1.0 + 0.00737 * ( xt - 1.0 + EXP( -xt ) ) )
!!$     cz   = 0.0
!!$     czz  = 0.0
!!$  END IF

    cr = 0.0
    gradc = [ cr, cz ]
    crz = 0.0
    crr = 0.0

    RETURN
  END SUBROUTINE Analytic

  !**********************************************************************!

  SUBROUTINE AnalyticCosh( x, c, cimag, gradc, crr, crz, czz, rho )

    REAL (KIND=8), INTENT( IN  ) :: x( 2 )
    REAL (KIND=8), INTENT( OUT ) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
    REAL (KIND=8)                :: cr, cz, A, B, D, z

    iSegz = 1
    rho   = 1.0
    
    D = 3000.0
    z = x( 2 )

    A = 1500.0
    B = 0.0003
    D = 1500.0

    c   =         A * COSH( B * ( Z - D ) )
    cz  =     B * A * SINH( B * ( Z - D ) )
    czz = B * B * A * COSH( B * ( Z - D ) )
      
    cimag = 0.

    cr = 0.0
    gradc = [ cr, cz ]
    crz = 0.0
    crr = 0.0

    RETURN
  END SUBROUTINE AnalyticCosh

!**********************************************************************!

SUBROUTINE Analytic3D( x, c, cimag, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho )

  REAL (KIND=8) :: x( 3 ), c, cimag, gradc( 3 ), cxx, cyy, czz, cxy, cxz, cyz, c0, W, Wz, epsilon, epsilon_y
  REAL (KIND=8) :: rho
  
  iSegz = 1
  c0    = 1500.0
  rho   = 1.0

!!$  IF ( x( 3 ) .LT. 5000.0 ) THEN
     epsilon   = 0.00737 + x( 2 ) / 100000.0 * 0.003
     epsilon_y = 0.003 / 100000.0

     W   = 2.0 * ( x( 3 ) - 1300.0 ) / 1300.0
     Wz  = 2.0 / 1300.0

     c           = c0 * ( 1.0 + epsilon * ( W - 1.0 + EXP( -W ) ) )
     cimag       = 0.
     gradc( 2 )  = c0 * epsilon_y * ( W - 1.0 + EXP( -W ) )
     gradc( 3 )  = c0 * epsilon * ( 1.0 - EXP( -W ) ) * Wz
     czz         = c0 * epsilon * EXP( -W ) * Wz **2
     cyz         = c0 * epsilon_y * ( 1.0 - EXP( -W ) ) * Wz
!!$  ELSE    ! HOMOGENEOUS HALF-SPACE
!!$     W           = 2.0 * ( 5000.0 - 1300.0 ) / 1300.0
!!$     c           = c0*( 1.0 + 0.00737 * ( W - 1.0 + EXP( -W ) ) )
!!$     gradc( 2 )  = 0.0
!!$     gradc( 3 )  = 0.0
!!$     czz         = 0.0
!!$     cyz         = 0.0
!!$  END IF

  gradc( 1 )  = 0.0
  cxx = 0.0
  cyy = 0.0
  cxz = 0.0
  cxy = 0.0

  RETURN
END SUBROUTINE Analytic3D

!**********************************************************************!

  SUBROUTINE ReadSSP( Depth, freq )

    ! reads the SSP data from the environmental file and convert to Nepers/m

    USE AttenMod

    REAL     (KIND=8), INTENT(IN) :: freq, Depth

    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Sound speed profile:'
    WRITE( PRTFile, "( '      z         alphaR      betaR     rho        alphaI     betaI'    )" )
    WRITE( PRTFile, "( '     (m)         (m/s)      (m/s)   (g/cm^3)      (m/s)     (m/s)', / )" )

    SSP%NPts = 1

    DO iz = 1, MaxSSP

       READ(  ENVFile, *    ) SSP%z( iz ), alphaR, betaR, rhoR, alphaI, betaI
       WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) SSP%z( iz ), alphaR, betaR, rhoR, alphaI, betaI

       SSP%c(   iz ) = CRCI( SSP%z( iz ), alphaR, alphaI, freq, freq, SSP%AttenUnit, betaPowerLaw, fT )
       SSP%rho( iz ) = rhoR

       ! verify that the depths are monotone increasing
       IF ( iz > 1 ) THEN
          IF ( SSP%z( iz ) .LE. SSP%z( iz - 1 ) ) THEN
              WRITE( PRTFile, * ) 'Bad depth in SSP: ', SSP%z( iz )
              CALL ERROUT( 'ReadSSP', 'The depths in the SSP must be monotone increasing' )
          END IF
       END IF

       ! compute gradient, cz
       IF ( iz > 1 ) SSP%cz( iz - 1 )  = ( SSP%c( iz ) - SSP%c( iz - 1 ) ) / &
                                         ( SSP%z( iz ) - SSP%z( iz - 1 ) )

       ! Did we read the last point?
       IF ( ABS( SSP%z( iz ) - Depth ) < 100. * EPSILON( 1.0e0 ) ) THEN
          SSP%Nz = SSP%NPts
          IF ( SSP%NPts == 1 ) THEN
              WRITE( PRTFile, * ) '#SSP points: ', SSP%NPts
              CALL ERROUT( 'ReadSSP', 'The SSP must have at least 2 points' )
          END IF

          RETURN
       ENDIF

       SSP%NPts = SSP%NPts + 1
    END DO
 
    ! Fall through means too many points in the profile
    WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
    CALL ERROUT( 'ReadSSP', 'Number of SSP points exceeds limit' )

  END SUBROUTINE ReadSSP

END MODULE sspmod
 
