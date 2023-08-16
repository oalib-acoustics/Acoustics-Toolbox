MODULE bdrymod

  ! Loads altimetry (top bdry) and bathymetry (bottom bdry) data

  !USE norms
  USE monotonicMod
  USE FatalError

  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: ATIFile = 40, BTYFile = 41, Number_to_Echo = 21
  INTEGER            :: IsegTop, IsegBot                ! indices that point to the current active segment
  INTEGER, PROTECTED :: NATIPts = 2, NBTYPts = 2
  INTEGER            :: ii, IOStat, IAllocStat, iSmallStepCtr = 0

  REAL (KIND=8)      :: rTopseg( 2 ), rBotseg( 2 )      ! range intervals defining the current active segment
  CHARACTER  (LEN=2) :: atiType= 'LS', btyType = 'LS'

  ! Halfspace properties
  TYPE HSInfo2
     REAL     (KIND=8) :: alphaR, alphaI, betaR, betaI  ! compressional and shear wave speeds/attenuations in user units
     COMPLEX  (KIND=8) :: cP, cS                        ! P-wave, S-wave speeds
     REAL     (KIND=8) :: rho, Depth                    ! density, depth
     CHARACTER (LEN=1) :: BC                            ! Boundary condition type
     CHARACTER (LEN=6) :: Opt
  END TYPE

  TYPE BdryPt
     REAL    (KIND=8) :: x( 2 ), t( 2 ), n( 2 )         ! coordinate, tangent, and outward normal for a segment
     REAL    (KIND=8) :: Nodet( 2 ), Noden( 2 )         ! tangent and normal at the node, if the curvilinear option is used
     REAL    (KIND=8) :: Len, Kappa                     ! length and curvature of a segment
     REAL    (KIND=8) :: Dx, Dxx, Dss                   ! first, second derivatives wrt depth; s is along tangent
     TYPE( HSInfo2 )  :: HS
  END TYPE

  TYPE(BdryPt), ALLOCATABLE :: Top( : ), Bot( : )

CONTAINS

  SUBROUTINE ReadATI( FileRoot, TopATI, DepthT, PRTFile )

    ! Reads in the top altimetry

    INTEGER,            INTENT( IN ) :: PRTFile
    CHARACTER (LEN= 1), INTENT( IN ) :: TopATI
    REAL      (KIND=8), INTENT( IN ) :: DepthT
    CHARACTER (LEN=80), INTENT( IN ) :: FileRoot

    SELECT CASE ( TopATI )
    CASE ( '~', '*' )
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Using top-altimetry file'

       OPEN( UNIT = ATIFile,   FILE = TRIM( FileRoot ) // '.ati', STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
       IF ( IOsTAT /= 0 ) THEN
          WRITE( PRTFile, * ) 'ATIFile = ', TRIM( FileRoot ) // '.ati'
          CALL ERROUT( 'ReadATI', 'Unable to open altimetry file' )
       END IF

       READ(  ATIFile, * ) atiType
       AltiType: SELECT CASE ( atiType( 1 : 1 ) )
       CASE ( 'C' )
          WRITE( PRTFile, * ) 'Curvilinear Interpolation'
       CASE ( 'L' )
          WRITE( PRTFile, * ) 'Piecewise linear interpolation'
       CASE DEFAULT
          CALL ERROUT( 'ReadATI', 'Unknown option for selecting altimetry interpolation' )
       END SELECT AltiType

       READ(  ATIFile, * ) NatiPts
       WRITE( PRTFile, * ) 'Number of altimetry points = ', NatiPts
       NatiPts = NatiPts + 2   ! we'll be extending the altimetry to infinity to the left and right

       ALLOCATE( Top( NatiPts ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
            CALL ERROUT( 'BELLHOP:ReadATI', 'Insufficient memory for altimetry data: reduce # ati points' )

       WRITE( PRTFile, * )
       AltiTypeB: SELECT CASE ( atiType( 2 : 2 ) )
       CASE ( 'S', '' )
          WRITE( PRTFile, * ) 'Short format (altimetry only)'
          WRITE( PRTFile, * ) ' Range (km)  Depth (m)'
       CASE ( 'L' )
          WRITE( PRTFile, * ) 'Long format (altimetry and geoacoustics)'
          WRITE( PRTFile, "( ' Range (km)  Depth (m)  alphaR (m/s)  betaR  rho (g/cm^3)  alphaI     betaI', / )" )
       CASE DEFAULT
          CALL ERROUT( 'ReadBTY', 'Unknown option for selecting altimetry interpolation' )
       END SELECT AltiTypeB

       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) ' Range (km)  Depth (m)'

       atiPt: DO ii = 2, NatiPts - 1

          SELECT CASE ( atiType( 2 : 2 ) )
          CASE ( 'S', '' )
             READ(  ATIFile, * ) Top( ii )%x
             IF ( ii < Number_to_Echo .OR. ii == NatiPts - 1 ) THEN   ! echo some values
                WRITE( PRTFile, FMT = "(2G11.3)" ) Top( ii )%x
             END IF
          CASE ( 'L' )
             READ(  ATIFile, * )                   Top( ii )%x, Top( ii )%HS%alphaR, Top( ii )%HS%betaR, Top( ii )%HS%rho, &
                                                                Top( ii )%HS%alphaI, Top( ii )%HS%betaI
             IF ( ii < Number_to_Echo .OR. ii == NatiPts - 1 ) THEN   ! echo some values
                WRITE( PRTFile, FMT = "(7G11.3)" ) Top( ii )%x, Top( ii )%HS%alphaR, Top( ii )%HS%betaR, Top( ii )%HS%rho, &
                                                                Top( ii )%HS%alphaI, Top( ii )%HS%betaI
             END IF
          CASE DEFAULT
             CALL ERROUT( 'ReadATI', 'Unknown option for selecting altimetry option' )
          END SELECT

          IF ( Top( ii )%x( 2 ) < DepthT ) THEN
             CALL ERROUT( 'BELLHOP:ReadATI', 'Altimetry rises above highest point in the sound speed profile' )
          END IF
       END DO atiPt

       CLOSE( ATIFile )

       Top( : )%x( 1 ) = 1000.0 * Top( : )%x( 1 )   ! Convert ranges in km to m

    CASE DEFAULT   ! no altimetry given, use SSP depth for flat top
       ALLOCATE( Top( 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP', 'Insufficient memory for altimetry data'  )
       Top( 1 )%x = [ -sqrt( huge( Top( 1 )%x( 1 ) ) ) / 1.0d5, DepthT ]
       Top( 2 )%x = [  sqrt( huge( Top( 1 )%x( 1 ) ) ) / 1.0d5, DepthT ]
    END SELECT

    CALL ComputeBdryTangentNormal( Top, 'Top' )

    IF ( .NOT. monotonic( Top%x( 1 ), NAtiPts ) ) THEN
       CALL ERROUT( 'BELLHOP:ReadATI', 'Altimetry ranges are not monotonically increasing' )
    END IF 
 
  END SUBROUTINE ReadATI

  ! **********************************************************************!

  SUBROUTINE ReadBTY( FileRoot, BotBTY, DepthB, PRTFile )

    ! Reads in the bottom bathymetry

    INTEGER,            INTENT( IN ) :: PRTFile
    CHARACTER (LEN= 1), INTENT( IN ) :: BotBTY
    REAL      (KIND=8), INTENT( IN ) :: DepthB
    CHARACTER (LEN=80), INTENT( IN ) :: FileRoot

    SELECT CASE ( BotBTY )
    CASE ( '~', '*' )
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Using bottom-bathymetry file'

       OPEN( UNIT = BTYFile,   FILE = TRIM( FileRoot ) // '.bty', STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
       IF ( IOsTAT /= 0 ) THEN
          WRITE( PRTFile, * ) 'BTYFile = ', TRIM( FileRoot ) // '.bty'
          CALL ERROUT( 'ReadBTY', 'Unable to open bathymetry file' )
       END IF
 
       READ( BTYFile, * ) btyType
 
       BathyType: SELECT CASE ( btyType( 1 : 1 ) )
       CASE ( 'C' )
          WRITE( PRTFile, * ) 'Curvilinear Interpolation'
       CASE ( 'L' )
          WRITE( PRTFile, * ) 'Piecewise linear interpolation'
       CASE DEFAULT
          CALL ERROUT( 'ReadBTY', 'Unknown option for selecting bathymetry interpolation' )
       END SELECT BathyType

       READ(  BTYFile, * ) NbtyPts
       WRITE( PRTFile, * ) 'Number of bathymetry points = ', NbtyPts

       NbtyPts = NbtyPts + 2   ! we'll be extending the bathymetry to infinity on both sides
       ALLOCATE( Bot( NbtyPts ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
            CALL ERROUT( 'BELLHOP:ReadBTY', 'Insufficient memory for bathymetry data: reduce # bty points' )

       WRITE( PRTFile, * )
       BathyTypeB: SELECT CASE ( btyType( 2 : 2 ) )
       CASE ( 'S', '' )
          WRITE( PRTFile, * ) 'Short format (bathymetry only)'
          WRITE( PRTFile, * ) ' Range (km)  Depth (m)'
       CASE ( 'L' )
          WRITE( PRTFile, * ) 'Long format (bathymetry and geoacoustics)'
          WRITE( PRTFile, "( ' Range (km)  Depth (m)  alphaR (m/s)  betaR  rho (g/cm^3)  alphaI     betaI', / )" )
       CASE DEFAULT
          CALL ERROUT( 'ReadBTY', 'Unknown option for selecting bathymetry interpolation' )
       END SELECT BathyTypeB

       btyPt: DO ii = 2, NbtyPts - 1

          SELECT CASE ( btyType( 2 : 2 ) )
          CASE ( 'S', '' )   ! short format
             READ(  BTYFile, * ) Bot( ii )%x
             IF ( ii < Number_to_Echo .OR. ii == NbtyPts - 1 ) THEN   ! echo some values
                WRITE( PRTFile, FMT = "(2G11.3)" ) Bot( ii )%x
             END IF
          CASE ( 'L' )       ! long format
             READ(  BTYFile, * )                   Bot( ii )%x, Bot( ii )%HS%alphaR, Bot( ii )%HS%betaR, Bot( ii )%HS%rho, &
                                                                Bot( ii )%HS%alphaI, Bot( ii )%HS%betaI
             IF ( ii < Number_to_Echo .OR. ii == NbtyPts - 1 ) THEN   ! echo some values
                WRITE( PRTFile, FMT="( F10.2, F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
                                                   Bot( ii )%x, Bot( ii )%HS%alphaR, Bot( ii )%HS%betaR, Bot( ii )%HS%rho, &
                                                                Bot( ii )%HS%alphaI, Bot( ii )%HS%betaI
             END IF
          CASE DEFAULT
             CALL ERROUT( 'ReadBTY', 'Unknown option for selecting bathymetry option' )
          END SELECT

          IF ( Bot( ii )%x( 2 ) > DepthB ) THEN
             CALL ERROUT( 'BELLHOP:ReadBTY', 'Bathymetry drops below lowest point in the sound speed profile' )
          END IF
 
       END DO btypt

       CLOSE( BTYFile )

       Bot( : )%x( 1 ) = 1000.0 * Bot( : )%x( 1 )   ! Convert ranges in km to m

    CASE DEFAULT   ! no bathymetry given, use SSP depth for flat bottom
       ALLOCATE( Bot( 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP', 'Insufficient memory for bathymetry data'  )
       Bot( 1 )%x = [ -sqrt( huge( Bot( 1 )%x( 1 ) ) ) / 1.0d5, DepthB ]
       Bot( 2 )%x = [  sqrt( huge( Bot( 1 )%x( 1 ) ) ) / 1.0d5, DepthB ]
    END SELECT

    CALL ComputeBdryTangentNormal( Bot, 'Bot' )

    IF ( .NOT. monotonic( Bot%x( 1 ), NBtyPts ) ) THEN
       CALL ERROUT( 'BELLHOP:ReadBTY', 'Bathymetry ranges are not monotonically increasing' )
    END IF 

  END SUBROUTINE ReadBTY

  ! **********************************************************************!

  SUBROUTINE ComputeBdryTangentNormal( Bdry, BotTop )

    ! Does some pre-processing on the boundary points to pre-compute segment
    ! lengths  (%Len),
    ! tangents (%t, %nodet),
    ! normals  (%n, %noden), and
    ! curvatures (%kappa)
    !
    ! The boundary is also extended with a constant depth to infinity to cover cases where the ray
    ! exits the domain defined by the user

    INTEGER                          :: NPts = 0
    REAL      (KIND=8), ALLOCATABLE  :: phi( : )
    REAL      (KIND=8)               :: sss
    TYPE(BdryPt)                     :: Bdry( : )
    CHARACTER (LEN=3),  INTENT( IN ) :: BotTop           ! Flag indicating bottom or top reflection
    CHARACTER (LEN=2)                :: CurvilinearFlag = '-'

    SELECT CASE ( BotTop )
    CASE ( 'Bot' )
       NPts = NbtyPts
       CurvilinearFlag = btyType
    CASE ( 'Top' )
       NPts = NatiPts
       CurvilinearFlag = atiType
    END SELECT

    ! extend the bathymetry to +/- infinity in a piecewise constant fashion

    Bdry( 1    )%x( 1 ) = -sqrt( huge( Bdry( 1 )%x( 1 ) ) ) / 1.0d5
    Bdry( 1    )%x( 2 ) = Bdry( 2        )%x( 2 )
    Bdry( 1    )%HS     = Bdry( 2        )%HS
    Bdry( NPts )%x( 1 ) = +sqrt( huge( Bdry( 1 )%x( 1 ) ) ) / 1.0d5
    Bdry( NPts )%x( 2 ) = Bdry( NPts - 1 )%x( 2 )
    Bdry( NPts )%HS     = Bdry( NPts - 1 )%HS

    ! compute tangent and outward-pointing normal to each bottom segment
    ! tBdry( 1, : ) = xBdry( 1, 2:NPts ) - xBdry( 1, 1:NPts - 1 )
    ! tBdry( 2, : ) = xBdry( 2, 2:NPts ) - xBdry( 2, 1:NPts - 1 )
    ! above caused compiler problems

    BoundaryPt: DO ii = 1, NPts - 1
       Bdry( ii )%t   = Bdry( ii + 1 )%x      - Bdry( ii )%x
       Bdry( ii )%Dx  = Bdry( ii )%t( 2 ) / Bdry( ii )%t( 1 )   ! first derivative
       ! write( *, * ) 'Dx, t', Bdry( ii )%Dx, Bdry( ii )%x, 1 / ( Bdry( ii )%x( 2 ) / 500 )

       ! normalize the tangent vector
       Bdry( ii )%Len = NORM2( Bdry( ii )%t )
       Bdry( ii )%t   = Bdry( ii )%t / Bdry( ii )%Len

       SELECT CASE ( BotTop )
       CASE ( 'Bot' )
          Bdry( ii )%n( 1 ) = -Bdry( ii )%t( 2 )
          Bdry( ii )%n( 2 ) = +Bdry( ii )%t( 1 )
       CASE ( 'Top' )
          Bdry( ii )%n( 1 ) = +Bdry( ii )%t( 2 )
          Bdry( ii )%n( 2 ) = -Bdry( ii )%t( 1 )
       END SELECT

    END DO BoundaryPt

    IF ( CurvilinearFlag( 1 : 1 ) == 'C' ) THEN ! curvilinear option: compute tangent and normal at node by averaging normals on adjacent segments
       ! averaging two centered differences is equivalent to forming a single centered difference of two steps ...
       DO ii = 2, NPts - 1
          sss = Bdry( ii - 1 )%Len / ( Bdry( ii - 1 )%Len + Bdry( ii )%Len )
          sss = 0.5
          Bdry( ii )%Nodet = ( 1.0 - sss ) * Bdry( ii - 1 )%t + sss * Bdry( ii )%t
       END DO

       Bdry( 1    )%Nodet = [ 1.0, 0.0 ]   ! tangent left-end  node
       Bdry( NPts )%Nodet = [ 1.0, 0.0 ]   ! tangent right-end node

       SELECT CASE ( BotTop )
       CASE ( 'Bot' )
          Bdry( : )%Noden( 1 ) = -Bdry( : )%Nodet( 2 )
          Bdry( : )%Noden( 2 ) = +Bdry( : )%Nodet( 1 )
       CASE ( 'Top' )
          Bdry( : )%Noden( 1 ) = +Bdry( : )%Nodet( 2 )
          Bdry( : )%Noden( 2 ) = -Bdry( : )%Nodet( 1 )
       END SELECT

       ! compute curvature in each segment
       ALLOCATE( phi( NPts ), Stat = IAllocStat )
       phi = atan2( Bdry( : )%Nodet( 2 ), Bdry( : )%Nodet( 1 ) )   ! this is the angle at each node

       DO ii = 1, NPts - 1
          Bdry( ii )%kappa = ( phi( ii + 1 ) - phi( ii ) ) / Bdry( ii )%Len ! this is curvature = dphi/ds
          Bdry( ii )%Dxx   = ( Bdry( ii + 1 )%Dx     - Bdry( ii )%Dx     ) / &   ! second derivative
                             ( Bdry( ii + 1 )%x( 1 ) - Bdry( ii )%x( 1 ) )
          Bdry( ii )%Dss   = Bdry( ii )%Dxx * Bdry( ii )%t( 1 ) ** 3   ! derivative in direction of tangent
          !write( *, * ) 'kappa, Dss, Dxx', Bdry( ii )%kappa, Bdry( ii )%Dss, Bdry( ii )%Dxx, &
               ! 1 / ( ( 8 / 1000 ** 2 ) * ABS( Bdry( ii )%x( 2 ) ) ** 3 ), Bdry( ii )%x( 2 )
               !  -1 / ( 4 * ( Bdry( ii )%x( 2 ) ) ** 3 / 1000000 ), Bdry( ii )%x( 2 )

          Bdry( ii )%kappa = Bdry( ii )%Dss   !over-ride kappa !!!!!
       END DO

       DEALLOCATE( phi )   ! should happen automatically when phi goes out of scope ...
    ELSE
       Bdry%kappa = 0
    END IF

  END SUBROUTINE ComputeBdryTangentNormal

  ! **********************************************************************!

  SUBROUTINE GetTopSeg( r )

    ! Get the Top segment info (index and range interval) for range, r

    INTEGER, PARAMETER :: PRTFile = 6
    INTEGER IsegTopT( 1 )
    REAL (KIND=8), INTENT( IN ) :: r

    IsegTopT = MAXLOC( Top( : )%x( 1 ), Top( : )%x( 1 ) < r )

    IF ( IsegTopT( 1 ) > 0 .AND. IsegTopT( 1 ) < NatiPts ) THEN  ! IsegTop MUST LIE IN [ 1, NatiPts-1 ]
       IsegTop = IsegTopT( 1 )
       rTopSeg = [ Top( IsegTop )%x( 1 ), Top( IsegTop + 1 )%x( 1 ) ]   ! segment limits in range
    ELSE
       WRITE( PRTFile, * ) 'r = ', r
       WRITE( PRTFile, * ) 'rLeft  = ', Top( 1       )%x( 1 )
       WRITE( PRTFile, * ) 'rRight = ', Top( NatiPts )%x( 1 )
       CALL ERROUT( 'GetTopSeg', 'Top altimetry undefined above the ray' )
    ENDIF

  END SUBROUTINE GetTopSeg

  ! **********************************************************************!

  SUBROUTINE GetBotSeg( r )

    ! Get the Bottom segment info (index and range interval) for range, r

    INTEGER, PARAMETER :: PRTFile = 6
    INTEGER IsegBotT( 1 )
    REAL (KIND=8), INTENT( IN ) :: r

    IsegBotT = MAXLOC( Bot( : )%x( 1 ), Bot( : )%x( 1 ) < r )

    IF ( IsegBotT( 1 ) > 0 .AND. IsegBotT( 1 ) < NbtyPts ) THEN  ! IsegBot MUST LIE IN [ 1, NbtyPts-1 ]
       IsegBot = IsegBotT( 1 )   
       rBotSeg = [ Bot( IsegBot )%x( 1 ), Bot( IsegBot + 1 )%x( 1 ) ]   ! segment limits in range
    ELSE
       WRITE( PRTFile, * ) 'r = ', r
       WRITE( PRTFile, * ) 'rLeft  = ', Bot( 1       )%x( 1 )
       WRITE( PRTFile, * ) 'rRight = ', Bot( NbtyPts )%x( 1 )
       CALL ERROUT( 'GetBotSeg', 'Bottom bathymetry undefined below the source' )
    ENDIF

  END SUBROUTINE GetBotSeg

END MODULE bdrymod
