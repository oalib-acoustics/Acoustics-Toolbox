MODULE bdry3Dmod

  ! Loads
  ! altimetry (top bdry) and bathymetry (bottom bdry) data
  ! This version is for BELLHOP3D
  !
  ! x = coordinate of boundary
  ! t = tangent for a facet
  ! n = normal  for a facet (outward pointing)
  ! n1, n2 are normals for each of the triangles in a pair, n is selected from those
  ! Len = length of tangent (temporary variable to normalize tangent)

  USE monotonicMod
  USE SubTabulate
  USE FatalError

  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: ATIFile = 40, BTYFile = 41, Number_to_Echo = 21
  INTEGER            :: IsegTopT( 1 ), IsegBotT( 1 ), IsegTopx, IsegTopy, IsegBotx, IsegBoty, &
       NATIPts( 2 ), NBTYPts( 2 )
  INTEGER            :: ix, iy, IOStat, IAllocStat, iSmallStepCtr = 0
  REAL (KIND=8) :: xTopseg( 2 ), yTopseg( 2 ), xBotseg( 2 ), yBotseg( 2 ), &
       Topx( 3 ), Botx( 3 ), &                            ! coordinates of corner of active rectangle
       Topn( 3 ), Botn( 3 )                               ! tangent and normal    of active triangle
  REAL (KIND=8), PROTECTED :: Top_tri_n( 2 ), Bot_tri_n( 2 )   ! triangle normals
  REAL (KIND=8), PROTECTED :: Top_deltax, Top_deltay, Bot_deltax, Bot_deltay   ! lengths of sides for active rectangle
  REAL (KIND=8), PARAMETER :: big = 1E25                  ! large number used for domain termination when no altimetry/bathymetry given

  CHARACTER  (LEN=1) :: atiType, btyType

  REAL (KIND=8), ALLOCATABLE :: BotGlobalx( : ), BotGlobaly( : ), TopGlobalx( : ), TopGlobaly( : )
  TYPE BdryPt
     REAL (KIND=8) :: x( 3 ), t( 3 ), n( 3 ), n1( 3 ), n2( 3 ), Len, Noden( 3 ), Noden_unscaled( 3 ), z_xx, z_xy, z_yy, &
          phi_xx, phi_xy, phi_yy, kappa_xx, kappa_xy, kappa_yy
  END TYPE BdryPt

  TYPE(BdryPt), ALLOCATABLE :: Bot( :, : ), Top( :, : )

CONTAINS

  SUBROUTINE ReadATI3D( FileRoot, TopATI, DepthT, PRTFile )

    CHARACTER (LEN= 1), INTENT( IN ) :: TopATI        ! Set to '~' if altimetry is not flat
    INTEGER,            INTENT( IN ) :: PRTFile       ! unit number for print file
    REAL      (KIND=8), INTENT( IN ) :: DepthT        ! Nominal top depth
    CHARACTER (LEN=80), INTENT( IN ) :: FileRoot
    REAL (KIND=8), ALLOCATABLE :: Temp( : )

    SELECT CASE ( TopATI )
    CASE ( '~', '*' )
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Using top-altimetry file'

       OPEN( UNIT = ATIFile, FILE = TRIM( FileRoot ) // '.ati', STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
       IF ( IOsTAT /= 0 ) THEN
          WRITE( PRTFile, * ) 'ATIFile = ', TRIM( FileRoot ) // '.ati'
          CALL ERROUT( 'ReadATI', 'Unable to open altimetry file' )
       END IF

       READ(  ATIFile, * ) atiType
       SELECT CASE ( atiType )
       CASE ( 'R' )
          WRITE( PRTFile, * ) 'Regular grid for a 3D run'
       CASE ( 'C' )
          WRITE( PRTFile, * ) 'Regular grid for a 3D run (curvilinear)'
       CASE DEFAULT
          CALL ERROUT( 'ReadATI3D', 'Unknown option for selecting altimetry interpolation' )
       END SELECT

       ! x values
       READ(  ATIFile, * ) NatiPts( 1 )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of altimetry points in x', NatiPts( 1 )

       ALLOCATE( TopGlobalx( MAX( NatiPts( 1 ), 3 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Insufficient memory for altimetry data: reduce # ati points' )

       TopGlobalx( 3 ) = -999.9
       READ(  ATIFile, * ) TopGlobalx( 1 : NatiPts( 1 ) )
       CALL SubTab( TopGlobalx, NatiPts( 1 ) )
       WRITE( PRTFile, "( 5G14.6 )" ) ( TopGlobalx( ix ), ix = 1, MIN( NatiPts( 1 ), Number_to_Echo ) )
       IF ( NatiPts( 1 ) > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )" ) ' ... ', TopGlobalx( NatiPts( 1 ) )

       ! y values
       READ(  ATIFile, * ) NatiPts( 2 )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of altimetry points in y', NatiPts( 2 )

       ALLOCATE( TopGlobaly( MAX( NatiPts( 2 ), 3 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Insufficient memory for altimetry data: reduce # ati points' )

       TopGlobaly( 3 ) = -999.9
       READ(  ATIFile, * ) TopGlobaly( 1 : NatiPts( 2 ) )
       CALL SubTab( TopGlobaly, NatiPts( 2 ) )
       WRITE( PRTFile, "( 5G14.6 )" ) ( TopGlobaly( iy ), iy = 1, MIN( NatiPts( 2 ), Number_to_Echo ) )
       IF ( NatiPts( 2 ) > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )"  ) ' ... ', TopGlobaly( NatiPts( 2 ) )

       TopGlobalx = 1000. * TopGlobalx   ! convert km to m
       TopGlobaly = 1000. * TopGlobaly

       ! z values
       ALLOCATE( Top( NatiPts( 1 ), NatiPts( 2 ) ), Temp( NatiPts( 1 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Insufficient memory for altimetry data: reduce # ati points' )

       WRITE( PRTFile, * )

       DO iy = 1, NatiPts( 2 )
          READ( ATIFile, * ) Top( :, iy )%x( 3 )   ! read a row of depths

          ! IF ( iy < Number_to_Echo .OR. iy == NatiPts( 2 ) ) THEN   ! echo some values
          !    WRITE( PRTFile, FMT = "(G11.3)" ) Top( :, iy )%x( 3 )
          ! END IF
          ! IF ( ANY( Top( :, iy )%x( 3 ) > DepthB ) ) THEN
          !    CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Altimetry drops below lowest point in the sound speed profile' )
          ! END IF
       END DO

       CLOSE( ATIFile )

       IF ( ANY( ISNAN( Top( :, : )%x( 3 ) ) ) ) THEN
          WRITE( PRTFile, * ) 'Warning in BELLHOP3D - ReadATI3D : The altimetry file contains a NaN'
       END IF
 
       DO ix = 1, NatiPts( 1 ) 
          DO iy = 1, NatiPts( 2 )
             Top( ix, iy )%x( 1 ) = TopGlobalx( ix )
             Top( ix, iy )%x( 2 ) = TopGlobaly( iy )
          END DO
       END DO

       CALL ComputeBdryTangentNormal( Top, 'Top' )

       IF ( .NOT. monotonic( TopGlobalx, NAtiPts( 1 ) ) ) THEN
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Altimetry x-coordinates are not monotonically increasing' )
       END IF

       IF ( .NOT. monotonic( TopGlobaly, NAtiPts( 2 ) ) ) THEN
          CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Altimetry y-coordinates are not monotonically increasing' )
       END IF

    CASE DEFAULT   ! no altimetry given, use SSP depth for flat top
       atiType = 'R'
       NatiPts = [ 2, 2 ]
       ALLOCATE( TopGlobalx( 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Insufficient memory' )
       ALLOCATE( TopGlobaly( 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Insufficient memory' )
       ALLOCATE( Top( 2, 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP3D:ReadATI3D', 'Insufficient memory'  )

       !big = sqrt( huge( Top( 1, 1 )%x ) ) / 1.0d5

       TopGlobalx( 1 ) = -big
       TopGlobalx( 2 ) = +big

       TopGlobaly( 1 ) = -big
       TopGlobaly( 2 ) = +big

       Top_deltax = 2.0 * big
       Top_deltay = 2.0 * big

       Top( 1, 1 )%x = [ -big, -big, DepthT ]
       Top( 1, 2 )%x = [ -big,  big, DepthT ]
       Top( 2, 1 )%x = [  big, -big, DepthT ]
       Top( 2, 2 )%x = [  big,  big, DepthT ]

       Top( 1, 1 )%t  = [ 1.0, 0.0,  0.0 ]   ! tangent to top
       Top( 1, 1 )%n1 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 1, 1 )%n2 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 1, 2 )%t  = [ 1.0, 0.0,  0.0 ]   ! tangent to top
       Top( 1, 1 )%n1 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 1, 1 )%n2 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 2, 1 )%t  = [ 1.0, 0.0,  0.0 ]   ! tangent to top
       Top( 2, 1 )%n1 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 2, 1 )%n2 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 2, 2 )%t  = [ 1.0, 0.0,  0.0 ]   ! tangent to top
       Top( 2, 2 )%n1 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
       Top( 2, 2 )%n2 = [ 0.0, 0.0, -1.0 ]   ! outward-pointing normal
    END SELECT

    ! dummy TopSeg info to force GetTopSeg to search for the active segment on first call
    xTopSeg = [ +big, -big ]
    yTopSeg = [ +big, -big ]

  END SUBROUTINE ReadATI3D

  ! **********************************************************************!

  SUBROUTINE ReadBTY3D( FileRoot, BotBTY, DepthB, PRTFile )

    ! Reads in the bottom bathymetry

    CHARACTER (LEN= 1), INTENT( IN ) :: BotBTY        ! Set to '~' if bathymetry is not flat
    INTEGER,            INTENT( IN ) :: PRTFile       ! unit number for print file
    REAL      (KIND=8), INTENT( IN ) :: DepthB        ! Nominal bottom depth
    CHARACTER (LEN=80), INTENT( IN ) :: FileRoot
    REAL (KIND=8), ALLOCATABLE :: Temp( : )
 
    SELECT CASE ( BotBTY )
    CASE ( '~', '*' )
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Using bottom-bathymetry file'

       OPEN( UNIT = BTYFile, FILE = TRIM( FileRoot ) // '.bty', STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
       IF ( IOStat /= 0 ) THEN
         WRITE( PRTFile, * ) 'BTYFile = ', TRIM( FileRoot ) // '.bty'
         CALL ERROUT( 'ReadBTY3D', 'Unable to open bathymetry file' )
       END IF
 
       READ( BTYFile, * ) btyType

       SELECT CASE ( btyType )
       CASE ( 'R' )
          WRITE( PRTFile, * ) 'Regular grid for a 3D run'
       CASE ( 'C' )
          WRITE( PRTFile, * ) 'Regular grid for a 3D run (curvilinear)'
       CASE DEFAULT
          CALL ERROUT( 'ReadBTY3D', 'Unknown option for selecting bathymetry interpolation' )
       END SELECT

       ! x values
       READ(  BTYFile, * ) NbtyPts( 1 )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of bathymetry points in x', NbtyPts( 1 )

       ALLOCATE( BotGlobalx( MAX( NbtyPts( 1 ), 3 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Insufficient memory for bathymetry data: reduce # bty points' )

       BotGlobalx( 3 ) = -999.9
       READ(  BTYFile, * ) BotGlobalx( 1 : NbtyPts( 1 ) )
       CALL SubTab( BotGlobalx, NbtyPts( 1 ) )
       WRITE( PRTFile, "( 5G14.6 )" ) ( BotGlobalx( ix ), ix = 1, MIN( NbtyPts( 1 ), Number_to_Echo ) )
       IF ( NbtyPts( 1 ) > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )" ) ' ... ', BotGlobalx( NbtyPts( 1 ) )

       ! y values
       READ(  BTYFile, * ) NbtyPts( 2 )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of bathymetry points in y', NbtyPts( 2 )

       ALLOCATE( BotGlobaly( MAX( NbtyPts( 2 ), 3 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Insufficient memory for bathymetry data: reduce # bty points' )

       BotGlobaly( 3 ) = -999.9
       READ(  BTYFile, * ) BotGlobaly( 1 : NbtyPts( 2 ) )
       CALL SubTab( BotGlobaly, NbtyPts( 2 ) )
       WRITE( PRTFile, "( 5G14.6 )" ) ( BotGlobaly( iy ), iy = 1, MIN( NbtyPts( 2 ), Number_to_Echo ) )
       IF ( NbtyPts( 2 ) > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )" ) ' ... ', BotGlobaly( NbtyPts( 2 ) )

       BotGlobalx = 1000. * BotGlobalx   ! convert km to m
       BotGlobaly = 1000. * BotGlobaly

       ! z values
       ALLOCATE( Bot( NbtyPts( 1 ), NbtyPts( 2 ) ), Temp( NbtyPts( 1 ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Insufficient memory for bathymetry data: reduce # bty points' )

       WRITE( PRTFile, * )

       DO iy = 1, NbtyPts( 2 )
          READ( BTYFile, * ) Bot( :, iy )%x( 3 )    ! read a row of depths
          ! IF ( iy < Number_to_Echo .OR. iy == NbtyPts( 2 ) ) THEN   ! echo some values
          !    WRITE( PRTFile, FMT = "(G11.3)" ) Bot( :, iy )%x( 3 )
          ! END IF
          ! IF ( ANY( Bot( :, iy )%x( 3 ) > DepthB ) ) THEN
          !    CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Bathymetry drops below lowest point in the sound speed profile' )
          ! END IF
       END DO

       CLOSE( BTYFile )

       IF ( ANY( ISNAN( Bot( :, : )%x( 3 ) ) ) ) THEN
          WRITE( PRTFile, * ) 'Warning in BELLHOP3D - ReadBTY3D : The bathymetry file contains a NaN'
       END IF
 
       DO ix = 1, NbtyPts( 1 ) 
          DO iy = 1, NbtyPts( 2 )
             Bot( ix, iy )%x( 1 ) = BotGlobalx( ix )
             Bot( ix, iy )%x( 2 ) = BotGlobaly( iy )
          END DO
       END DO

       CALL ComputeBdryTangentNormal( Bot, 'Bot' )
       
       IF ( .NOT. monotonic( BotGlobalx, NBtyPts( 1 ) ) ) THEN
          CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Bathymetry x-coordinates are not monotonically increasing' )
       END IF

       IF ( .NOT. monotonic( BotGlobaly, NBTYPts( 2 ) ) ) THEN
          CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Bathymetry y-coordinates are not monotonically increasing' )
       END IF

    CASE DEFAULT   ! no bathymetry given, use SSP depth for flat bottom
       btyType = 'R'
       NbtyPts = [ 2, 2 ]
       ALLOCATE( BotGlobalx( 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Insufficient memory' )
       ALLOCATE( BotGlobaly( 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Insufficient memory' )
       ALLOCATE( Bot( 2, 2 ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP3D:ReadBTY3D', 'Insufficient memory'  )

       ! big = sqrt( huge( Bot( 1, 1 )%x ) ) / 1.0d5

       BotGlobalx( 1 ) = -big
       BotGlobalx( 2 ) = +big

       BotGlobaly( 1 ) = -big
       BotGlobaly( 2 ) = +big

       Bot_deltax = 2.0 * big
       Bot_deltay = 2.0 * big

       Bot( 1, 1 )%x = [ -big, -big, DepthB ]
       Bot( 1, 2 )%x = [ -big,  big, DepthB ]
       Bot( 2, 1 )%x = [  big, -big, DepthB ]
       Bot( 2, 2 )%x = [  big,  big, DepthB ]

       Bot( 1, 1 )%t  = [ 1.0, 0.0, 0.0 ]   ! tangent to bottom
       Bot( 1, 1 )%n1 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 1, 1 )%n2 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 1, 2 )%t  = [ 1.0, 0.0, 0.0 ]   ! tangent to bottom
       Bot( 1, 2 )%n1 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 1, 2 )%n2 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 2, 1 )%t  = [ 1.0, 0.0, 0.0 ]   ! tangent to bottom
       Bot( 2, 1 )%n1 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 2, 1 )%n2 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 2, 2 )%t  = [ 1.0, 0.0, 0.0 ]   ! tangent to bottom
       Bot( 2, 2 )%n1 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal
       Bot( 2, 2 )%n2 = [ 0.0, 0.0, 1.0 ]   ! outward-pointing normal

    END SELECT

    ! dummy BotSeg info to force GetBotSeg to search for the active segment on first call
    xBotSeg = [ +big, -big ]
    yBotSeg = [ +big, -big ]

  END SUBROUTINE ReadBTY3D

  ! **********************************************************************!

  SUBROUTINE GetTopSeg3D( x )

    ! Get the Top segment info (index and range interval) for range, r
    ! sets Topx and Topn

    INTEGER, PARAMETER :: PRTFile = 6
    REAL (KIND=8), INTENT( IN ) :: x( 3 )
    LOGICAL       :: Goodx = .TRUE., Goody = .TRUE.
    INTEGER       :: IsegTopT( 1 )

    ! x coordinate

    IF ( x( 1 ) < xTopSeg( 1 ) .OR. x( 1 ) > xTopSeg( 2 ) ) THEN

       IsegTopT = MAXLOC( Top( :, 1 )%x( 1 ), Top( :, 1 )%x( 1 ) < x( 1 ) )

       IF ( IsegTopT( 1 ) > 0 .AND. IsegTopT( 1 ) < NatiPts( 1 ) ) THEN  ! IsegTop MUST LIE IN [ 1, NatiPts-1 ]
          Goodx = .TRUE.
          IsegTopx = IsegTopT( 1 )
          xTopSeg  = [ Top( IsegTopx, 1 )%x( 1 ), Top( IsegTopx + 1, 1 )%x( 1 ) ]   ! segment limits in range
          Top_deltax = xTopSeg( 2 ) - xTopSeg( 1 )
       ELSE
          Goodx = .FALSE.
          WRITE( PRTFile, * ) 'x = ', x( 1 )
          WRITE( PRTFile, * ) 'xMin = ', Top( 1           , 1 )%x( 1 )
          WRITE( PRTFile, * ) 'xMax = ', Top( NatiPts( 1 ), 1 )%x( 1 )
          WRITE( PRTFile, * ) 'Warning in GetTopSeg3D : Altimetry undefined above the ray'
       ENDIF
    END IF

    ! y coordinate

    IF ( x( 2 ) < yTopSeg( 1 ) .OR. x( 2 ) > yTopSeg( 2 ) ) THEN

       IsegTopT = MAXLOC( Top( 1, : )%x( 2 ), Top( 1, : )%x( 2 ) < x( 2 ) )

       IF ( IsegTopT( 1 ) > 0 .AND. IsegTopT( 1 ) < NatiPts( 2 ) ) THEN  ! IsegTop MUST LIE IN [ 1, NatiPts-1 ]
          GoodY = .TRUE.
          IsegTopy = IsegTopT( 1 )
          yTopSeg  = [ Top( 1, IsegTopy )%x( 2 ), Top( 1, IsegTopy + 1 )%x( 2 ) ]   ! segment limits in range
          Top_deltay = yTopSeg( 2 ) - yTopSeg( 1 )
       ELSE
          GoodY = .FALSE.
          WRITE( PRTFile, * ) 'y = ', x( 2 )
          WRITE( PRTFile, * ) 'yMin = ', Top( 1, 1            )%x( 2 )
          WRITE( PRTFile, * ) 'yMax = ', Top( 1, NatiPts( 2 ) )%x( 2 )
          WRITE( PRTFile, * ) 'Warning in GetTopSeg3D : Altimetry undefined above the ray'
       ENDIF
    END IF

    IF ( Goodx .AND. Goody ) THEN
       Topx = Top( IsegTopx, IsegTopy )%x

       ! identify the normal based on the active triangle of a pair

       Top_tri_n = [ -Top_deltay, Top_deltax ]  ! normal of triangle side pointing up and to the left

       IF ( DOT_PRODUCT( x( 1 : 2 ) - Top( IsegTopx, IsegTopy )%x( 1 : 2 ), Top_tri_n ) < 0 ) THEN
          Topn = Top( IsegTopx, IsegTopy )%n1
       ELSE
          Topn = Top( IsegTopx, IsegTopy )%n2
       END IF

       ! if the Top depth is bad (a NaN) then set the segment flags to indicate that
       IF ( ISNAN( Topx( 3 ) ) .OR. ANY( ISNAN( Topn ) ) ) THEN
          IsegTopx = 0
          IsegTopy = 0
       END IF
    END IF
  END SUBROUTINE GetTopSeg3D

  ! **********************************************************************!

  SUBROUTINE GetBotSeg3D( x )

    ! Get the Bottom segment info (index and range interval) for range, r
    ! sets Botx and Botn

    INTEGER, PARAMETER :: PRTFile = 6
    REAL (KIND=8), INTENT( IN ) :: x( 3 )
    LOGICAL       :: Goodx = .TRUE., Goody = .TRUE.
    INTEGER       :: IsegBotT( 1 )

    ! x coordinate

    IF ( x( 1 ) < xBotSeg( 1 ) .OR. x( 1 ) > xBotSeg( 2 ) ) THEN   ! are we outside the segment?

       ! calculate index of bracketing segment
       IsegBotT = MAXLOC( Bot( :, 1 )%x( 1 ), Bot( :, 1 )%x( 1 ) < x( 1 ) )

!!$       ! The above MAXLOC is concise, but it's unnecessarily testing every segment
!!$       ! Only need to test adjacent segments; however, this doesn't seem to be much faster
!!$       ! we reached the same conclusion in SSPMod/Quad

!!$       IF ( xBotSeg( 1 ) == big ) THEN   ! for first call we do a full search in a lazy way
!!$          IsegBotT = MAXLOC( Bot( :, 1 )%x( 1 ), Bot( :, 1 )%x( 1 ) < x( 1 ) )
!!$       ELSE
!!$          ! search left
!!$          DO WHILE ( x( 1 ) < xBotSeg( 1 ) .AND. IsegBotx > 1 )
!!$             IsegBotx = IsegBotx - 1
!!$             xBotSeg  = [ Bot( IsegBotx, 1 )%x( 1 ), Bot( IsegBotx + 1, 1 )%x( 1 ) ]   ! segment limits in range
!!$          END DO
!!$
!!$          ! search right
!!$          DO WHILE ( x( 1 ) > xBotSeg( 2 ) .AND. IsegBotx < NbtyPts( 1 ) - 1 )
!!$             IsegBotx = IsegBotx + 1
!!$             xBotSeg  = [ Bot( IsegBotx, 1 )%x( 1 ), Bot( IsegBotx + 1, 1 )%x( 1 ) ]   ! segment limits in range
!!$          END DO
!!$
!!$          IsegBotT( 1 ) = IsegBotx
!!$       END IF
!!$
!!$       ! if we didn't find a bracketing segment set the segment to 0 indicating a failure
!!$       IsegBotx = IsegBotT( 1 )
!!$       xBotSeg  = [ Bot( IsegBotx, 1 )%x( 1 ), Bot( IsegBotx + 1, 1 )%x( 1 ) ] 
!!$       IF ( x( 1 ) < xBotSeg( 1 ) .OR. x( 1 ) > xBotSeg( 2 ) ) IsegBotT( 1 ) = 0

!!$       ! Here's yet another way to do the search for the segment
!!$       ! on the first call, need to do a full search
!!$       IF ( xBotSeg( 1 ) == big ) THEN
!!$          IsegBotT = MAXLOC( Bot( :, 1 )%x( 1 ), Bot( :, 1 )%x( 1 ) < x( 1 ) )
!!$       ELSE
!!$          DO
!!$             IF (      x( 1 ) < xBotSeg( 1 ) ) THEN
!!$                IF ( IsegBotx <= 0            ) EXIT
!!$                IsegBotx = IsegBotx - 1
!!$             ELSE IF ( x( 1 ) > xBotSeg( 2 ) ) THEN
!!$                IF ( IsegBotx >= NBtyPts( 1 ) ) EXIT
!!$                IsegBotx = IsegBotx + 1
!!$             ELSE
!!$                EXIT
!!$             END IF
!!$
!!$             xBotSeg  = [ Bot( IsegBotx, 1 )%x( 1 ), Bot( IsegBotx + 1, 1 )%x( 1 ) ]   ! segment limits in range
!!$          END DO
!!$
!!$          IsegBotT( 1 ) = IsegBotx
!!$       END IF

       IF ( IsegBotT( 1 ) > 0 .AND. IsegBotT( 1 ) < NbtyPts( 1 ) ) THEN  ! IsegBot MUST LIE IN [ 1, NbtyPts-1 ]
          Goodx = .TRUE.
          IsegBotx = IsegBotT( 1 )   
          xBotSeg  = [ Bot( IsegBotx, 1 )%x( 1 ), Bot( IsegBotx + 1, 1 )%x( 1 ) ]   ! segment limits in range
          Bot_deltax = xBotSeg( 2 ) - xBotSeg( 1 )
          !write( *, * ) 'New x', IsegBotx, xBotSeg
       ELSE
          Goodx = .FALSE.
          IsegBotx = 0
          WRITE( PRTFile, * ) 'x = ', x( 1 )
          WRITE( PRTFile, * ) 'xMin = ', Bot( 1           , 1 )%x( 1 )
          WRITE( PRTFile, * ) 'xMax = ', Bot( NbtyPts( 1 ), 1 )%x( 1 )
          WRITE( PRTFile, * ) 'Warning in GetBotSeg3D : Bathymetry undefined below the ray'
       ENDIF

    END IF

    ! y coordinate

    IF ( x( 2 ) < yBotSeg( 1 ) .OR. x( 2 ) > yBotSeg( 2 ) ) THEN

       IsegBotT = MAXLOC( Bot( 1, : )%x( 2 ), Bot( 1, : )%x( 2 ) < x( 2 ) )

       !write( *, * )
       !write( *, * ) 'ISeg  ', IsegBotT, Bot( 1, IsegBotT )%x( 2 )

!!$       ! The above MAXLOC is concise, but it's unnecessarily testing every segment
!!$       ! Only need to test adjacent segments; however, this doesn't seem to be much faster
!!$       ! on the first call, need to do a full search
!!$       IF ( yBotSeg( 1 ) == big ) THEN
!!$          IsegBotT = MAXLOC( Bot( 1, : )%x( 2 ), Bot( 1, : )%x( 2 ) < x( 2 ) )
!!$       ELSE
!!$          DO
!!$             IF (      x( 2 ) < yBotSeg( 1 ) ) THEN
!!$                IF ( IsegBoty <= 0            ) EXIT
!!$                IsegBoty = IsegBoty - 1
!!$             ELSE IF ( x( 2 ) > yBotSeg( 2 ) ) THEN
!!$                IF ( IsegBoty >= NBtyPts( 2 ) ) EXIT
!!$                IsegBoty = IsegBoty + 1
!!$             ELSE
!!$                EXIT
!!$             END IF
!!$             yBotSeg  = [ Bot( 1, IsegBoty )%x( 2 ), Bot( 1, IsegBoty + 1 )%x( 2 ) ]   ! segment limits in range
!!$          END DO
!!$
!!$          IsegBotT( 1 ) = IsegBoty
!!$       END IF

       IF ( IsegBotT( 1 ) > 0 .AND. IsegBotT( 1 ) < NbtyPts( 2 ) ) THEN  ! IsegBot MUST LIE IN [ 1, NbtyPts-1 ]
          Goody = .TRUE.
          IsegBoty = IsegBotT( 1 )
          yBotSeg  = [ Bot( 1, IsegBoty )%x( 2 ), Bot( 1, IsegBoty + 1 )%x( 2 ) ]   ! segment limits in range
          Bot_deltay = yBotSeg( 2 ) - yBotSeg( 1 )
          !write( *, * ) 'New y', IsegBoty, yBotSeg

       ELSE
          Goody = .FALSE.
          IsegBoty = 0
          WRITE( PRTFile, * ) 'y = ', x( 2 )
          WRITE( PRTFile, * ) 'yMin = ', Bot( 1, 1            )%x( 2 )
          WRITE( PRTFile, * ) 'yMax = ', Bot( 1, NbtyPts( 2 ) )%x( 2 )
          WRITE( PRTFile, * ) 'Warning in GetBotSeg3D : Bathymetry undefined below the ray'
       ENDIF

    END IF

    IF ( Goodx .AND. Goody ) THEN
       Botx = Bot( IsegBotx, IsegBoty )%x

       ! identify the normal based on the active triangle of a pair

       Bot_tri_n = [ -Bot_deltay, Bot_deltax ]  ! normal of triangle side pointing up and to the left

       IF ( DOT_PRODUCT( x( 1 : 2 ) - Bot( IsegBotx, IsegBoty )%x( 1 : 2 ), Bot_tri_n ) < 0 ) THEN
          Botn = Bot( IsegBotx, IsegBoty )%n1
       ELSE
          Botn = Bot( IsegBotx, IsegBoty )%n2
       END IF

       ! if the Bot depth is bad (a NaN) then set the segment flags to indicate that
       IF ( ISNAN( Botx( 3 ) ) .OR. ANY( ISNAN( Botn ) ) ) THEN
          IsegBotx = 0
          IsegBoty = 0
       END IF
    END IF
  END SUBROUTINE GetBotSeg3D

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

    INTEGER                          :: NPts( 2 ) = [ 0, 0 ]
    REAL      (KIND=8)               :: p1( 3 ), p2( 3 ), p3( 3 ), p4( 3 ), U( 3 ), V( 3 )
    REAL      (KIND=8)               :: n1( 3 ), n2( 3 )      ! normal vectors to the pair of triangles
    REAL      (KIND=8)               :: tvec( 3 ), Len
    TYPE(BdryPt)                     :: Bdry( :, : )
    CHARACTER (LEN=3),  INTENT( IN ) :: BotTop           ! Flag indicating bottom or top reflection
    CHARACTER (LEN=2)                :: CurvilinearFlag = '-'
    REAL      (KIND=8)               :: mx, my, n( 3 )

    SELECT CASE ( BotTop )
    CASE ( 'Bot' )
       NPts = NbtyPts
       CurvilinearFlag = btyType
    CASE ( 'Top' )
       NPts = NatiPts
       CurvilinearFlag = atiType
    END SELECT

    ! normals on triangle faces
    DO ix = 1, NPts( 1 ) - 1
       DO iy = 1, NPts( 2 ) - 1
          ! coordinates of corner nodes, moving counter-clockwise around the rectangle
          p1 = Bdry( ix,     iy     )%x
          p2 = Bdry( ix + 1, iy     )%x
          p3 = Bdry( ix + 1, iy + 1 )%x
          p4 = Bdry( ix,     iy + 1 )%x

          ! edges for triangle 1
          U = p2 - p1   ! tangent along one edge
          V = p3 - p1   ! tangent along another edge

          ! normal vector is the cross-product of the edge tangents
          n1( 1 ) = U( 2 ) * V( 3 ) - U( 3 ) * V( 2 )
          n1( 2 ) = U( 3 ) * V( 1 ) - U( 1 ) * V( 3 )
          n1( 3 ) = U( 1 ) * V( 2 ) - U( 2 ) * V( 1 )
          IF ( BotTop == 'Top' ) n1 = -n1

          Bdry( ix, iy )%n1 = n1 / NORM2( n1 )   ! scale to make it a unit normal

          ! edges for triangle 2
          U = p3 - p1   ! tangent along one edge
          V = p4 - p1   ! tangent along another edge

          ! normal vector is the cross-product of the edge tangents
          n2( 1 ) = U( 2 ) * V( 3 ) - U( 3 ) * V( 2 )
          n2( 2 ) = U( 3 ) * V( 1 ) - U( 1 ) * V( 3 )
          n2( 3 ) = U( 1 ) * V( 2 ) - U( 2 ) * V( 1 )
          IF ( BotTop == 'Top' ) n2 = -n2
          
          Bdry( ix, iy )%n2 = n2 / NORM2( n2 )   ! scale to make it a unit normal
       
       END DO
    END DO

    ! normals at nodes
    ! use forward, centered, or backward difference formulas
    DO ix = 1, NPts( 1 )
       DO iy = 1, NPts( 2 )
          IF ( ix == 1 ) THEN
             mx = ( Bdry( ix + 1, iy     )%x( 3 ) - Bdry( ix    , iy     )%x( 3 ) ) / &
                  ( Bdry( ix + 1, iy     )%x( 1 ) - Bdry( ix    , iy     )%x( 1 ) )
          ELSE IF ( ix == Npts( 1 ) ) THEN
             mx = ( Bdry( ix    , iy     )%x( 3 ) - Bdry( ix - 1, iy     )%x( 3 ) ) / &
                  ( Bdry( ix    , iy     )%x( 1 ) - Bdry( ix - 1, iy     )%x( 1 ) )
          ELSE
             mx = ( Bdry( ix + 1, iy     )%x( 3 ) - Bdry( ix - 1, iy     )%x( 3 ) ) / &
                  ( Bdry( ix + 1, iy     )%x( 1 ) - Bdry( ix - 1, iy     )%x( 1 ) )
          END IF

          IF ( iy == 1 ) THEN
             my = ( Bdry( ix    , iy + 1 )%x( 3 ) - Bdry( ix    , iy     )%x( 3 ) ) / &
                  ( Bdry( ix    , iy + 1 )%x( 2 ) - Bdry( ix    , iy     )%x( 2 ) )
          ELSE IF ( iy == Npts( 2 ) ) THEN
             my = ( Bdry( ix    , iy     )%x( 3 ) - Bdry( ix    , iy - 1 )%x( 3 ) ) / &
                  ( Bdry( ix    , iy     )%x( 2 ) - Bdry( ix    , iy - 1 )%x( 2 ) )
          ELSE
             my = ( Bdry( ix    , iy + 1 )%x( 3 ) - Bdry( ix    , iy - 1 )%x( 3 ) ) / &
                  ( Bdry( ix    , iy + 1 )%x( 2 ) - Bdry( ix    , iy - 1 )%x( 2 ) )
          END IF

          n = [ -mx, -my, 1.0D0 ]   ! this a normal to the surface

          IF ( ix < NPts( 1 ) .AND. iy < NPts( 2 ) ) THEN
             ! xx term
             Bdry( ix, iy )%phi_xx = atan2( n( 3 ), n( 1 ) )   ! this is the angle at each node

             ! xy term
             tvec = Bdry( ix + 1, iy + 1 )%x - Bdry( ix, iy )%x
             Len  = SQRT( tvec( 1 ) ** 2 + tvec( 2 ) ** 2 )
             tvec = tvec / Len
             Bdry( ix, iy )%phi_xy = atan2( n( 3 ), n( 1 ) * tvec( 1 ) + n( 2 ) * tvec( 2 ) )   ! this is the angle at each node

             ! yy term
             Bdry( ix, iy )%phi_yy = atan2( n( 3 ), n( 2 ) )   ! this is the angle at each node
          END IF

          Bdry( ix, iy )%Noden_unscaled = n
          Bdry( ix, iy )%Noden = n / NORM2( n )          
       END DO
    END DO

    IF ( CurvilinearFlag( 1 : 1 ) == 'C' ) THEN ! curvilinear option: compute derivative as centered difference between two nodes

       ! compute curvatures in each segment

       ! - sign below because the node normal = ( -mx, -my, 1 )
       DO ix = 1, NPts( 1 ) - 1
          DO iy = 1, NPts( 2 ) - 1
             ! z_xx (difference in x of z_x)

             Bdry( ix, iy )%z_xx = -( Bdry( ix + 1, iy     )%Noden_unscaled( 1 ) - Bdry( ix, iy )%Noden_unscaled( 1 ) ) / &
                                    ( Bdry( ix + 1, iy     )%x(              1 ) - Bdry( ix, iy )%x(              1 ) )

             tvec = Bdry( ix + 1, iy )%x - Bdry( ix, iy )%x
             Len  = SQRT( tvec( 1 ) ** 2 + tvec( 3 ) ** 2 )
             Bdry( ix, iy )%kappa_xx = ( Bdry( ix + 1, iy )%phi_xx - Bdry( ix, iy )%phi_xx ) / Len ! this is curvature = dphi/ds

             ! z_xy (difference in y of z_x)

             Bdry( ix, iy )%z_xy = -( Bdry( ix    , iy + 1 )%Noden_unscaled( 1 ) - Bdry( ix, iy )%Noden_unscaled( 1 ) ) / &
                                    ( Bdry( ix    , iy + 1 )%x(              2 ) - Bdry( ix, iy )%x(              2 ) )

             tvec = Bdry( ix + 1, iy + 1 )%x - Bdry( ix, iy )%x
             Len  = SQRT( tvec( 1 ) ** 2 + tvec( 2 ) ** 2 + tvec( 3 ) ** 2 )
             Bdry( ix, iy )%kappa_xy = ( Bdry( ix + 1, iy + 1 )%phi_xy - Bdry( ix, iy )%phi_xy ) / Len ! this is curvature = dphi/ds

             ! new
             tvec = Bdry( ix, iy + 1 )%x - Bdry( ix, iy )%x
             Len  = SQRT( tvec( 2 ) ** 2 + tvec( 3 ) ** 2 )
             Bdry( ix, iy )%kappa_xy = ( Bdry( ix, iy + 1 )%phi_xx - Bdry( ix, iy )%phi_xx ) / Len ! this is curvature = dphi/ds

             ! z_yy (difference in y of z_y)

             Bdry( ix, iy )%z_yy = -( Bdry( ix    , iy + 1 )%Noden_unscaled( 2 ) - Bdry( ix, iy )%Noden_unscaled( 2 ) ) / &
                                    ( Bdry( ix    , iy + 1 )%x(              2 ) - Bdry( ix, iy )%x(              2 ) )

             tvec = Bdry( ix, iy + 1 )%x - Bdry( ix, iy )%x
             Len  = SQRT( tvec( 2 ) ** 2 + tvec( 3 ) ** 2 )
             Bdry( ix, iy )%kappa_yy = ( Bdry( ix, iy + 1 )%phi_yy - Bdry( ix, iy )%phi_yy ) / Len ! this is curvature = dphi/ds

             ! introduce Len factor per Eq. 4.4.18 in Cerveny's book
             Len = NORM2( Bdry( ix, iy )%Noden_unscaled )
             Bdry( ix, iy )%z_xx = Bdry( ix, iy )%z_xx / Len
             Bdry( ix, iy )%z_xy = Bdry( ix, iy )%z_xy / Len
             Bdry( ix, iy )%z_yy = Bdry( ix, iy )%z_yy / Len
          END DO
       END DO
    ELSE
       Bdry%z_xx = 0
       Bdry%z_xy = 0
       Bdry%z_yy = 0

       Bdry%kappa_xx = 0
       Bdry%kappa_xy = 0
       Bdry%kappa_yy = 0
    END IF
!!$    write( *, * ) 'ix=1', Bdry( 1, : )%kappa_xx
!!$    write( *, * ) 'ix=1', Bdry( 1, : )%kappa_xy
!!$    write( *, * ) 'ix=1', Bdry( 1, : )%kappa_yy
!!$    write( *, * ) 'iy=1', Bdry( :, 1 )%kappa_xx
!!$    write( *, * ) 'iy=1', Bdry( :, 1 )%kappa_xy
!!$    write( *, * ) 'iy=1', Bdry( :, 1 )%kappa_yy
!!$    write( *, * ) 'D'
!!$    write( *, * ) Bdry( :, : )%z_xx
!!$    write( *, * ) Bdry( :, : )%z_xy
!!$    write( *, * ) Bdry( :, : )%z_yy
!!$
!!$    write( *, * ) 'kappa'
!!$    write( *, * ) Bdry( :, : )%kappa_xx
!!$    write( *, * ) Bdry( :, : )%kappa_xy
!!$    write( *, * ) Bdry( :, : )%kappa_yy

  END SUBROUTINE ComputeBdryTangentNormal

END MODULE bdry3Dmod

 
