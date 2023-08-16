PROGRAM FIELD3D

  ! Generate file of replica vectors
  ! Multiple receivers depths are handled inefficiently.

  USE SourceReceiverPositions
  USE BeamPattern
  USE ElementMod
  USE ReadModes
  USE EvaluateGBMod
  USE Evaluate3DMod
  USE EvaluatepdqMod
  USE RWSHDFile
  USE FatalError

  IMPLICIT NONE
  !INTEGER, PARAMETER   :: SHDFile = 25, MaxM = 2380, MaxSet = 40000 ! MaxSet = 112000 for Channel Islands
  INTEGER, PARAMETER   :: SHDFile = 25, MaxM = 7000, MaxSet = 40000 ! MaxSet = 112000 for Channel Islands
  INTEGER              :: M( MaxSet ), MS( 3 ), Msource, Mlimit, iostat, ifreq = 1
  INTEGER              :: iCorner1, isx, isy, isz1, irz1, itheta, NSets, iRec, IElementSource, INode
  REAL        (KIND=8) :: freq, freq0, atten = 0
  REAL        (KIND=8) :: Rmin, Rmax   ! Receiver ranges
  REAL        (KIND=8) :: TStart, TEnd
  !COMPLEX              :: k( MaxM, MaxSet ), phiR( MaxM, MaxSet ), kS( MaxM, 3 ), phiST( MaxM, 3 )                   
  COMPLEX              :: kS( MaxM, 3 ), phiST( MaxM, 3 )                   
  CHARACTER   (LEN= 7) :: Option
  CHARACTER   (LEN=80) :: Title, TitleEnv, FileRoot, SHDFileName, PRTFileName
  CHARACTER   (LEN=10) :: PlotType = 'TL        '
  COMPLEX, ALLOCATABLE :: k( :, : ), P( :, : ), phiS( :, :, : ), phiR( :, : )

  CALL CPU_TIME( Tstart )   ! log the start time

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )
  SHDFileName = TRIM( FileRoot ) // '.shd'
  PRTFileName = TRIM( FileRoot ) // '.prt'

  ! Open the print file
  ! OPEN( UNIT = PRTFile, FILE = 'field3d.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )
  OPEN( UNIT = PRTFile, FILE = PRTFileName, STATUS = 'UNKNOWN', IOSTAT = iostat )

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) '__________________________________________________________________________'
  WRITE( PRTFile, * ) 'Running FIELD3D'
  WRITE( PRTFile, * ) 'Sums modes, producing pressure'
  WRITE( PRTFile, * )

  CALL READIN( FileRoot, Title, Option, Mlimit, Rmin, Rmax )  ! Read in all the user input
  ALLOCATE( phiS( MaxM, Pos%NSz, 3 ) )
  ALLOCATE( k( MaxM, MaxSet ), phiR( MaxM, MaxSet ) )   ! Allocated dynamically avoid limit on large matrices (mcmodel=small)

  SBPFlag = Option( 7 : 7 )
  CALL ReadPAT( FileRoot, PRTFile )  ! Source Beam Pattern

  CALL BuildAdjacentElementTable

  ALLOCATE( P( Pos%Ntheta, Pos%NRr ) )
  Nfreq = 1000                   !!! maximum number of frequencies
  allocate( freqVec( Nfreq ) )   !!! must match dimension in getmodes.f90

  ! MAIN Loop: over receiver depth

  RcvrDepth: DO irz1 = 1, Pos%NRz
     CALL ReadAllModeSets( MaxSet, MaxM, Pos%rz( irz1 ), NSets, M, k, phiR )   ! Get modes at the receiver depth   

     Source_x: DO isx = 1, Pos%NSx

        Source_y: DO isy = 1, Pos%NSy
           IElementSource = IdentifySourceElement( Pos%sx( isx ), Pos%sy( isy ) )
           ! write( PRTFile, * ) 'source position, element', Pos%sx( isx ), Pos%sy( isy ), IElementSource

           ! Read modes at all the source depths, for the corner nodes of the source triangle

           Corner: DO ICorner1 = 1, 3 
              INode = Node( ICorner1, IElementSource ) 
              IF ( ModeFileName( INode )( 1 : 5 ) /= 'DUMMY' ) THEN
                 CALL GetModes( ModeFileName( INode ), 1, ifreq, MaxM, Pos%sz, Pos%NSz, 'N', &
                      kS( 1, ICorner1 ), phiS( 1, 1, ICorner1 ), MS( ICorner1 ), freqVec, Nfreq, TitleEnv )
              ELSE
                 MS( ICorner1 ) = 0   ! Dummy mode file, set number of modes to zero
              END IF
           END DO Corner

           Msource = MIN( MS( 1 ), MS( 2 ), MS( 3 ) )   ! limited by the # at the corners of the surrounding triangle

           Source_z: DO isz1 = 1, Pos%NSz 

              P = 0
              IF ( Msource > 0 ) THEN

                 phiST( 1 : Msource, : ) = phiS( 1 : Msource, isz1, : ) ! Get modes at the source depth for this particular index

                 ! Call the appropriate routine to evaluate the field
                 SELECT CASE ( Option(1:3) )
                 CASE ( 'STD' )  ! STandDard (ignores hor. refraction)
                    CALL Evaluate3D(  freq, k, phiR, phiST, M, MaxM, IElementSource, Pos%sx( isx ), Pos%sy( isy ), &
                         Pos%theta, Pos%Ntheta, Rmin, Rmax, Pos%NRr, Mlimit, P )       
                 CASE ( 'PDQ' )  ! Pretty Darn Quick (ignores and ignores)
                    CALL EvaluatePDQ( freq, k, phiR, phiST, M, MaxM, IElementSource, Pos%sx( isx ), Pos%sy( isy ), &
                         Pos%theta, Pos%Ntheta, Rmin, Rmax, Pos%NRr, Mlimit, P )       
                 CASE ( 'GBT' )  ! Gaussian Beam Tracing (include hor. refraction)
                    CALL EvaluateGB(  k, phiR, phiST, M, MaxM, IElementSource, Pos%sx( isx ), Pos%sy( isy ), &
                         Pos%theta, Pos%Ntheta, Rmin, Rmax, Pos%NRr, Mlimit, Option, P, freq )       
                 CASE DEFAULT
                    CALL ERROUT( 'FIELD3D', 'Unknown option' ) 
                 END SELECT
              END IF

              ! Write header (we do it here, because we don't know the frequency until the first mode set is read)
              IF ( irz1 == 1 .AND. isx == 1 .AND. isy == 1 .AND. isz1 == 1 ) THEN
                 CALL WriteHeader( SHDFileName, Title, freq0, atten, PlotType )
              END IF

              ! Write out the field

              RcvrBearing: DO itheta = 1, Pos%Ntheta
                 iRec = 10 + &
                      ( isx    - 1 ) * Pos%NSy * Pos%Ntheta * Pos%NSz * Pos%NRz + &
                      ( isy    - 1 ) *           Pos%Ntheta * Pos%NSz * Pos%NRz + &
                      ( itheta - 1 ) *                        Pos%NSz * Pos%NRz + &
                      ( isz1   - 1 ) *                        Pos%NRz + irz1
                 WRITE( SHDFile, REC = iRec ) P( itheta, 1 : Pos%NRr )
              END DO RcvrBearing
              FLUSH( SHDFile )

           END DO Source_z
        END DO Source_y
     END DO Source_x
  END DO RcvrDepth

  CALL CPU_TIME( Tend )   ! check elapsed time
  WRITE( PRTFile, '( ''CPU time ='', G15.3, ''s'' )' ) Tend - Tstart
  WRITE( PRTFile, * ) 'Field3D completed successfully'

  CLOSE( PRTFile )
  CLOSE( SHDFile )

CONTAINS

  !**********************************************************************C

  SUBROUTINE READIN( FileRoot, Title, Option, Mlimit, Rmin, Rmax )

    ! Reads in all the user input

    INTEGER, PARAMETER :: FLPFile = 5
    CHARACTER (LEN=80), INTENT( IN  ) :: FileRoot
    INTEGER,            INTENT( OUT ) :: Mlimit
    REAL      (KIND=8), INTENT( OUT ) :: Rmin, Rmax
    CHARACTER (LEN= 7), INTENT( OUT ) :: Option
    CHARACTER (LEN=80), INTENT( OUT ) :: Title
    INTEGER :: I, IAllocStat, iElt, iostat

    ! Open the field parameter file
    OPEN( UNIT = FLPFile, FILE = TRIM( FileRoot ) // '.flp', STATUS = 'OLD', IOSTAT = iostat, ACTION = 'READ' )
    IF ( IOSTAT /= 0 ) THEN   ! successful open?
       WRITE( PRTFile, * ) 'FLPFile = ', TRIM( FileRoot ) // '.flp'
       CALL ERROUT( 'FIELD3D - READIN', 'Unable to open the field parameter file' )
    END IF

    READ(  FLPFile, * ) Title 
    WRITE( PRTFile, * ) Title 

    READ(  FLPFile, * ) Option 
    WRITE( PRTFile, * ) 'Option = ', Option

    READ(  FLPFile, * ) Mlimit 
    WRITE( PRTFile, * ) 'Number of modes = ', Mlimit 
    Mlimit = MIN( Mlimit, MaxM ) 

    ! Read source/receiver information
    CALL ReadVector( Pos%NSx, Pos%Sx, 'Source   x-coordinates, Sx', 'km' )
    CALL ReadVector( Pos%NSy, Pos%Sy, 'Source   y-coordinates, Sy', 'km' )
    CALL ReadSzRz( 0.0, 1.0E6 )  ! source/rcvr depths, Sz, Rz
    CALL ReadRcvrRanges          ! receiver    ranges, Rr
    CALL ReadRcvrBearings        ! receiver angles for radials

    Rmin = Pos%rr(  1 )
    Rmax = Pos%rr( Pos%NRr )

    ! Read nodal coordinates
    READ(  FLPFile, * ) NNodes 
    WRITE( PRTFile, * ) 'NNodes = ', NNodes
    ALLOCATE( x( NNodes ), y( NNodes ), ModeFileName( NNodes ), Iset( NNodes ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) THEN 
       CALL ERROUT( 'FIELD3D - READIN', 'Too many nodes' ) 
    ENDIF

    NodeLoop: DO I = 1, NNodes 
       READ( FLPFile, * ) x( I ), y( I ), ModeFileName( I ) 
    END DO NodeLoop
 
    x = 1000.0 * x   ! convert km to m
    y = 1000.0 * y

    ! Read in element definitions
    READ(  FLPFile, * ) NElts 
    WRITE( PRTFile, * ) 'NElts  = ', NElts 
    ALLOCATE( AdjacentElement( 3, NElts ), Node( 3, NElts ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) THEN 
       CALL ERROUT( 'FIELD3D - READIN', 'Too many elements' ) 
    ENDIF

    EltLoop: DO iElt = 1, NElts 
       READ( FLPFile, * ) Node( 1 : 3, IELt ) 
    END DO EltLoop

    IF ( Option( 4 : 4 ) == 'T' ) THEN 
       WRITE( PRTFile, * ) 'Performing a Tesselation check' 
       CALL TesselationCheck
       WRITE( PRTFile, * ) 'Passed the Tesselation check' 
    ENDIF

  END SUBROUTINE READIN

  !**********************************************************************C

  SUBROUTINE BuildAdjacentElementTable

    ! Constructs a table AdjacentElement( iElt, iside ) which gives the
    ! element number that shares iside with element number iElt

    INTEGER iElt, iEltT, iside, isideT, Node1, Node1T, Node2, Node2T

    AdjacentElement = 0

    Element1: DO iElt = 1, NElts  ! Loop over each trianglular element

       Side1: DO iside = 1, 3 ! Loop over triangle sides

          IF ( AdjacentElement( iside, iElt ) == 0 ) THEN 

             Node1 = Node( ICorner( iside, 1 ), iElt ) 
             Node2 = Node( ICorner( iside, 2 ), iElt ) 

             ! Search other elements to find common side
             Element2: DO iEltT = 1, NElts 
                IF ( iEltT /= iElt ) THEN 
                   Side2: DO isideT = 1, 3 
                      Node1T = Node( ICorner( isideT, 1 ), iEltT ) 
                      Node2T = Node( ICorner( isideT, 2 ), iEltT ) 

                      ! Do iElt and iEltT share this side?
                      IF ( ( Node1 == Node1T .AND. Node2 == Node2T ) .OR. &
                           ( Node1 == Node2T .AND. Node2 == Node1T ) ) THEN           
                         AdjacentElement( iside,  iElt  ) = iEltT 
                         AdjacentElement( isideT, iEltT ) = iElt 
                         CYCLE Side1
                      ENDIF

                   END DO Side2
                ENDIF
             END DO Element2
          ENDIF

       END DO Side1

    END DO Element1

  END SUBROUTINE BuildAdjacentElementTable

  !**********************************************************************C

  FUNCTION IdentifySourceElement( xs, ys ) 

    !     Identifies the element containing the source at ( xs, ys )          

    !     We define a function Enclosure( xs, ys ) which is 1 when ( xs, ys )     
    !     is inside a given triangle and decreases from 1 the further the   
    !     source moves away from the triangle.

    !     The element with the highest value of Enclosure is identified as the source element.
    !     If several elements enclose, the highest numbered element is given posession.

    REAL (KIND=8), INTENT( IN ) :: xs, ys
    INTEGER            :: iElt, Node1, Node2, Node3, IdentifySourceElement
    REAL (KIND=8)      :: x1, x2, x3, y1, y2, y3, A1, A2, A3, Delta, Enclosure, EnclosureMax

    IdentifySourceElement = 0
    EnclosureMax = 0.0 

    Element: DO iElt = 1, NElts 
       Node1 = Node( 1, iElt )
       Node2 = Node( 2, iElt ) 
       Node3 = Node( 3, iElt ) 

       x1 = x( Node1 )   ;   y1 = y( Node1 ) 
       x2 = x( Node2 )   ;   y2 = y( Node2 ) 
       x3 = x( Node3 )   ;   y3 = y( Node3 ) 

       ! Compute areas of triangles
       Delta = ( x2 * y3 - y2 * x3 ) - ( x1 * y3 - y1 * x3 ) + ( x1 * y2 - y1 * x2 )
       A1    = ( x2 * y3 - y2 * x3 ) - ( xs * y3 - ys * x3 ) + ( xs * y2 - ys * x2 )
       A2    = ( xs * y3 - ys * x3 ) - ( x1 * y3 - y1 * x3 ) + ( x1 * ys - y1 * xs )
       A3    = ( x2 * ys - y2 * xs ) - ( x1 * ys - y1 * xs ) + ( x1 * y2 - y1 * x2 )

       Enclosure = ABS( Delta ) / ( ABS( A1 ) + ABS( A2 ) + ABS( A3 ) ) 

       IF ( Enclosure > EnclosureMax ) THEN 
          IdentifySourceElement  = iElt 
          EnclosureMax = Enclosure 
       ENDIF

    END DO Element

    IF ( IdentifySourceElement == 0 ) CALL ERROUT( 'FIELD3D - IdentifySourceElement', 'Source not inside the grid' ) 

  END FUNCTION IdentifySourceElement

  !**********************************************************************C

  SUBROUTINE TesselationCheck 

    !     Checks to see that triangulation is a tesselation                 
    !     that is, that there are no overlapping triangles                  
    !     (holes or triangles inside triangles are still possible)          

    INTEGER  :: iElt1, iElt2, iside1, iside2
    REAL     :: x1, x2, x3, x4, y1, y2, y3, y4, Ux, Uy, Vx, Vy, Wx, Wy, Delta, S1, S2

    ! ICorner maps a side (1, 2 or 3) and a local node (1 or 2) to a
    ! corner (1, 2, OR 3) of the triangle

    Element1: DO iElt1 = 1, NElts-1 
       Side1: DO iside1 = 1, 3 
          x1 = x( Node( ICorner( iside1, 1 ), iElt1 ) ) 
          y1 = y( Node( ICorner( iside1, 1 ), iElt1 ) ) 

          x2 = x( Node( ICorner( iside1, 2 ), iElt1 ) )
          y2 = y( Node( ICorner( iside1, 2 ), iElt1 ) ) 

          Ux = x2 - x1   ;   Uy = y2 - y1 

          Element2: DO iElt2 = iElt1 + 1, NElts 
             Side2: DO iside2 = 1, 3 
                x3 = x( Node( ICorner( iside2, 1 ), iElt2 ) ) 
                y3 = y( Node( ICorner( iside2, 1 ), iElt2 ) ) 

                x4 = x( Node( ICorner( iside2, 2 ), iElt2 ) ) 
                y4 = y( Node( ICorner( iside2, 2 ), iElt2 ) )

                Vx = x4 - x3
                Vy = y4 - y3

                Wx = x4 - x2
                Wy = y4 - y2

                Delta = Ux * Vy - Uy * Vx 

                IF ( ABS( Delta ) > MAX( ABS( Ux * Uy ), ABS( Vx * Vy ) ) / 10000.0 ) THEN 
                   S1 = ( Ux * Wy - Uy * Wx ) / Delta 
                   S2 = ( Vx * Wy - Vy * Wx ) / Delta 

                   IF ( 0.001 < S1 .AND. S1 < 0.999 .AND.       &
                        &   0.001 < S2 .AND. S2 < 0.999 ) THEN      
                      WRITE( PRTFile, * ) 'Sides cross for elts = ', iElt1, iElt2      
                      WRITE( PRTFile, * ) 'Side 1:' 
                      WRITE( PRTFile, * ) '   (', x1, y1, ')' 
                      WRITE( PRTFile, * ) '   (', x2, y2, ')' 
                      WRITE( PRTFile, * ) '   S = ', S1 

                      WRITE( PRTFile, * ) 'Side 2:' 
                      WRITE( PRTFile, * ) '   (', x3, y3, ')' 
                      WRITE( PRTFile, * ) '   (', x4, y4, ')' 
                      WRITE( PRTFile, * ) '   S = ', S2 
                      STOP 
                   ENDIF
                ENDIF

             END DO Side2
          END DO Element2

       END DO Side1
    END DO Element1

  END SUBROUTINE TesselationCheck

  !**********************************************************************C

  SUBROUTINE ReadAllModeSets( MaxSet, MaxM, rz, NSets, M, k, phiR )      

    ! Reads in the values of the modes at the given receiver depth
    ! for every node in the triangulation

    INTEGER, INTENT( IN  ) :: MaxSet, MaxM
    INTEGER, INTENT( OUT ) :: M( * ), NSets
    COMPLEX, INTENT( OUT ) :: k( MaxM, * ), phiR( MaxM, * ) 
    REAL,    INTENT( IN  ) :: rz( 1 )
    INTEGER                :: INode, JNode, ifreq = 1, Nfreq
    REAL  (KIND=8)         :: freqVec( 1000 )      ! !!! frequency for which modes were calculated
    CHARACTER     (LEN=80) :: TitleEnv

    NSets = 0 

    NodeLoop: DO INode = 1, NNodes 
       ! Check if the modes have already been read
       IF ( INode >= 2 ) THEN 
          DO JNode = 1, INode-1 
             IF ( ModeFileName( INode ) == ModeFileName( JNode ) ) THEN 
                !write( PRTFile, * ) 'Copy previously read modes'
                Iset( INode ) = Iset( JNode )  ! Copy previously read modes
                CYCLE NodeLoop
             ENDIF
          END DO
       ENDIF

       NSets = NSets + 1 
       IF ( NSets > MaxSet ) THEN 
          WRITE( PRTFile, * ) 'MaxSet = ', MaxSet 
          CALL ERROUT( 'FIELD3D - ReadAllModeSets', 'Too many mode sets' )               
       ENDIF
       Iset( INode ) = NSets 

       ! Check for 'DUMMY' elts (acoustic absorbers)
       IF ( ModeFileName( INode )( 1 : 5 ) == 'DUMMY' ) THEN 
          M( NSets ) = 0 
       ELSE  ! Necessary to read in modes
          CALL GetModes( ModeFileName( INode ), 1, ifreq, MaxM, rz, 1, 'N', k( 1, NSets ), phiR( 1, NSets ), M( NSets ), &
               freqVec, Nfreq, TitleEnv )
       ENDIF
       FLUSH( PRTFile )
    END DO NodeLoop

    ! WRITE( PRTFile, * ) 'Number of distinct sets of modes = ', NSets

  END SUBROUTINE ReadAllModeSets

END PROGRAM FIELD3D
