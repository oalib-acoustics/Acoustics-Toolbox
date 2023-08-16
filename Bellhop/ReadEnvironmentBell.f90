MODULE ReadEnvironmentBell

  ! mbp 12/2018, based on much older subroutine

  USE BellhopMod
  USE MathConstants
  USE sspmod
  USE AttenMod
  USE FatalError

  IMPLICIT NONE

CONTAINS
  SUBROUTINE ReadEnvironment( FileRoot, ThreeD )

    ! Routine to read in and echo all the input data
    ! Note that default values of SSP, DENSITY, Attenuation will not work

    USE anglemod
    USE SourceReceiverPositions

    REAL      (KIND=8), PARAMETER   :: c0 = 1500.0
    LOGICAL,            INTENT(IN ) :: ThreeD
    CHARACTER (LEN=80), INTENT(IN ) :: FileRoot
    INTEGER            :: NPts, NMedia, iostat
    REAL               :: ZMin, ZMax
    REAL      (KIND=8) :: x( 2 ), c, cimag, gradc( 2 ), crr, crz, czz, rho, sigma, Depth
    CHARACTER (LEN= 2) :: AttenUnit
    CHARACTER (LEN=10) :: PlotType

    WRITE( PRTFile, * ) 'BELLHOP/BELLHOP3D'
    WRITE( PRTFile, * )

    ! Open the environmental file
    OPEN( UNIT = ENVFile, FILE = TRIM( FileRoot ) // '.env', STATUS = 'OLD', IOSTAT = iostat, ACTION = 'READ' )
    IF ( IOSTAT /= 0 ) THEN   ! successful open?
       WRITE( PRTFile, * ) 'ENVFile = ', TRIM( FileRoot ) // '.env'
       CALL ERROUT( 'BELLHOP - READIN', 'Unable to open the environmental file' )
    END IF

    ! Prepend model name to title
    IF ( ThreeD ) THEN
       Title( 1 : 11 ) = 'BELLHOP3D- '
       READ(  ENVFile, * ) Title( 12 : 80 )
    ELSE
       Title( 1 :  9 ) = 'BELLHOP- '
       READ(  ENVFile, * ) Title( 10 : 80 )
    END IF

    WRITE( PRTFile, * ) Title

    READ(  ENVFile, *    ) freq
    WRITE( PRTFile, '('' frequency = '', G11.4, '' Hz'', / )' ) freq

    READ(  ENVFile, * ) NMedia
    WRITE( PRTFile, * ) 'Dummy parameter NMedia = ', NMedia
    IF ( NMedia /= 1 ) CALL ERROUT( 'READIN', &
         'Only one medium or layer is allowed in BELLHOP; sediment layers must be handled using a reflection coefficient' )

    CALL ReadTopOpt( Bdry%Top%HS%Opt, Bdry%Top%HS%BC, AttenUnit, FileRoot )

    ! *** Top BC ***

    IF ( Bdry%Top%HS%BC == 'A' ) THEN
       WRITE( PRTFile, "( //, '      z         alphaR      betaR     rho        alphaI     betaI'    )" )
       WRITE( PRTFile, "(     '     (m)         (m/s)      (m/s)   (g/cm^3)      (m/s)     (m/s)', / )" )
    END IF

    CALL TopBot( freq, AttenUnit, Bdry%Top%HS )

    ! ****** Read in ocean SSP data ******

    READ(  ENVFile, * ) NPts, Sigma, Bdry%Bot%HS%Depth
    WRITE( PRTFile, * )
    WRITE( PRTFile, FMT = "( ' Depth = ', F10.2, ' m' )" ) Bdry%Bot%HS%Depth

    IF ( Bdry%Top%HS%Opt( 1 : 1 ) == 'A' ) THEN
       WRITE( PRTFile, * ) 'Analytic SSP option'
       ! following is hokey, should be set in Analytic routine
       SSP%NPts = 2
       SSP%z( 1 ) = 0.0
       SSP%z( 2 ) = Bdry%Bot%HS%Depth
    ELSE
       x = [ 0.0D0, Bdry%Bot%HS%Depth ]   ! tells SSP Depth to read to
       CALL EvaluateSSP( x, c, cimag, gradc, crr, crz, czz, rho, freq, 'INI' )
    ENDIF

    Bdry%Top%HS%Depth = SSP%z( 1 )   ! Depth of top boundary is taken from first SSP point
    ! bottom depth should perhaps be set the same way?

    ! *** Bottom BC ***
    Bdry%Bot%HS%Opt = '  '   ! initialize to blanks
    READ(  ENVFile, * ) Bdry%Bot%HS%Opt, Sigma
    WRITE( PRTFile, * )
    WRITE( PRTFile, FMT = "(33X, '( RMS roughness = ', G10.3, ' )' )" ) Sigma

    SELECT CASE ( Bdry%Bot%HS%Opt( 2 : 2 ) )
    CASE ( '~', '*' )
       WRITE( PRTFile, * ) '    Bathymetry file selected'
    CASE( '-', '_', ' ' )
    CASE DEFAULT
       CALL ERROUT( 'READIN', 'Unknown bottom option letter in second position' )
    END SELECT

    Bdry%Bot%HS%BC = Bdry%Bot%HS%Opt( 1 : 1 )
    CALL TopBot( freq, AttenUnit, Bdry%Bot%HS )

    ! *** source and receiver locations ***

    CALL ReadSxSy( ThreeD )     ! Read source/receiver x-y coordinates

    ZMin = SNGL( Bdry%Top%HS%Depth )
    ZMax = SNGL( Bdry%Bot%HS%Depth )
    ! CALL ReadSzRz( ZMin + 100 * SPACING( ZMin ), ZMax - 100 * SPACING( ZMax ) )   ! not sure why I had this
    CALL ReadSzRz( ZMin, ZMax )
    CALL ReadRcvrRanges
    IF ( ThreeD ) CALL ReadRcvrBearings
    CALL ReadfreqVec( freq,  Bdry%Top%HS%Opt( 6:6 ) )
    CALL ReadRunType( Beam%RunType, PlotType )

    Depth = Zmax - Zmin   ! water depth
    CALL ReadRayElevationAngles( freq, Depth, Bdry%Top%HS%Opt, Beam%RunType )
    IF ( ThreeD ) CALL ReadRayBearingAngles(   freq, Bdry%Top%HS%Opt, Beam%RunType )

    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) '__________________________________________________________________________'
    WRITE( PRTFile, * )

    ! Limits for tracing beams

    IF ( ThreeD ) THEN
       READ(  ENVFile, * ) Beam%deltas, Beam%Box%x, Beam%Box%y, Beam%Box%z
       Beam%Box%x = 1000.0 * Beam%Box%x   ! convert km to m
       Beam%Box%y = 1000.0 * Beam%Box%y   ! convert km to m

       IF ( Beam%deltas == 0.0 ) Beam%deltas = ( Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth ) / 10.0   ! Automatic step size selection
       ! WRITE( PRTFile, '('' frequency = '', G11.4, '' Hz'', / )' ) freq

       WRITE( PRTFile, * )
       WRITE( PRTFile, fmt = '(  '' Step length,       deltas = '', G11.4, '' m'' )' ) Beam%deltas
       WRITE( PRTFile, * )
       WRITE( PRTFile, fmt = '(  '' Maximum ray x-range, Box%x = '', G11.4, '' m'' )' ) Beam%Box%x
       WRITE( PRTFile, fmt = '(  '' Maximum ray y-range, Box%y = '', G11.4, '' m'' )' ) Beam%Box%y
       WRITE( PRTFile, fmt = '(  '' Maximum ray z-range, Box%z = '', G11.4, '' m'' )' ) Beam%Box%z
    ELSE
       READ(  ENVFile, * ) Beam%deltas, Beam%Box%z, Beam%Box%r

       WRITE( PRTFile, * )
       WRITE( PRTFile, fmt = '(  '' Step length,       deltas = '', G11.4, '' m'' )' ) Beam%deltas
       WRITE( PRTFile, * )
       WRITE( PRTFile, fmt = '(  '' Maximum ray depth, Box%z  = '', G11.4, '' m'' )' ) Beam%Box%z
       WRITE( PRTFile, fmt = '(  '' Maximum ray range, Box%r  = '', G11.4, ''km'' )' ) Beam%Box%r

       Beam%Box%r = 1000.0 * Beam%Box%r   ! convert km to m
    END IF

    ! *** Beam characteristics ***

    Beam%Type( 4 : 4 ) = Beam%RunType( 7 : 7 )   ! selects beam shift option
          
    SELECT CASE ( Beam%Type( 4 : 4 ) )
    CASE ( 'S' )
       WRITE( PRTFile, * ) 'Beam shift in effect'
    CASE DEFAULT
       WRITE( PRTFile, * ) 'No beam shift in effect'
    END SELECT

    IF ( Beam%RunType( 1 : 1 ) /= 'R' ) THEN   ! no worry about the beam type if this is a ray trace run

       ! Beam%Type( 1 : 1 ) is
       !   'G' or '^' Geometric hat beams in Cartesian coordinates
       !   'g' Geometric hat beams in ray-centered coordinates
       !   'B' Geometric Gaussian beams in Cartesian coordinates
       !   'b' Geometric Gaussian beams in ray-centered coordinates
       !   'S' Simple Gaussian beams
       !   'C' Cerveny Gaussian beams in Cartesian coordinates
       !   'R' Cerveny Gaussian beams in Ray-centered coordinates
       ! Beam%Type( 2 : 2 ) controls the setting of the beam width
       !   'F' space Filling
       !   'M' minimum width
       !   'W' WKB beams
       ! Beam%Type( 3 : 3 ) controls curvature changes on boundary reflections
       !   'D' Double
       !   'S' Single
       !   'Z' Zero
       ! Beam%Type( 4 : 4 ) selects whether beam shifts are implemented on boundary reflection
       !   'S' yes
       !   'N' no

       ! Curvature change can cause overflow in grazing case
       ! Suppress by setting BeamType( 3 : 3 ) = 'Z'

       Beam%Type( 1 : 1 ) = Beam%RunType( 2 : 2 )
       SELECT CASE ( Beam%Type( 1 : 1 ) )
       CASE ( 'G', 'g' , '^', 'B', 'b', 'S' )   ! geometric hat beams, geometric Gaussian beams, or simple Gaussian beams
       CASE ( 'R', 'C' )   ! Cerveny Gaussian Beams; read extra lines to specify the beam options
          READ(  ENVFile, * ) Beam%Type( 2 : 3 ), Beam%epsMultiplier, Beam%rLoop
          WRITE( PRTFile, * )
          WRITE( PRTFile, * )
          WRITE( PRTFile, * ) 'Type of beam = ', Beam%Type( 1 : 1 )

          SELECT CASE ( Beam%Type( 3 : 3 ) )
          CASE ( 'D' )
             WRITE( PRTFile, * ) 'Curvature doubling invoked'
          CASE ( 'Z' )
             WRITE( PRTFile, * ) 'Curvature zeroing invoked'
          CASE ( 'S' )
             WRITE( PRTFile, * ) 'Standard curvature condition'
          CASE DEFAULT
             CALL ERROUT( 'READIN', 'Unknown curvature condition' )
          END SELECT

          WRITE( PRTFile, * ) 'Epsilon multiplier', Beam%epsMultiplier
          WRITE( PRTFile, * ) 'Range for choosing beam width', Beam%rLoop

          ! Images, windows
          READ(  ENVFile, * ) Beam%Nimage, Beam%iBeamWindow, Beam%Component
          WRITE( PRTFile, * )
          WRITE( PRTFile, * ) 'Number of images, Nimage  = ', Beam%Nimage
          WRITE( PRTFile, * ) 'Beam windowing parameter  = ', Beam%iBeamWindow
          WRITE( PRTFile, * ) 'Component                 = ', Beam%Component
       CASE DEFAULT
          CALL ERROUT( 'READIN', 'Unknown beam type (second letter of run type)' )
       END SELECT
    END IF

    WRITE( PRTFile, * )
    CLOSE( ENVFile )

  END SUBROUTINE ReadEnvironment

  !**********************************************************************!

  SUBROUTINE ReadTopOpt( TopOpt, BC, AttenUnit, FileRoot )

    CHARACTER (LEN= 6), INTENT( OUT ) :: TopOpt
    CHARACTER (LEN= 1), INTENT( OUT ) :: BC                     ! Boundary condition type
    CHARACTER (LEN= 2), INTENT( OUT ) :: AttenUnit
    CHARACTER (LEN=80), INTENT( IN  ) :: FileRoot
    INTEGER            :: iostat

    TopOpt = '      '   ! initialize to blanks
    READ(  ENVFile, * ) TopOpt
    WRITE( PRTFile, * )

    SSP%Type  = TopOpt( 1 : 1 )
    BC        = TopOpt( 2 : 2 )
    AttenUnit = TopOpt( 3 : 4 )
    SSP%AttenUnit = AttenUnit

    ! SSP approximation options

    SELECT CASE ( SSP%Type )
    CASE ( 'N' )
       WRITE( PRTFile, * ) '    N2-linear approximation to SSP'
    CASE ( 'C' )
       WRITE( PRTFile, * ) '    C-linear approximation to SSP'
    CASE ( 'P' )
       WRITE( PRTFile, * ) '    PCHIP approximation to SSP'
    CASE ( 'S' )
       WRITE( PRTFile, * ) '    Spline approximation to SSP'
    CASE ( 'Q' )
       WRITE( PRTFile, * ) '    Quad approximation to SSP'
       OPEN ( FILE = TRIM( FileRoot ) // '.ssp', UNIT = SSPFile, FORM = 'FORMATTED', STATUS = 'OLD', IOSTAT = iostat )
       IF ( IOSTAT /= 0 ) THEN   ! successful open?
          WRITE( PRTFile, * ) 'SSPFile = ', TRIM( FileRoot ) // '.ssp'
          CALL ERROUT( 'BELLHOP - READIN', 'Unable to open the SSP file' )
       END IF
    CASE ( 'H' )
       WRITE( PRTFile, * ) '    Hexahedral approximation to SSP'
       OPEN ( FILE = TRIM( FileRoot ) // '.ssp', UNIT = SSPFile, FORM = 'FORMATTED', STATUS = 'OLD', IOSTAT = iostat )
       IF ( IOSTAT /= 0 ) THEN   ! successful open?
          WRITE( PRTFile, * ) 'SSPFile = ', TRIM( FileRoot ) // '.ssp'
          CALL ERROUT( 'BELLHOP - READIN', 'Unable to open the SSP file' )
       END IF
    CASE ( 'A' )
       WRITE( PRTFile, * ) '    Analytic SSP option'
    CASE DEFAULT
       CALL ERROUT( 'READIN', 'Unknown option for SSP approximation' )
    END SELECT

    ! Attenuation options

    SELECT CASE ( AttenUnit( 1 : 1 ) )
    CASE ( 'N' )
       WRITE( PRTFile, * ) '    Attenuation units: nepers/m'
    CASE ( 'F' )
       WRITE( PRTFile, * ) '    Attenuation units: dB/mkHz'
    CASE ( 'M' )
       WRITE( PRTFile, * ) '    Attenuation units: dB/m'
    CASE ( 'W' )
       WRITE( PRTFile, * ) '    Attenuation units: dB/wavelength'
    CASE ( 'Q' )
       WRITE( PRTFile, * ) '    Attenuation units: Q'
    CASE ( 'L' )
       WRITE( PRTFile, * ) '    Attenuation units: Loss parameter'
    CASE DEFAULT
       CALL ERROUT( 'READIN', 'Unknown attenuation units' )
    END SELECT

    ! optional addition of volume attenuation using standard formulas

    SELECT CASE ( AttenUnit( 2 : 2 ) )
    CASE ( 'T' )
       WRITE( PRTFile, * ) '    THORP volume attenuation added'
    CASE ( 'F' )
       WRITE( PRTFile, * ) '    Francois-Garrison volume attenuation added'
       READ(  ENVFile, * ) T, Salinity, pH, z_bar
       WRITE( PRTFile, "( ' T = ', G11.4, 'degrees   S = ', G11.4, ' psu   pH = ', G11.4, ' z_bar = ', G11.4, ' m' )" ) &
            T, Salinity, pH, z_bar
    CASE ( 'B' )
       WRITE( PRTFile, * ) '    Biological attenaution'
       READ( ENVFile, *  ) NBioLayers
       WRITE( PRTFile, * ) '      Number of Bio Layers = ', NBioLayers

       DO iBio = 1, NBioLayers
          READ( ENVFile, *  ) bio( iBio )%Z1, bio( iBio )%Z2, bio( iBio )%f0, bio( iBio )%Q, bio( iBio )%a0
          WRITE( PRTFile, * ) '      Top    of layer = ', bio( iBio )%Z1, ' m'
          WRITE( PRTFile, * ) '      Bottom of layer = ', bio( iBio )%Z2, ' m'
          WRITE( PRTFile, * ) '      Resonance frequency = ', bio( iBio )%f0, ' Hz'
          WRITE( PRTFile, * ) '      Q  = ', bio( iBio )%Q
          WRITE( PRTFile, * ) '      a0 = ', bio( iBio )%a0
       END DO
    CASE ( ' ' )
    CASE DEFAULT
       CALL ERROUT( 'READIN', 'Unknown top option letter in fourth position' )
    END SELECT

    SELECT CASE ( TopOpt( 5 : 5 ) )
    CASE ( '~', '*' )
       WRITE( PRTFile, * ) '    Altimetry file selected'
    CASE ( '-', '_', ' ' )
    CASE DEFAULT
       CALL ERROUT( 'READIN', 'Unknown top option letter in fifth position' )
    END SELECT

    SELECT CASE ( TopOpt( 6 : 6 ) )
    CASE ( 'I' )
       WRITE( PRTFile, * ) '    Development options enabled'
    CASE ( ' ' )
    CASE DEFAULT
       CALL ERROUT( 'READIN', 'Unknown top option letter in sixth position' )
    END SELECT

  END SUBROUTINE ReadTopOpt

  !**********************************************************************!

  SUBROUTINE ReadRunType( RunType, PlotType )

    ! Read the RunType variable and echo with explanatory information to the print file

    USE SourceReceiverPositions

    CHARACTER (LEN= 7), INTENT( OUT ) :: RunType
    CHARACTER (LEN=10), INTENT( OUT ) :: PlotType

    READ(  ENVFile, * ) RunType
    WRITE( PRTFile, * )

    SELECT CASE ( RunType( 1 : 1 ) )
    CASE ( 'R' )
       WRITE( PRTFile, * ) 'Ray trace run'
    CASE ( 'E' )
       WRITE( PRTFile, * ) 'Eigenray trace run'
    CASE ( 'I' )
       WRITE( PRTFile, * ) 'Incoherent TL calculation'
    CASE ( 'S' )
       WRITE( PRTFile, * ) 'Semi-coherent TL calculation'
    CASE ( 'C' )
       WRITE( PRTFile, * ) 'Coherent TL calculation'
    CASE ( 'A' )
       WRITE( PRTFile, * ) 'Arrivals calculation, ASCII  file output'
    CASE ( 'a' )
       WRITE( PRTFile, * ) 'Arrivals calculation, binary file output'
    CASE DEFAULT
       CALL ERROUT( 'READIN', 'Unknown RunType selected' )
    END SELECT

    SELECT CASE ( RunType( 2 : 2 ) )
    CASE ( 'C' )
       WRITE( PRTFile, * ) 'Cartesian beams'
    CASE ( 'R' )
       WRITE( PRTFile, * ) 'Ray centered beams'
    CASE ( 'S' )
       WRITE( PRTFile, * ) 'Simple gaussian beams'
    CASE ( 'b' )
       WRITE( PRTFile, * ) 'Geometric gaussian beams in ray-centered coordinates'
    CASE ( 'B' )
       WRITE( PRTFile, * ) 'Geometric gaussian beams in Cartesian coordinates'
    CASE ( 'g' )
       WRITE( PRTFile, * ) 'Geometric hat beams in ray-centered coordinates'
    CASE DEFAULT
       RunType( 2 : 2 ) = 'G'
       WRITE( PRTFile, * ) 'Geometric hat beams in Cartesian coordinates'
    END SELECT

    SELECT CASE ( RunType( 4 : 4 ) )
    CASE ( 'R' )
       WRITE( PRTFile, * ) 'Point source (cylindrical coordinates)'
    CASE ( 'X' )
       WRITE( PRTFile, * ) 'Line source (Cartesian coordinates)'
    CASE DEFAULT
       RunType( 4 : 4 ) = 'R'
       WRITE( PRTFile, * ) 'Point source (cylindrical coordinates)'
    END SELECT

    SELECT CASE ( RunType( 5 : 5 ) )
    CASE ( 'R' )
       WRITE( PRTFile, * ) 'Rectilinear receiver grid: Receivers at ( Rr( ir ), Rz( ir ) ) )'
       PlotType = 'rectilin  '
    CASE ( 'I' )
       WRITE( PRTFile, * ) 'Irregular grid: Receivers at Rr( : ) x Rz( : )'
       IF ( Pos%NRz /= Pos%NRr ) CALL ERROUT( 'READIN', 'Irregular grid option selected with NRz not equal to Nr' )
       PlotType = 'irregular '
    CASE DEFAULT
       WRITE( PRTFile, * ) 'Rectilinear receiver grid: Receivers at Rr( : ) x Rz( : )'
       RunType( 5 : 5 ) = 'R'
       PlotType = 'rectilin  '
    END SELECT

    SELECT CASE ( RunType( 6 : 6 ) )
    CASE ( '2' )
       WRITE( PRTFile, * ) 'N x 2D calculation (neglects horizontal refraction)'
    CASE ( '3' )
       WRITE( PRTFile, * ) '3D calculation'
    CASE DEFAULT
       RunType( 6 : 6 ) = '2'
    END SELECT

  END SUBROUTINE ReadRunType

  !**********************************************************************!

  SUBROUTINE TopBot( freq, AttenUnit, HS )

    ! Handles top and bottom boundary conditions

    REAL     (KIND=8), INTENT( IN    ) :: freq               ! frequency
    CHARACTER (LEN=2), INTENT( IN    ) :: AttenUnit
    TYPE( HSInfo ),    INTENT( INOUT ) :: HS
    REAL     (KIND=8)                  :: Mz, vr, alpha2_f          ! values related to grain size

    ! Echo to PRTFile user's choice of boundary condition

    SELECT CASE ( HS%BC )
    CASE ( 'V' )
       WRITE( PRTFile, * ) '    VACUUM'
    CASE ( 'R' )
       WRITE( PRTFile, * ) '    Perfectly RIGID'
    CASE ( 'A' )
       WRITE( PRTFile, * ) '    ACOUSTO-ELASTIC half-space'
    CASE ( 'G' )
       WRITE( PRTFile, * ) '    Grain size to define half-space'
    CASE ( 'F' )
       WRITE( PRTFile, * ) '    FILE used for reflection loss'
    CASE ( 'W' )
       WRITE( PRTFile, * ) '    Writing an IRC file'
    CASE ( 'P' )
       WRITE( PRTFile, * ) '    reading PRECALCULATED IRC'
    CASE DEFAULT
       CALL ERROUT( 'TopBot', 'Unknown boundary condition type' )
    END SELECT

    ! ****** Read in BC parameters depending on particular choice ******

    HS%cp  = 0.0
    HS%cs  = 0.0
    HS%rho = 0.0

    SELECT CASE ( HS%BC )
    CASE ( 'A' )                  ! *** Half-space properties ***
       zTemp = 0.0
       READ(  ENVFile, *    ) zTemp, alphaR, betaR, rhoR, alphaI, betaI
       WRITE( PRTFile, FMT = "( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
            zTemp, alphaR, betaR, rhoR, alphaI, betaI

       ! dummy parameters for a layer with a general power law for attenuation
       ! these are not in play because the AttenUnit for this is not allowed yet
       !freq0         = freq
       betaPowerLaw  = 1.0
       ft            = 1000.0

       HS%cp  = CRCI( zTemp, alphaR, alphaI, freq, freq, AttenUnit, betaPowerLaw, ft )
       HS%cs  = CRCI( zTemp, betaR,  betaI,  freq, freq, AttenUnit, betaPowerLaw, ft )

       HS%rho = rhoR
    CASE ( 'G' )                  ! *** Grain size (formulas from UW-APL HF Handbook)

       ! These formulas are from the UW-APL Handbook
       ! The code is taken from older Matlab and is unnecesarily verbose
       ! vr   is the sound speed ratio
       ! rhor is the density ratio
       READ(  ENVFile, *    ) zTemp, Mz
       WRITE( PRTFile, FMT = "( F10.2, 3X, F10.2 )" ) zTemp, Mz

       IF ( Mz >= -1 .AND. Mz < 1 ) THEN
          vr   = 0.002709 * Mz ** 2 - 0.056452 * Mz + 1.2778
          rhor = 0.007797 * Mz ** 2 - 0.17057  * Mz + 2.3139
       ELSE IF ( Mz >= 1 .AND. Mz < 5.3 ) THEN
          vr   = -0.0014881 * Mz ** 3 + 0.0213937 * Mz ** 2 - 0.1382798 * Mz + 1.3425
          rhor = -0.0165406 * Mz ** 3 + 0.2290201 * Mz ** 2 - 1.1069031 * Mz + 3.0455
       ELSE
          vr   = -0.0024324 * Mz + 1.0019
          rhor = -0.0012973 * Mz + 1.1565
       END IF

       IF ( Mz >= -1 .AND. Mz < 0 ) THEN
          alpha2_f = 0.4556
       ELSE IF ( Mz >= 0 .AND. Mz < 2.6 ) THEN
          alpha2_f = 0.4556 + 0.0245 * Mz
       ELSE IF( Mz >= 2.6 .AND. Mz < 4.5 ) THEN
          alpha2_f = 0.1978 + 0.1245 * Mz
       ELSE IF( Mz >= 4.5 .AND. Mz < 6.0 ) THEN
          alpha2_f = 8.0399 - 2.5228 * Mz + 0.20098 * Mz ** 2
       ELSE IF( Mz >= 6.0 .AND. Mz < 9.5 ) THEN
          alpha2_f = 0.9431 - 0.2041 * Mz + 0.0117 * Mz ** 2
       ELSE
          alpha2_f =  0.0601
       END IF

       ! AttenUnit = 'L'   ! loss parameter
!!! following uses a reference sound speed of 1500 ???
!!! should be sound speed in the water, just above the sediment
       ! the term vr / 1000 converts vr to units of m per ms 
       alphaR = vr * 1500.0
       alphaI = alpha2_f * ( vr / 1000 ) * 1500.0 * log( 10.0 ) / ( 40.0 * pi )   ! loss parameter Sect. IV., Eq. (4) of handbook

       HS%cp  = CRCI( zTemp, alphaR, alphaI, freq, freq, 'L ', betaPowerLaw, ft )
       HS%cs  = 0.0
       HS%rho = rhoR
       WRITE( PRTFile, FMT = "( 'Converted sound speed =', 2F10.2, 3X, 'density = ', F10.2, 3X, 'loss parm = ', F10.4 )" ) &
            HS%cp, rhor, alphaI

    END SELECT

  END SUBROUTINE TopBot

  ! **********************************************************************!

  SUBROUTINE OpenOutputFiles( FileRoot, ThreeD )
    ! Write appropriate header information

    USE SourceReceiverPositions
    USE angleMod
    USE bdryMod
    USE RWSHDFile

    LOGICAL,            INTENT( IN ) :: ThreeD
    CHARACTER (LEN=80), INTENT( IN ) :: FileRoot
    REAL      (KIND=8)               :: atten
    CHARACTER (LEN=10)               :: PlotType

    SELECT CASE ( Beam%RunType( 1 : 1 ) )
    CASE ( 'R', 'E' )   ! Ray trace or Eigenrays
       OPEN ( FILE = TRIM( FileRoot ) // '.ray', UNIT = RAYFile, FORM = 'FORMATTED' )
       WRITE( RAYFile, * ) '''', Title( 1 : 50 ), ''''
       WRITE( RAYFile, * ) freq
       WRITE( RAYFile, * ) Pos%NSx, Pos%NSy, Pos%NSz
       WRITE( RAYFile, * ) Angles%Nalpha, Angles%Nbeta
       WRITE( RAYFile, * ) Bdry%Top%HS%Depth
       WRITE( RAYFile, * ) Bdry%Bot%HS%Depth

       IF ( ThreeD ) THEN
          WRITE( RAYFile, * ) '''xyz'''
       ELSE
          WRITE( RAYFile, * ) '''rz'''
       END IF

    CASE ( 'A' )        ! arrival file in ascii format
       OPEN ( FILE = TRIM( FileRoot ) // '.arr', UNIT = ARRFile, FORM = 'FORMATTED' )
       ! Note that some compilers will wrap an ascii file at 80 characters
       ! You can add a RECL specifier in the open or use a compiler switch to over-ride

       ! write source locations

       IF ( ThreeD ) THEN
          WRITE( ARRFile, * ) '''3D'''
          WRITE( ARRFile, * ) freq
          WRITE( ARRFile, * ) Pos%NSx,    Pos%Sx(    1 : Pos%NSx )
          WRITE( ARRFile, * ) Pos%NSy,    Pos%Sy(    1 : Pos%NSy )
          WRITE( ARRFile, * ) Pos%NSz,    Pos%Sz(    1 : Pos%NSz )
       ELSE
          WRITE( ARRFile, * ) '''2D'''
          WRITE( ARRFile, * ) freq
          WRITE( ARRFile, * ) Pos%NSz,    Pos%Sz(    1 : Pos%NSz )
       END IF

       ! write receiver locations

       WRITE( ARRFile, *    ) Pos%NRz,    Pos%Rz(    1 : Pos%NRz )
       WRITE( ARRFile, *    ) Pos%NRr,    Pos%Rr(    1 : Pos%NRr )
       IF ( ThreeD ) THEN
          WRITE( ARRFile, * ) Pos%Ntheta, Pos%theta( 1 : Pos%Ntheta )
       END IF

    CASE ( 'a' )        ! arrival file in binary format
       OPEN ( FILE = TRIM( FileRoot ) // '.arr', UNIT = ARRFile, FORM = 'UNFORMATTED' )

       ! write source locations

       IF ( ThreeD ) THEN
          WRITE( ARRFile ) '''3D'''
          WRITE( ARRFile    ) SNGL( freq )
          WRITE( ARRFile    ) Pos%NSx,    Pos%Sx(    1 : Pos%NSx )
          WRITE( ARRFile    ) Pos%NSy,    Pos%Sy(    1 : Pos%NSy )
          WRITE( ARRFile    ) Pos%NSz,    Pos%Sz(    1 : Pos%NSz )
       ELSE
          WRITE( ARRFile ) '''2D'''
          WRITE( ARRFile    ) SNGL( freq )
          WRITE( ARRFile    ) Pos%NSz,    Pos%Sz(    1 : Pos%NSz )
       END IF

       ! write receiver locations

       WRITE( ARRFile       ) Pos%NRz,    Pos%Rz(    1 : Pos%NRz )
       WRITE( ARRFile       ) Pos%NRr,    Pos%Rr(    1 : Pos%NRr )
       IF ( ThreeD ) THEN
          WRITE( ARRFile    ) Pos%Ntheta, Pos%theta( 1 : Pos%Ntheta )
       END IF

    CASE DEFAULT
       atten = 0.0

       ! following to set PlotType has already been done in READIN if that was used for input
       SELECT CASE ( Beam%RunType( 5 : 5 ) )
       CASE ( 'R' )
          PlotType = 'rectilin  '
       CASE ( 'I' )
          PlotType = 'irregular '
       CASE DEFAULT
          PlotType = 'rectilin  '
       END SELECT

       CALL WriteHeader( TRIM( FileRoot ) // '.shd', Title, freq, atten, PlotType )
    END SELECT

  END SUBROUTINE OpenOutputFiles

END MODULE ReadEnvironmentBell
