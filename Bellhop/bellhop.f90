PROGRAM BELLHOP

  ! BELLHOP Beam tracing for ocean acoustics

  ! Copyright (C) 2009 Michael B. Porter

  ! This program is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.

  ! You should have received a copy of the GNU General Public License
  ! along with this program.  If not, see <http://www.gnu.org/licenses/>.

  ! First version (1983) originally developed with Homer Bucker, Naval Ocean Systems Center
  
  USE ReadEnvironmentBell
  USE BeamPattern
  USE bdryMod
  USE RefCoef
  USE SourceReceiverPositions
  USE angleMod
  USE sspMod
  USE influence
  USE FatalError

  IMPLICIT NONE
  
  LOGICAL, PARAMETER   :: Inline = .FALSE.
  INTEGER              :: jj
  CHARACTER ( LEN=2  ) :: AttenUnit
  CHARACTER ( LEN=80 ) :: FileRoot

  ThreeD = .FALSE.

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )
  ! Open the print file
  OPEN( UNIT = PRTFile, FILE = TRIM( FileRoot ) // '.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )

  ! Read in or otherwise initialize inline all the variables used by BELLHOP 

  IF ( Inline ) THEN
     ! NPts, Sigma not used by BELLHOP
     Title = 'BELLHOP- Calibration case with envfil passed as parameters'
     freq  = 250
     ! NMedia variable is not used by BELLHOP

     ! *** Boundary information (type of boundary condition and, if a halfspace, then halfspace info)

     AttenUnit         = 'W'
     Bdry%Top%HS%BC    = 'V'
     Bdry%Top%HS%Depth = 0
     Bdry%Bot%HS%Depth = 100
     Bdry%Bot%HS%Opt   = 'A_'
     Bdry%Bot%HS%BC    = 'A'
     Bdry%Bot%HS%cp    = CRCI( 1D20, 1590D0,    0.5D0, freq, freq, AttenUnit, betaPowerLaw, fT )   ! compressional wave speed
     Bdry%Bot%HS%cs    = CRCI( 1D20,    0D0,      0D0, freq, freq, AttenUnit, betaPowerLaw, fT )   ! shear         wave speed
     Bdry%Bot%HS%rho   = 1.2                               ! density

     ! *** sound speed in the water column ***

     SSP%Type = 'C'   ! interpolation method for SSP
     SSP%NPts = 2     ! number of SSP points
     SSP%z(  1 : 2 ) = [    0,  100 ]
     SSP%c(  1 : 2 ) = [ 1500, 1500 ]
     SSP%cz( 1 : 2 ) = [    0,    0 ]   ! user should really not have to supply this ...

     ! *** source and receiver positions ***

     Pos%NSz = 1
     Pos%NRz = 100
     Pos%NRr  = 500

     ALLOCATE( Pos%sz( Pos%NSz ), Pos%ws( Pos%NSz ), Pos%isz( Pos%NSz ) )
     ALLOCATE( Pos%rz( Pos%NRz ), Pos%wr( Pos%NRz ), Pos%irz( Pos%NRz ) )
     ALLOCATE( Pos%Rr(  Pos%NRr  ) )

     Pos%Sz( 1 ) = 50.
     !Pos%Rz     = [ 0, 50, 100 ]
     !Pos%Rr     = 1000. * [ 1, 2, 3, 4, 5 ]   ! meters !!!
     Pos%Rz      = [ ( jj, jj = 1, Pos%NRz ) ]
     Pos%Rr      = 10. * [ ( jj, jj = 1 , Pos%NRr ) ]   ! meters !!!

     Beam%RunType = 'C'
     Beam%Type    = 'G   '
     Beam%deltas  = 0
     Beam%Box%z   = 101.
     Beam%Box%r   = 5100 ! meters

     Angles%Nalpha = 1789
     ! Angles%alpha  = [ -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80 ] ! -89 89
     Angles%alpha  = ( 180. / Angles%Nalpha ) * [ ( jj, jj = 1, Angles%Nalpha ) ] - 90.

     ! *** altimetry ***

     ALLOCATE( Top( 2 ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP', 'Insufficient memory for altimetry data'  )
     Top( 1 )%x = [ -sqrt( huge( Top( 1 )%x( 1 ) ) ) / 1.0d5, 0.d0 ]
     Top( 2 )%x = [  sqrt( huge( Top( 1 )%x( 1 ) ) ) / 1.0d5, 0.d0 ]

     CALL ComputeBdryTangentNormal( Top, 'Top' )

     ! *** bathymetry ***

     ALLOCATE( Bot( 2 ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP', 'Insufficient memory for bathymetry data'  )
     Bot( 1 )%x = [ -sqrt( huge( Bot( 1 )%x( 1 ) ) ) / 1.0d5, 5000.d0 ]
     Bot( 2 )%x = [  sqrt( huge( Bot( 1 )%x( 1 ) ) ) / 1.0d5, 5000.d0 ]

     CALL ComputeBdryTangentNormal( Bot, 'Bot' )

     ALLOCATE( RBot( 1 ), Stat = IAllocStat )   ! bottom reflection coefficient
     ALLOCATE( RTop( 1 ), Stat = iAllocStat )   ! top    reflection coefficient

     NBotPts = 1
     NTopPts = 1

     ! *** Source Beam Pattern ***
     NSBPPts = 2
     ALLOCATE( SrcBmPat( 2, 2 ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP-ReadPat', 'Insufficient memory'  )
     SrcBmPat( 1, : ) = [ -180.0, 0.0 ]
     SrcBmPat( 2, : ) = [  180.0, 0.0 ]
     SrcBmPat( :, 2 ) = 10 ** ( SrcBmPat( :, 2 ) / 20 )  ! convert dB to linear scale !!!
  ELSE
     CALL ReadEnvironment(  FileRoot, ThreeD )
     CALL ReadATI( FileRoot, Bdry%Top%HS%Opt( 5 : 5 ), Bdry%Top%HS%Depth, PRTFile )                          ! AlTImetry
     CALL ReadBTY( FileRoot, Bdry%Bot%HS%Opt( 2 : 2 ), Bdry%Bot%HS%Depth, PRTFile )                          ! BaThYmetry
     CALL ReadReflectionCoefficient( FileRoot, Bdry%Bot%HS%Opt( 1 : 1 ), Bdry%Top%HS%Opt( 2 : 2 ), PRTFile ) ! (top and bottom)
     SBPFlag = Beam%RunType( 3 : 3 )
     CALL ReadPat( FileRoot, PRTFile )   ! Source Beam Pattern
     ! dummy bearing angles
     Pos%Ntheta = 1
     ALLOCATE( Pos%theta( Pos%Ntheta ), Stat = IAllocStat )
     Pos%theta( 1 ) = 0.
  END IF

  CALL OpenOutputFiles( FileRoot, ThreeD )
  CALL BellhopCore

  CONTAINS

! **********************************************************************!

SUBROUTINE BellhopCore

  USE ArrMod
  USE AttenMod

  INTEGER, PARAMETER   :: ArrivalsStorage = 20000000, MinNArr = 10
  INTEGER              :: IBPvec( 1 ), ibp, is, iBeamWindow2, Irz1, Irec, NalphaOpt, iSeg
  REAL                 :: Tstart, Tstop
  REAL        (KIND=8) :: Amp0, DalphaOpt, xs( 2 ), RadMax, s, &
                          c, cimag, gradc( 2 ), crr, crz, czz, rho
  COMPLEX, ALLOCATABLE :: U( :, : )
  COMPLEX     (KIND=8) :: epsilon

  CALL CPU_TIME( Tstart )

  omega = 2.0 * pi * freq

  IF ( Beam%deltas == 0.0 ) THEN
     Beam%deltas = ( Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth ) / 10.0   ! Automatic step size selection
     WRITE( PRTFile, * )
     WRITE( PRTFile, fmt = '(  '' Step length,       deltas = '', G11.4, '' m (automatically selected)'' )' ) Beam%deltas
  END IF

  Angles%alpha  = DegRad * Angles%alpha   ! convert to radians
  Angles%Dalpha = 0.0
  IF ( Angles%Nalpha /= 1 ) &
       Angles%Dalpha = ( Angles%alpha( Angles%Nalpha ) - Angles%alpha( 1 ) ) / ( Angles%Nalpha - 1 )  ! angular spacing between beams

  ! convert range-dependent geoacoustic parameters from user to program units
  IF ( atiType( 2 : 2 ) == 'L' ) THEN
     DO iSeg = 1, NatiPts
        Top( iSeg )%HS%cp = CRCI( 1D20, Top( iSeg )%HS%alphaR, Top( iSeg )%HS%alphaI, freq, freq, 'W ', &
             betaPowerLaw, ft )   ! compressional wave speed
        Top( iSeg )%HS%cs = CRCI( 1D20, Top( iSeg )%HS%betaR,  Top( iSeg )%HS%betaI,  freq, freq, 'W ', &
             betaPowerLaw, ft )   ! shear         wave speed
     END DO
  END IF
   
  IF ( btyType( 2 : 2 ) == 'L' ) THEN
     DO iSeg = 1, NbtyPts
        Bot( iSeg )%HS%cp = CRCI( 1D20, Bot( iSeg )%HS%alphaR, Bot( iSeg )%HS%alphaI, freq, freq, 'W ', &
             betaPowerLaw, ft )   ! compressional wave speed 
        Bot( iSeg )%HS%cs = CRCI( 1D20, Bot( iSeg )%HS%betaR,  Bot( iSeg )%HS%betaI,  freq, freq, 'W ', &
             betaPowerLaw, ft )   ! shear         wave speed
     END DO
  END IF

  SELECT CASE ( Beam%RunType( 5 : 5 ) )
  CASE ( 'I' )
     NRz_per_range = 1         ! irregular grid
  CASE DEFAULT
     NRz_per_range = Pos%NRz   ! rectilinear grid
  END SELECT

  ! for a TL calculation, allocate space for the pressure matrix
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'C', 'S', 'I' )        ! TL calculation
     ALLOCATE ( U( NRz_per_range, Pos%NRr ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP', 'Insufficient memory for TL matrix: reduce Nr * NRz'  )
  CASE ( 'A', 'a', 'R', 'E' )   ! Arrivals calculation
     ALLOCATE ( U( 1, 1 ), Stat = IAllocStat )   ! open a dummy variable
  END SELECT

  ! for an arrivals run, allocate space for arrivals matrices
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'A', 'a' )
     MaxNArr = MAX( ArrivalsStorage / ( NRz_per_range * Pos%NRr ), MinNArr )   ! allow space for at least 10 arrivals
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) '( Maximum # of arrivals = ', MaxNArr, ')'

     ALLOCATE ( Arr( NRz_per_range, Pos%NRr, MaxNArr ), NArr( NRz_per_range, Pos%NRr ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP', &
          'Insufficient memory to allocate arrivals matrix; reduce parameter ArrivalsStorage' )
  CASE DEFAULT
     MaxNArr = 1
     ALLOCATE ( Arr( NRz_per_range, Pos%NRr, 1 ), NArr( NRz_per_range, Pos%NRr ), Stat = IAllocStat )
  END SELECT

  WRITE( PRTFile, * )

  SourceDepth: DO is = 1, Pos%NSz
     xs = [ 0.0, Pos%sz( is ) ]   ! source coordinate

     SELECT CASE ( Beam%RunType( 1 : 1 ) )
     CASE ( 'C', 'S', 'I' ) ! TL calculation, zero out pressure matrix
        U = 0.0
     CASE ( 'A', 'a' )      ! Arrivals calculation, zero out arrival matrix
        NArr = 0
     END SELECT

     iSegr = 1   ! this is not really necessary, but ensures the segment search process begins in the same segment for each ray
     iSegz = 1

     CALL EvaluateSSP( xs, c, cimag, gradc, crr, crz, czz, rho, freq, 'TAB' )
     RadMax = 5 * c / freq  ! 5 wavelength max radius

     ! Are there enough beams?
     DalphaOpt = SQRT( c / ( 6.0 * freq * Pos%Rr( Pos%NRr ) ) )
     NalphaOpt = 2 + INT( ( Angles%alpha( Angles%Nalpha ) - Angles%alpha( 1 ) ) / DalphaOpt )

     IF ( Beam%RunType( 1 : 1 ) == 'C' .AND. Angles%Nalpha < NalphaOpt ) THEN
        WRITE( PRTFile, * ) 'Warning in BELLHOP : Too few beams'
        WRITE( PRTFile, * ) 'Nalpha should be at least = ', NalphaOpt
     ENDIF

     ! Trace successive beams

     DeclinationAngle: DO ialpha = 1, Angles%Nalpha
        SrcDeclAngle = RadDeg * Angles%alpha( ialpha )          ! take-off declination angle in degrees

        IF ( Angles%iSingle_alpha == 0 .OR. ialpha == Angles%iSingle_alpha ) THEN    ! Single beam run?

           IBPvec = maxloc( SrcBmPat( :, 1 ), mask = SrcBmPat( :, 1 ) < SrcDeclAngle )  ! index of ray angle in beam pattern
           IBP    = IBPvec( 1 )
           IBP    = MAX( IBP, 1 )               ! don't go before beginning of table
           IBP    = MIN( IBP, NSBPPts - 1 )     ! don't go past end of table

           ! linear interpolation to get amplitude
           s    = ( SrcDeclAngle  - SrcBmPat( IBP, 1 ) ) / ( SrcBmPat( IBP + 1, 1 ) - SrcBmPat( IBP, 1 ) )
           Amp0 = ( 1 - s ) * SrcBmPat( IBP, 2 ) + s * SrcBmPat( IBP + 1, 2 )

           ! Lloyd mirror pattern for semi-coherent option
           IF ( Beam%RunType( 1 : 1 ) == 'S' ) &
              Amp0 = Amp0 * SQRT( 2.0 ) * ABS( SIN( omega / c * xs( 2 ) * SIN( Angles%alpha( ialpha ) ) ) )

           ! show progress ...
           IF ( MOD( ialpha - 1, max( Angles%Nalpha / 50, 1 ) ) == 0 ) THEN
              WRITE( PRTFile, FMT = "( 'Tracing beam ', I7, F10.2 )" ) ialpha, SrcDeclAngle
              FLUSH( PRTFile )
           END IF

           CALL TraceRay2D( xs, Angles%alpha( ialpha ), Amp0 )   ! Trace a ray

           IF ( Beam%RunType( 1 : 1 ) == 'R' ) THEN   ! Write the ray trajectory to RAYFile
              CALL WriteRay2D( SrcDeclAngle, Beam%Nsteps )
           ELSE                                       ! Compute the contribution to the field

              epsilon = PickEpsilon( Beam%Type( 1 : 2 ), omega, c, gradc, Angles%alpha( ialpha ), &
                   Angles%Dalpha, Beam%rLoop, Beam%epsMultiplier ) ! 'optimal' beam constant

              SELECT CASE ( Beam%Type( 1 : 1 ) )
              CASE ( 'R' )
                 iBeamWindow2 = Beam%iBeamWindow **2
                 RadMax       = 50 * c / freq  ! 50 wavelength max radius
                 CALL InfluenceCervenyRayCen(   U, epsilon, Angles%alpha( ialpha ), iBeamWindow2, RadMax )
              CASE ( 'C' )
                 iBeamWindow2 = Beam%iBeamWindow **2
                 RadMax       = 50 * c / freq  ! 50 wavelength max radius
                 CALL InfluenceCervenyCart(     U, epsilon, Angles%alpha( ialpha ), iBeamWindow2, RadMax )
              CASE ( 'g' )
                 CALL InfluenceGeoHatRayCen(    U, Angles%alpha( ialpha ), Angles%Dalpha )
              CASE ( 'S' )
                 CALL InfluenceSGB(             U, Angles%alpha( ialpha ), Angles%Dalpha )
              CASE ( 'B' )
                 CALL InfluenceGeoGaussianCart( U, Angles%alpha( ialpha ), Angles%Dalpha )
             CASE DEFAULT
                 CALL InfluenceGeoHatCart(      U, Angles%alpha( ialpha ), Angles%Dalpha )
              END SELECT

           END IF
        END IF
     END DO DeclinationAngle

     ! write results to disk

     SELECT CASE ( Beam%RunType( 1 : 1 ) )
     CASE ( 'C', 'S', 'I' )   ! TL calculation
        CALL ScalePressure( Angles%Dalpha, ray2D( 1 )%c, Pos%Rr, U, NRz_per_range, Pos%NRr, Beam%RunType, freq )
        IRec = 10 + NRz_per_range * ( is - 1 )
        RcvrDepth: DO Irz1 = 1, NRz_per_range
           IRec = IRec + 1
           WRITE( SHDFile, REC = IRec ) U( Irz1, 1 : Pos%NRr )
        END DO RcvrDepth
     CASE ( 'A' )             ! arrivals calculation, ascii
        CALL WriteArrivalsASCII(  Pos%Rr, NRz_per_range, Pos%NRr, Beam%RunType( 4 : 4 ) )
     CASE ( 'a' )             ! arrivals calculation, binary
        CALL WriteArrivalsBinary( Pos%Rr, NRz_per_range, Pos%NRr, Beam%RunType( 4 : 4 ) )
     END SELECT

  END DO SourceDepth

  ! Display run time
  CALL CPU_TIME( Tstop )
  WRITE( PRTFile, "( /, ' CPU Time = ', G15.3, 's' )" ) Tstop - Tstart

  ! close all files
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'C', 'S', 'I' )      ! TL calculation
     CLOSE( SHDFile )
  CASE ( 'A', 'a' )           ! arrivals calculation
     CLOSE( ARRFile )
  CASE ( 'R', 'E' )           ! ray trace
     CLOSE( RAYFile )
  END SELECT

  CLOSE( PRTFile )

END SUBROUTINE BellhopCore

! **********************************************************************!

COMPLEX (KIND=8 ) FUNCTION PickEpsilon( BeamType, omega, c, gradc, alpha, Dalpha, rLoop, EpsMultiplier )

  ! Picks the optimum value for epsilon

  REAL      (KIND=8), INTENT( IN  ) :: omega, c, gradc( 2 ) ! angular frequency, sound speed and gradient
  REAL      (KIND=8), INTENT( IN  ) :: alpha, Dalpha        ! angular spacing for ray fan
  REAL      (KIND=8), INTENT( IN  ) :: epsMultiplier, Rloop ! multiplier, loop range
  CHARACTER (LEN= 2), INTENT( IN  ) :: BeamType
  LOGICAL, SAVE      :: INIFlag = .TRUE.
  REAL      (KIND=8) :: HalfWidth = 0.0
  REAL      (KIND=8) :: cz
  COMPLEX   (KIND=8) :: epsilonOpt = 0.0
  CHARACTER (LEN=40) :: TAG

  SELECT CASE ( BeamType( 1 : 1 ) )
  CASE ( 'C', 'R' )
     TAG    = 'Paraxial beams'
     SELECT CASE ( BeamType( 2 : 2 ) )
     CASE ( 'F' )
        TAG        = 'Space filling beams'
        halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
        epsilonOpt = i * 0.5 * omega * halfwidth ** 2
     CASE ( 'M' )
        TAG        = 'Minimum width beams'
        halfwidth  = SQRT( 2.0 * c * 1000.0 * rLoop / omega )
        epsilonOpt = i * 0.5 * omega * halfwidth ** 2
     CASE ( 'W' )
        TAG       = 'WKB beams'
        halfwidth = HUGE( halfwidth )
        cz        = gradc( 2 )
        IF ( cz == 0.0 ) THEN
           epsilonOpt = 1.0D10
        ELSE
           epsilonOpt = ( -SIN( alpha ) / COS( alpha ** 2 ) ) * c * c / cz
        ENDIF
     END SELECT

  CASE ( 'G', 'g', '^' )
     TAG        = 'Geometric hat beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2

  CASE ( 'B' )
     TAG        = 'Geometric Gaussian beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2

  CASE ( 'b' )
     CALL ERROUT( 'BELLHOP', 'Geo Gaussian beams in ray-centered coords. not implemented in BELLHOP' )

  CASE ( 'S' )
     TAG        = 'Simple Gaussian beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2
  END SELECT

  PickEpsilon = EpsMultiplier * epsilonOpt

  ! On first call write info to prt file
  IF ( INIFlag ) THEN
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) TAG
     WRITE( PRTFile, * ) 'halfwidth  = ', halfwidth
     WRITE( PRTFile, * ) 'epsilonOpt = ', epsilonOpt
     WRITE( PRTFile, * ) 'EpsMult    = ', EpsMultiplier
     WRITE( PRTFile, * )
     INIFlag = .FALSE.
  END IF

END FUNCTION PickEpsilon

! **********************************************************************!

SUBROUTINE TraceRay2D( xs, alpha, Amp0 )

  ! Traces the beam corresponding to a particular take-off angle

  USE Step
  USE WriteRay

  REAL     (KIND=8), INTENT( IN ) :: xs( 2 )      ! x-y coordinate of the source
  REAL     (KIND=8), INTENT( IN ) :: alpha, Amp0  ! initial angle, amplitude
  INTEGER           :: is, is1                    ! index for a step along the ray
  REAL     (KIND=8) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
  REAL     (KIND=8) :: dEndTop( 2 ), dEndBot( 2 ), TopnInt( 2 ), BotnInt( 2 ), ToptInt( 2 ), BottInt( 2 )
  REAL     (KIND=8) :: DistBegTop, DistEndTop, DistBegBot, DistEndBot ! Distances from ray beginning, end to top and bottom
  REAL     (KIND=8) :: sss

  ! Initial conditions

  iSegr = 1   ! this is not really necessary, but ensures the segment search process begins in the same segment for each ray
  iSegz = 1
  
  iSmallStepCtr = 0
  CALL EvaluateSSP( xs, c, cimag, gradc, crr, crz, czz, rho, freq, 'TAB' )
  ray2D( 1 )%c         = c
  ray2D( 1 )%x         = xs
  ray2D( 1 )%t         = [ COS( alpha ), SIN( alpha ) ] / c
  ray2D( 1 )%p         = [ 1.0, 0.0 ]
  ray2D( 1 )%q         = [ 0.0, 1.0 ]
  ray2D( 1 )%tau       = 0.0
  ray2D( 1 )%Amp       = Amp0
  ray2D( 1 )%Phase     = 0.0
  ray2D( 1 )%NumTopBnc = 0
  ray2D( 1 )%NumBotBnc = 0

  ! second component of qv is not used in geometric beam tracing
  ! set I.C. to 0 in hopes of saving run time
  IF ( Beam%RunType( 2 : 2 ) == 'G' ) ray2D( 1 )%q = [ 0.0, 0.0 ]

  CALL GetTopSeg( xs( 1 ) )   ! identify the top    segment above the source
  CALL GetBotSeg( xs( 1 ) )   ! identify the bottom segment below the source

  ! convert range-dependent geoacoustic parameters from user to program units
  ! compiler is not accepting the copy of the whole structure at once ...
  IF ( atiType( 2 : 2 ) == 'L' ) THEN
     Bdry%Top%HS%cp  = Top( IsegTop )%HS%cp   ! grab the geoacoustic info for the new segment
     Bdry%Top%HS%cs  = Top( IsegTop )%HS%cs
     Bdry%Top%HS%rho = Top( IsegTop )%HS%rho
  END IF

  IF ( btyType( 2 : 2 ) == 'L' ) THEN
     Bdry%Bot%HS%cp  = Bot( IsegBot )%HS%cp
     Bdry%Bot%HS%cs  = Bot( IsegBot )%HS%cs
     Bdry%Bot%HS%rho = Bot( IsegBot )%HS%rho
  END IF

  ! Trace the beam (note that Reflect alters the step index, is)
  is = 0
  CALL Distances2D( ray2D( 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
       Top( IsegTop )%n, Bot( IsegBot )%n, DistBegTop, DistBegBot )

  IF ( DistBegTop <= 0 .OR. DistBegBot <= 0 ) THEN
     Beam%Nsteps = 1
     WRITE( PRTFile, * ) 'Terminating the ray trace because the source is on or outside the boundaries'
     RETURN       ! source must be within the medium
  END IF

  Stepping: DO istep = 1, MaxN - 1
     is  = is + 1
     is1 = is + 1

     CALL Step2D( ray2D( is ), ray2D( is1 ),  &
          Top( IsegTop )%x, Top( IsegTop )%n, &
          Bot( IsegBot )%x, Bot( IsegBot )%n )

     ! New altimetry segment?
     IF ( ray2D( is1 )%x( 1 ) < rTopSeg( 1 ) .OR. &
          ray2D( is1 )%x( 1 ) > rTopSeg( 2 ) ) THEN
        CALL GetTopSeg( ray2D( is1 )%x( 1 ) )
        IF ( atiType( 2 : 2 ) == 'L' ) THEN
           Bdry%Top%HS%cp  = Top( IsegTop )%HS%cp   ! grab the geoacoustic info for the new segment
           Bdry%Top%HS%cs  = Top( IsegTop )%HS%cs
           Bdry%Top%HS%rho = Top( IsegTop )%HS%rho
        END IF
     END IF

     ! New bathymetry segment?
     IF ( ray2D( is1 )%x( 1 ) < rBotSeg( 1 ) .OR. &
          ray2D( is1 )%x( 1 ) > rBotSeg( 2 ) ) THEN
        CALL GetBotSeg( ray2D( is1 )%x( 1 ) )
        IF ( btyType( 2 : 2 ) == 'L' ) THEN
           Bdry%Bot%HS%cp  = Bot( IsegBot )%HS%cp   ! grab the geoacoustic info for the new segment
           Bdry%Bot%HS%cs  = Bot( IsegBot )%HS%cs
           Bdry%Bot%HS%rho = Bot( IsegBot )%HS%rho
        END IF
     END IF

     ! Reflections?
     ! Tests that ray at step is is inside, and ray at step is+1 is outside
     ! to detect only a crossing from inside to outside
     ! DistBeg is the distance at step is,   which is saved
     ! DistEnd is the distance at step is+1, which needs to be calculated

     CALL Distances2D( ray2D( is1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
          Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     IF      ( DistBegTop > 0.0d0 .AND. DistEndTop <= 0.0d0 ) THEN  ! test top reflection

        IF ( atiType == 'C' ) THEN
           sss     = DOT_PRODUCT( dEndTop, Top( IsegTop )%t ) / Top( IsegTop )%Len   ! proportional distance along segment
           TopnInt = ( 1 - sss ) * Top( IsegTop )%Noden + sss * Top( 1 + IsegTop )%Noden
           ToptInt = ( 1 - sss ) * Top( IsegTop )%Nodet + sss * Top( 1 + IsegTop )%Nodet
        ELSE
           TopnInt = Top( IsegTop )%n   ! normal is constant in a segment
           ToptInt = Top( IsegTop )%t
        END IF

        CALL Reflect2D( is, Bdry%Top%HS, 'TOP', ToptInt, TopnInt, Top( IsegTop )%kappa, RTop, NTopPTS )
        ray2D( is + 1 )%NumTopBnc = ray2D( is )%NumTopBnc + 1

        CALL Distances2D( ray2D( is + 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
             Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     ELSE IF ( DistBegBot > 0.0d0 .AND. DistEndBot <= 0.0d0 ) THEN  ! test bottom reflection

        IF ( btyType == 'C' ) THEN
           sss     = DOT_PRODUCT( dEndBot, Bot( IsegBot )%t ) / Bot( IsegBot )%Len   ! proportional distance along segment
           BotnInt = ( 1 - sss ) * Bot( IsegBot )%Noden + sss * Bot( 1 + IsegBot )%Noden
           BottInt = ( 1 - sss ) * Bot( IsegBot )%Nodet + sss * Bot( 1 + IsegBot )%Nodet
        ELSE
           BotnInt = Bot( IsegBot )%n   ! normal is constant in a segment
           BottInt = Bot( IsegBot )%t
        END IF

        CALL Reflect2D( is, Bdry%Bot%HS, 'BOT', BottInt, BotnInt, Bot( IsegBot )%kappa, RBot, NBotPTS )
        ray2D( is + 1 )%NumBotBnc = ray2D( is )%NumBotBnc + 1
        CALL Distances2D( ray2D( is + 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
             Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     END IF

     ! Has the ray left the box, lost its energy, escaped the boundaries, or exceeded storage limit?
     ! There is a test here for when two adjacent ray points are outside the boundaries
     ! BELLHOP3D handles this differently using a counter-limit for really small steps.
     IF ( ABS( ray2D( is + 1 )%x( 1 ) ) > Beam%Box%r .OR. &
          ABS( ray2D( is + 1 )%x( 2 ) ) > Beam%Box%z .OR. ray2D( is + 1 )%Amp < 0.005 .OR. &
          ( DistBegTop < 0.0 .AND. DistEndTop < 0.0 ) .OR. &
          ( DistBegBot < 0.0 .AND. DistEndBot < 0.0 ) ) THEN
          ! ray2D( is + 1 )%t( 1 ) < 0 ) THEN ! this last test kills off a backward traveling ray
        Beam%Nsteps = is + 1
        EXIT Stepping
     ELSE IF ( is >= MaxN - 3 ) THEN
        WRITE( PRTFile, * ) 'Warning in TraceRay2D : Insufficient storage for ray trajectory'
        Beam%Nsteps = is
        EXIT Stepping
     END IF

     DistBegTop = DistEndTop
     DistBegBot = DistEndBot

  END DO Stepping
END SUBROUTINE TraceRay2D

! **********************************************************************!

SUBROUTINE Distances2D( rayx, Topx, Botx, dTop, dBot, Topn, Botn, DistTop, DistBot )

  ! Calculates the distances to the boundaries
  ! Formula differs from JKPS because code uses outward pointing normals

  REAL (KIND=8), INTENT( IN  ) :: rayx( 2 )              ! ray coordinate
  REAL (KIND=8), INTENT( IN  ) :: Topx( 2 ), Botx( 2 )   ! top, bottom coordinate
  REAL (KIND=8), INTENT( IN  ) :: Topn( 2 ), Botn( 2 )   ! top, bottom normal vector (outward)
  REAL (KIND=8), INTENT( OUT ) :: dTop( 2 ), dBot( 2 )   ! vector pointing from top, bottom bdry to ray
  REAL (KIND=8), INTENT( OUT ) :: DistTop, DistBot       ! distance (normal to bdry) from the ray to top, bottom boundary

  dTop    = rayx - Topx  ! vector pointing from top    to ray
  dBot    = rayx - Botx  ! vector pointing from bottom to ray
  DistTop = -DOT_PRODUCT( Topn, dTop )
  DistBot = -DOT_PRODUCT( Botn, dBot )

END SUBROUTINE Distances2D

! **********************************************************************!

SUBROUTINE Reflect2D( is, HS, BotTop, tBdry, nBdry, kappa, RefC, Npts )

  INTEGER,              INTENT( IN ) :: Npts
  REAL     (KIND=8),    INTENT( IN ) :: tBdry( 2 ), nBdry( 2 )  ! Tangent and normal to the boundary
  REAL     (KIND=8),    INTENT( IN ) :: kappa                   ! Boundary curvature
  CHARACTER (LEN=3),    INTENT( IN ) :: BotTop                  ! Flag indicating bottom or top reflection
  TYPE( HSInfo ),       INTENT( IN ) :: HS                      ! half-space properties
  TYPE(ReflectionCoef), INTENT( IN ) :: RefC( NPts )            ! reflection coefficient
  INTEGER,              INTENT( INOUT ) :: is
  INTEGER           :: is1
  REAL     (KIND=8) :: c, cimag, gradc( 2 ), crr, crz, czz, rho                  ! derivatives of sound speed
  REAL     (KIND=8) :: RM, RN, Tg, Th, rayt( 2 ), rayn( 2 ), rayt_tilde( 2 ), rayn_tilde( 2 ), cnjump, csjump  ! for curvature change
  REAL     (KIND=8) :: ck, co, si, cco, ssi, pdelta, rddelta, sddelta, theta_bot ! for beam shift
  COMPLEX  (KIND=8) :: kx, kz, kzP, kzS, kzP2, kzS2, mu, f, g, y2, y4, Refl      ! for tabulated reflection coef.
  COMPLEX  (KIND=8) :: ch, a, b, d, sb, delta, ddelta                            ! for beam shift
  TYPE(ReflectionCoef) :: RInt

  is  = is + 1
  is1 = is + 1

  Tg = DOT_PRODUCT( ray2D( is )%t, tBdry )  ! component of ray tangent, along boundary
  Th = DOT_PRODUCT( ray2D( is )%t, nBdry )  ! component of ray tangent, normal to boundary

  ray2D( is1 )%NumTopBnc = ray2D( is )%NumTopBnc
  ray2D( is1 )%NumBotBnc = ray2D( is )%NumBotBnc
  ray2D( is1 )%x         = ray2D( is )%x
  ray2D( is1 )%t         = ray2D( is )%t - 2.0 * Th * nBdry  ! changing the ray direction

  ! Calculate the change in curvature
  ! Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).

  CALL EvaluateSSP( ray2D( is )%x, c, cimag, gradc, crr, crz, czz, rho, freq, 'TAB' )   ! just to get c

  ! incident unit ray tangent and normal
  rayt = c * ray2D( is )%t                              ! unit tangent to ray
  rayn = [ -rayt( 2 ), rayt( 1 ) ]                      ! unit normal  to ray

  ! reflected unit ray tangent and normal (the reflected tangent, normal system has a different orientation)
  rayt_tilde = c * ray2D( is1 )%t                       ! unit tangent to ray
  rayn_tilde = -[ -rayt_tilde( 2 ), rayt_tilde( 1 ) ]   ! unit normal  to ray

  RN = 2 * kappa / c ** 2 / Th    ! boundary curvature correction

  ! get the jumps (this could be simplified, e.g. jump in rayt is roughly 2 * Th * nbdry
  cnjump = -DOT_PRODUCT( gradc, rayn_tilde - rayn )
  csjump = -DOT_PRODUCT( gradc, rayt_tilde - rayt )

  IF ( BotTop == 'TOP' ) THEN
     cnjump = -cnjump   ! this is because the (t,n) system of the top boundary has a different sense to the bottom boundary
     RN     = -RN
  END IF

  RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
  RN = RN + RM * ( 2 * cnjump - RM * csjump ) / c ** 2

  SELECT CASE ( Beam%Type( 3 : 3 ) )
  CASE ( 'D' )
     RN = 2.0 * RN
  CASE ( 'Z' )
     RN = 0.0
  END SELECT
  
  ray2D( is1 )%c   = c
  ray2D( is1 )%tau = ray2D( is )%tau
  ray2D( is1 )%p   = ray2D( is )%p + ray2D( is )%q * RN
  ray2D( is1 )%q   = ray2D( is )%q

  ! account for phase change

  SELECT CASE ( HS%BC )
  CASE ( 'R' )                 ! rigid
     ray2D( is1 )%Amp   = ray2D( is )%Amp
     ray2D( is1 )%Phase = ray2D( is )%Phase
  CASE ( 'V' )                 ! vacuum
     ray2D( is1 )%Amp   = ray2D( is )%Amp
     ray2D( is1 )%Phase = ray2D( is )%Phase + pi
  CASE ( 'F' )                 ! file
     RInt%theta = RadDeg * ABS( ATAN2( Th, Tg ) )           ! angle of incidence (relative to normal to bathymetry)
     IF ( RInt%theta > 90 ) RInt%theta = 180. - RInt%theta  ! reflection coefficient is symmetric about 90 degrees
     CALL InterpolateReflectionCoefficient( RInt, RefC, Npts, PRTFile )
     ray2D( is1 )%Amp   = ray2D( is )%Amp * RInt%R
     ray2D( is1 )%Phase = ray2D( is )%Phase + RInt%phi
  CASE ( 'A', 'G' )            ! half-space
     kx = omega * Tg     ! wavenumber in direction parallel      to bathymetry
     kz = omega * Th     ! wavenumber in direction perpendicular to bathymetry (in ocean)

     ! notation below is a bit misleading
     ! kzS, kzP is really what I called gamma in other codes, and differs by a factor of +/- i
     IF ( REAL( HS%cS ) > 0.0 ) THEN
        kzS2 = kx ** 2 - ( omega / HS%cS ) ** 2
        kzP2 = kx ** 2 - ( omega / HS%cP ) ** 2
        kzS  = SQRT( kzS2 )
        kzP  = SQRT( kzP2 )
        mu   = HS%rho * HS%cS ** 2

        y2 = ( ( kzS2 + kx ** 2 ) ** 2 - 4.0D0 * kzS * kzP * kx ** 2 ) * mu
        y4 = kzP * ( kx ** 2 - kzS2 )

        f = omega ** 2 * y4
        g = y2
     ELSE
        kzP = SQRT( kx ** 2 - ( omega / HS%cP ) ** 2 )

        ! Intel and GFortran compilers return different branches of the SQRT for negative reals
        IF ( REAL( kzP ) == 0.0D0 .AND. AIMAG( kzP ) < 0.0D0 ) kzP = -kzP
        f   = kzP
        g   = HS%rho
     ENDIF

     Refl =  - ( rho * f - i * kz * g ) / ( rho * f + i * kz * g )   ! complex reflection coef.

     IF ( ABS( Refl ) < 1.0E-5 ) THEN   ! kill a ray that has lost its energy in reflection
        ray2D( is1 )%Amp   = 0.0
        ray2D( is1 )%Phase = ray2D( is )%Phase
     ELSE
        ray2D( is1 )%Amp   = ABS( Refl ) * ray2D(  is )%Amp
        ray2D( is1 )%Phase = ray2D( is )%Phase + ATAN2( AIMAG( Refl ), REAL( Refl ) )

        ! compute beam-displacement Tindle, Eq. (14)
        ! needs a correction to beam-width as well ...
        !  IF ( REAL( kz2Sq ) < 0.0 ) THEN
        !     rhoW   = 1.0   ! density of water
        !     rhoWSq  = rhoW  * rhoW
        !     rhoHSSq = rhoHS * rhoHS
        !     DELTA = 2 * GK * rhoW * rhoHS * ( kz1Sq - kz2Sq ) /
        ! &( kz1 * i * kz2 *
        ! &( -rhoWSq * kz2Sq + rhoHSSq * kz1Sq ) )
        !     RV( is + 1 ) = RV( is + 1 ) + DELTA
        !  END IF

        if ( Beam%Type( 4 : 4 ) == 'S' ) then   ! beam displacement & width change (Seongil's version)
           ch = ray2D( is )%c / conjg( HS%cP )
           co = ray2D( is )%t( 1 ) * ray2D( is )%c
           si = ray2D( is )%t( 2 ) * ray2D( is )%c
           ck = omega / ray2D( is )%c

           a   = 2 * HS%rho * ( 1 - ch * ch )
           b   = co * co - ch * ch
           d   = HS%rho * HS%rho * si * si + b
           sb  = sqrt( b )
           cco = co * co
           ssi = si * si

           IF ( si /= 0.0 ) THEN
              delta = a * co / si / ( ck * sb * d )   ! Do we need an abs() on this???
           ELSE
              delta = 0.0
           END IF

           pdelta  = real( delta ) / ( ray2D( is )%c / co)
           ddelta  = -a / ( ck*sb*d ) - a*cco / ssi / (ck*sb*d) + a*cco / (ck*b*sb*d) &
                     -a*co / si / (ck*sb*d*d) * (2* HS%rho * HS%rho *si*co-2*co*si)
           rddelta = -real( ddelta )
           sddelta = rddelta / abs( rddelta )        

           ! next 3 lines have an update by Diana McCammon to allow a sloping bottom
           ! I think the formulas are good, but this won't be reliable because it doesn't have the logic
           ! that tracks crossing into new segments after the ray displacement.
           
           theta_bot = datan( tBdry( 2 ) / tBdry( 1 ))  ! bottom angle
           ray2D( is1 )%x( 1 ) = ray2D( is1 )%x( 1 ) + real( delta ) * dcos( theta_bot )       ! range displacement
           ray2D( is1 )%x( 2 ) = ray2D( is1 )%x( 2 ) + real( delta ) * dsin( theta_bot )       ! depth displacement
           ray2D( is1 )%tau    = ray2D( is1 )%tau + pdelta                                     ! phase change
           ray2D( is1 )%q      = ray2D( is1 )%q + sddelta * rddelta * si * c * ray2D( is )%p   ! beam-width change
        endif

     ENDIF
  CASE DEFAULT
     WRITE( PRTFile, * ) 'HS%BC = ', HS%BC
     CALL ERROUT( 'Reflect2D', 'Unknown boundary condition type' )
  END SELECT

END SUBROUTINE Reflect2D

END PROGRAM BELLHOP
