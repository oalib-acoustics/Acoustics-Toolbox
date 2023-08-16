PROGRAM BELLHOP3D

  ! Gaussian beam tracing in three dimensions
  ! Michael B. Porter

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

  ! Original version done at NRL in 1986
  ! Published at SACLANT Undersea Research Center 1989
  !
  ! Converted to Fortran 2003 and greatly enhanced for outside users in 2010
  ! with the support of the Agency for Defense Development, Republic of Korea
  ! Supported also by the U.S. Office of Naval Research
  !
  ! Changes included:
  !    Use of standard BELLHOP input files
  !    Inclusion of multiple sources and receivers
  !    Integrated option for Nx2D or full 3D
  !    Implementation of several standard beam options for both the Nx2D and 3D cases
  !    Reading in:
  !       oceanography from a SSP file
  !       bathymetry or altimetry  from a BTY or ATI file
  !       a reflection coefficient from a BRC or TRC file
  !       a source beam pattern from a SBP file

  ! Loose ends:
  !    Nx2D version does not handle jumps in cx, cy (routine step2d)
  !    cannot specify isingle( 2 ) for alpha and beta
  !    Trilinear interpolation (hex) ignores cxy values.
  !    If the water depth for the lower half space is much larger than that for the SSP
  !    in the water column, it will use a step that is too large and exit the ray box.
 
  !    Cerveny beams (rarely used):
  !       influenceC calls SSP; need to select EvaluateSSP2D or EvaluateSSP3D for that to work in BELLHOP3D
  !       efficiency changes for Cerveny beams as well (and in BELLHOP)
  !       space filling or minimum width options are applied to both alpha and beta--- often only want space filling in azimuth

  !    Influend3DGeoHat writes no eigenray info if number of receiver ranges NR=1
  !    r( 1 ) = 1 m in BELLHOP plus logic for reversing
  !    geogaussian should pre-calculate dtauds and dqds like geohat?
  !    fix automatic deltas selection
  !    detect Sz below bottom with bty file

  !    Should probably use a Structure of Arrays instead of an Array of Structures for speed

  !    Desired additional features:
  !       Terrain following option for receiver depth
  !       Variable bottom type vs.lat/long

  USE bellhopMod
  USE ReadEnvironmentBell
  USE RefCoef
  USE bdry3DMod
  USE BeamPattern
  USE sspMod
  USE influence
  USE Influence3D
  USE FatalError

  IMPLICIT NONE
  CHARACTER ( LEN=80 ) :: FileRoot

  ThreeD = .TRUE.

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )
  ! Open the print file
  OPEN( UNIT = PRTFile, FILE = TRIM( FileRoot ) // '.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )

  ! Read in control data

  CALL ReadEnvironment(    FileRoot, ThreeD )  
  CALL ReadATI3D( FileRoot, Bdry%Top%HS%Opt( 5 : 5 ), Bdry%Top%HS%Depth, PRTFile )    ! AlTImetry
  CALL ReadBTY3D( FileRoot, Bdry%Bot%HS%Opt( 2 : 2 ), Bdry%Bot%HS%Depth, PRTFile )    ! BaThYmetry

  CALL ReadReflectionCoefficient(  FileRoot, Bdry%Bot%HS%Opt( 1 : 1 ), Bdry%Top%HS%Opt( 2 : 2 ), PRTFile )    ! (top and bottom)
  SBPFlag = Beam%RunType( 3 : 3 )
  CALL ReadPAT( FileRoot,                                 PRTFile )    ! Source Beam Pattern
  CALL OpenOutputFiles( FileRoot, ThreeD )

  CALL BellhopCore

  CONTAINS

! **********************************************************************!

SUBROUTINE BellhopCore

  USE angleMod
  USE SourceReceiverPositions
  USE ArrMod
  USE WriteRay

  INTEGER,   PARAMETER :: ArrivalsStorage = 400000000
  INTEGER              :: IBPvec( 1 ), ibp, iBeamWindow2, irz, itheta, isx, isy, isz, iRec, ir
  REAL      ( KIND=8 ) :: Tstart, Tstop
  REAL      ( KIND=8 ) :: Amp0, RadMax, S
  REAL      ( KIND=8 ) :: c0, cimag0, gradc( 3 ), cxx, cyy, czz, cxy, cxz, cyz, rho
  REAL      ( KIND=8 ), ALLOCATABLE :: x_rcvrMat( :, :, : ), t_rcvr( :, : )
  COMPLEX   ( KIND=8 ) :: epsilon( 2 )
  COMPLEX, ALLOCATABLE :: P( :, :, : ), U( :, : )

  CALL CPU_TIME( Tstart )

  omega = 2.0 * pi * freq

  IF ( Beam%deltas == 0.0 ) Beam%deltas = ( Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth ) / 10.0   ! Automatic step size selection

  Angles%alpha  = DegRad * Angles%alpha   ! convert to radians
  Angles%Dalpha = 0.0
  IF ( Angles%Nalpha /= 1 ) &
     Angles%Dalpha = ( Angles%alpha( Angles%Nalpha ) - Angles%alpha( 1 ) ) / ( Angles%Nalpha - 1 )  ! angular spacing between beams

  SELECT CASE ( Beam%RunType( 5 : 5 ) )
  CASE ( 'I' )
     NRz_per_range = 1         ! irregular grid
  CASE ( 'R' )
     NRz_per_range = Pos%NRz   ! rectilinear grid
  END SELECT

  ! for a TL calculation, allocate space for the pressure matrix
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'C', 'S', 'I' )        ! TL calculation
     ALLOCATE ( P( Pos%Ntheta, NRz_per_range, Pos%NRr ), Stat = IAllocStat )
     ALLOCATE ( U( NRz_per_range, Pos%NRr ), Stat = IAllocStat )    ! used for 2D option
     IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'BELLHOP', 'Insufficient memory for TL matrix: reduce Nr * NRz'  )
  CASE ( 'A', 'a', 'R', 'E' )   ! Arrivals calculation
     ALLOCATE ( P( 1, 1, 1 ), Stat = IAllocStat )   ! open a dummy variable
     ALLOCATE ( U( 1, 1 ),    Stat = IAllocStat )   ! open a dummy variable
  END SELECT

  ! for an arrivals run, allocate space for arrivals matrices
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'A', 'a' )
     MaxNArr = MAX( ArrivalsStorage / ( Pos%Ntheta * NRz_per_range * Pos%NRr ), 10 )   ! allow space for at least 10 arrivals
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) '( Maximum # of arrivals = ', MaxNArr, ')'

     ALLOCATE ( Arr3D( Pos%Ntheta, NRz_per_range, Pos%NRr, MaxNArr ), &
               NArr3D( Pos%Ntheta, NRz_per_range, Pos%NRr ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP', &
          'Insufficient memory to allocate arrivals matrix; reduce parameter ArrivalsStorage' )

     ! For a 2D Arrivals run, also need to allocate a structure for that
     IF ( Beam%RunType( 6 : 6 ) == '2' ) THEN
         ALLOCATE ( Arr( NRz_per_range, Pos%NRr, MaxNArr ), &
                   NArr( NRz_per_range, Pos%NRr ), Stat = IAllocStat )
         IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP', &
            'Insufficient memory to allocate arrivals matrix; reduce parameter ArrivalsStorage' )
         NArr( 1 : NRz_per_range, 1 : Pos%NRr ) = 0
     END IF
  CASE DEFAULT
     MaxNArr = 1
     ALLOCATE ( Arr3D( Pos%Ntheta, NRz_per_range, Pos%NRr, 1 ), NArr3D( Pos%Ntheta, NRz_per_range, Pos%NRr ), Stat = IAllocStat )
  END SELECT

  ALLOCATE( x_rcvrMat( 2, Pos%Ntheta, Pos%NRr ), t_rcvr( 2, Pos%Ntheta ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( 'BELLHOP', &
       'Insufficient memory to x_rcvrMat; reduce Nr or Ntheta' )

  ! tangent along receiver bearing line
  t_rcvr( 1, : ) = COS( DegRad * Pos%theta( 1 : Pos%Ntheta ) )
  t_rcvr( 2, : ) = SIN( DegRad * Pos%theta( 1 : Pos%Ntheta ) )

  WRITE( PRTFile, * )

  Source_z: DO isz = 1, Pos%NSz         ! loop over source z-coordinate

     Source_x: DO isx = 1, Pos%Nsx      ! loop over source x-coordinate

        Source_y: DO isy = 1, Pos%Nsy   ! loop over source y-coordinate
           P      = 0.0 ! zero out field matrix
           U      = 0.0
           NArr3D = 0   ! zero out arrival matrix by zeroing out number of arrivals

           ! IF ( r( 1 ) == 0.0 ) r( 1 ) = 1.0
           xs_3D = [ Pos%sx( isx ), Pos%sy( isy ), DBLE( Pos%sz( isz ) ) ]
           WRITE( PRTFile, * )
           WRITE( PRTFile, "( 'xs = ', G11.3, 2X, G11.3, 2X, G11.3 )" ) xs_3D

           ! positions of rcvrs in the x-y plane; this is pre-calculated for InfluenceGeoHatCart
           ! It is not clear that the pre-calculation saves time ...
           DO ir = 1, Pos%NRr
              DO itheta = 1, Pos%Ntheta
                 x_rcvrMat( 1 : 2, itheta, ir ) = xs_3D( 1 : 2 ) + Pos%Rr( ir ) * t_rcvr( :, itheta )  ! x-y coordinate of the receiver
              END DO
           END DO
  
           ! *** Compute 'optimal' beam constant ***

           CALL EvaluateSSP3D( xs_3D, c0, cimag0, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, 'TAB' )
           ray2D( 1 )%c = c0
           ray3D( 1 )%c = c0
           CALL PickEpsilon( Beam%Type( 1 : 2 ), omega, c0, Angles%Dalpha, Angles%Dbeta, Beam%rLoop, Beam%epsMultiplier, epsilon ) ! beam constant

           ! *** Trace successive beams ***

           AzimuthalAngle: DO ibeta = 1, Angles%Nbeta ! this is also the receiver bearing angle for a 2D run 
              SrcAzimAngle = RadDeg * Angles%beta( ibeta )           ! take-off azimuthal   angle in degrees
              ! if ( ibeta /= 134 ) cycle AzimuthalAngle
              IF ( Angles%iSingle_beta == 0 .OR. ibeta == Angles%iSingle_beta ) THEN    ! Single beam run?
              !IF ( ibeta == 2 ) THEN    ! Single beam run?
              !IF ( mod( ibeta+1, 2 ) == 0 ) THEN    ! Single beam run?
                WRITE( PRTFile, FMT = "( 'Tracing azimuthal beam ', I4, F10.2 )" ) ibeta, SrcAzimAngle
                FLUSH( PRTFile )
                 ! WRITE( *,       FMT = "( 'Tracing beam ', I4, F10.2 )" ) ibeta, RadDeg * Angles%beta( ibeta )

                 DeclinationAngle: DO ialpha = 1, Angles%Nalpha
                    SrcDeclAngle = RadDeg *Angles%alpha( ialpha )          ! take-off declination angle in degrees

                    IF ( Angles%iSingle_alpha == 0 .OR. ialpha == Angles%iSingle_alpha ) THEN    ! Single beam run?
                    !IF ( ialpha  == 11 .or. ialpha == 12 ) THEN    ! Single beam run?

                       !WRITE( PRTFile, FMT = "( '   Tracing declination beam ', I4, F10.2 )" ) ialpha, SrcDeclAngle
                       !flush( prtfile )

                       IBPvec = maxloc( SrcBmPat( :, 1 ), mask = SrcBmPat( :, 1 ) < SrcDeclAngle )  ! index of ray angle in beam pattern
                       IBP    = IBPvec( 1 )
                       IBP    = MAX( IBP, 1 )               ! don't go before beginning of table
                       IBP    = MIN( IBP, NSBPPts - 1 )     ! don't go past end of table

                       ! linear interpolation to get amplitude
                       s    = ( SrcDeclAngle  - SrcBmPat( IBP, 1 ) ) / ( SrcBmPat( IBP + 1, 1 ) - SrcBmPat( IBP, 1 ) )
                       Amp0 = ( 1             -                    s ) * SrcBmPat( IBP, 2 ) + s * SrcBmPat( IBP + 1, 2 )

                       ! Lloyd mirror pattern for semi-coherent option
                       IF ( Beam%RunType( 1 : 1 ) == 'S' ) &
                          Amp0 = Amp0 * SQRT( 2.0 ) * ABS( SIN( omega / c0 * xs_3D( 3 ) * SIN( Angles%alpha( ialpha ) ) ) )

                       SELECT CASE ( Beam%RunType( 6 : 6 ) )   ! flag for 2D or 3D calculation
                       CASE ( '2' )   ! Nx2D calculation, neglecting horizontal refraction
                          CALL TraceRay2D(  Angles%alpha( ialpha ), Angles%beta( ibeta ), Amp0 )
                          IF ( Beam%RunType( 1 : 1 ) /= 'R' ) THEN     ! If not a ray trace run, calculate the field
                             SELECT CASE ( Beam%Type( 1 : 1 ) )
                             CASE ( 'R' )
                                iBeamWindow2 = Beam%iBeamWindow ** 2
                                RadMax       = 50 * c0 / freq  ! 50 wavelength max radius
                                CALL InfluenceCervenyRayCen(   U, epsilon( 1 ), Angles%alpha( ialpha ), IBeamWindow2, RadMax )
                             CASE ( 'C' )
                                CALL ERROUT( 'BELLHOP3D', 'Run Type ''C'' not supported at this time' )
                                iBeamWindow2 = Beam%iBeamWindow ** 2
                                RadMax       = 50 * c0 / freq  ! 50 wavelength max radius
                                CALL InfluenceCervenyCart(     U, epsilon( 1 ), Angles%alpha( ialpha ), IBeamWindow2, RadMax )
                             CASE ( 'g' )
                                CALL InfluenceGeoHatRayCen(    U,       Angles%alpha( ialpha ), Angles%Dalpha )
                             CASE ( 'S' )
                                CALL InfluenceSGB(             U,       Angles%alpha( ialpha ), Angles%Dalpha )
                             CASE ( 'B' )
                                CALL InfluenceGeoGaussianCart( U,       Angles%alpha( ialpha ), Angles%Dalpha )
                             CASE ( 'G', '^', ' ' )
                                CALL InfluenceGeoHatCart(      U,       Angles%alpha( ialpha ), Angles%Dalpha )
                             CASE DEFAULT
                                CALL ERROUT( 'BELLHOP3D', 'Invalid Run Type' )
                             END SELECT
                          END IF

                       CASE ( '3' )   ! full 3D calculation
                          CALL TraceRay3D( Angles%alpha( ialpha ), Angles%beta( ibeta ), epsilon, Amp0 )

                          IF ( Beam%RunType( 1 : 1 ) /= 'R' ) THEN     ! If not a ray trace run, calculate the field

                             SELECT CASE ( Beam%RunType( 2 : 2 ) )
                             CASE ( 'C' )   ! Cerveny style beams
                                CALL ERROUT( 'BELLHOP3D', 'Run Type ''C'' not supported at this time' )

                                ! option: assemble f, g, h from p-q
                                !ray3D( 1 : Beam%Nsteps )%f    =    ray3D( 1 : Beam%Nsteps )%p_tilde( 1 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) - &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_hat(   1 )
                                !ray3D( 1 : Beam%Nsteps )%g    = -( ray3D( 1 : Beam%Nsteps )%p_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_hat(   1 ) - &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_hat(   2 ) )
                                !ray3D( 1 : Beam%Nsteps )%h    =    ray3D( 1 : Beam%Nsteps )%p_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) - &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_hat(   2 )
                                !ray3D( 1 : Beam%Nsteps )%DetP =    ray3D( 1 : Beam%Nsteps )%p_tilde( 1 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_hat(   2 ) - &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%p_hat(   1 )
                                !ray3D( 1 : Beam%Nsteps )%DetQ =    ray3D( 1 : Beam%Nsteps )%q_tilde( 1 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_hat(   2 ) - &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_tilde( 2 ) * &
                                !                                   ray3D( 1 : Beam%Nsteps )%q_hat(   1 )

                                !CALL Influence3D_3D( xs_3D, Angles%alpha( ialpha ), iBeamWindow, P )
                             CASE ( 'g' )        ! Geometric-beams with hat-shape
                                CALL Influence3DGeoHatRayCen(       Angles%alpha( ialpha ), Angles%beta( ibeta ), &
                                                              Angles%Dalpha, Angles%Dbeta, P )
                             CASE ( 'G', '^', ' ' )   ! Geometric-beams with hat-shape in Cartesian coordinates
                                CALL Influence3DGeoHatCart(         Angles%alpha( ialpha ), Angles%beta( ibeta ), &
                                                              Angles%Dalpha, Angles%Dbeta, P, x_rcvrMat, t_rcvr )
                             CASE ( 'b' )        ! Geometric-beams with Gaussian-shape
                                CALL Influence3DGeoGaussianRayCen(  Angles%alpha( ialpha ), Angles%beta( ibeta ), &
                                                              Angles%Dalpha, Angles%Dbeta, P )
                             CASE ( 'B' )        ! Geometric-beams with Gaussian-shape in Cartesian coordiantes
                                CALL Influence3DGeoGaussianCart(    Angles%alpha( ialpha ), Angles%beta( ibeta ), &
                                                              Angles%Dalpha, Angles%Dbeta, P, x_rcvrMat, t_rcvr )
                            CASE DEFAULT
                                CALL ERROUT( 'BELLHOP3D', 'Invalid Run Type' )
                             END SELECT
                          END IF
                       END SELECT

                       ! Optionally dump rays to a disk file
                       IF ( Beam%RunType( 1 : 1 ) == 'R' ) THEN
                          CALL WriteRay3D( Angles%alpha( ialpha ), Angles%beta( ibeta ), Beam%Nsteps )
                       ENDIF
                    ENDIF   ! closes iSingle test
                 END DO DeclinationAngle

                 ! for a 2D TL run, scale the pressure and copy the 2D slice into the radial of the 3D field
                 IF ( Beam%RunType( 6 : 6 ) == '2' ) THEN  ! 2D calculation
                    SELECT CASE ( Beam%RunType( 1 : 1 ) )
                    CASE ( 'C', 'S', 'I' )   ! TL calculation
                       CALL ScalePressure( Angles%Dalpha, ray2D( 1 )%c, Pos%Rr, U, NRz_per_range, Pos%NRr, Beam%RunType, freq )
                       P( ibeta, :, : ) = U
                       U = 0   ! clear out the pressure field on the radial in prep for next radial
                    CASE ( 'A' )             ! arrivals calculation, ascii
                       NArr3D( ibeta, :, :    ) = NArr( :, : )
                       Arr3D(  ibeta, :, :, : ) = Arr(  :, :, : )
                       Arr3D(  ibeta, :, :, : )%SrcAzimAngle  = SNGL( SrcAzimAngle )   ! angle
                       Arr3D(  ibeta, :, :, : )%RcvrAzimAngle = SNGL( SrcAzimAngle )   ! angle (rcvr angle is same as source angle)
                       Narr = 0   ! this clears out the 2D arrival structure
                    CASE ( 'a' )             ! arrivals calculation, binary
                       NArr3D( ibeta, :, :    ) = NArr( :, : )
                       Arr3D(  ibeta, :, :, : ) = Arr(  :, :, : )
                       Arr3D(  ibeta, :, :, : )%SrcAzimAngle  = SNGL( SrcAzimAngle )   ! angle
                       Arr3D(  ibeta, :, :, : )%RcvrAzimAngle = SNGL( SrcAzimAngle )   ! angle (rcvr angle is same as source angle)
                       Narr = 0   ! this clears out the 2D arrival structure
                    END SELECT
                 END IF
              ENDIF   ! closes iSingle test
           END DO AzimuthalAngle

           ! *** Scale the complex pressure field ***

           IF ( Beam%RunType( 6 : 6 ) == '3' ) THEN  ! 3D calculation
              SELECT CASE ( Beam%RunType( 1 : 1 ) )
              CASE ( 'C', 'S', 'I' )   ! TL calculation
                 CALL ScalePressure3D( Angles%Dalpha, Angles%Dbeta, ray2D( 1 )%c, epsilon, P, &
                                       Pos%Ntheta, NRz_per_range, Pos%NRr, Beam%RunType, freq )
              END SELECT
           END IF

           ! Write out the field

           SELECT CASE ( Beam%RunType( 1 : 1 ) )
           CASE ( 'C', 'S', 'I' )   ! TL calculation
              DO irz = 1, Pos%NRz
                 RcvrBearing: DO itheta = 1, Pos%Ntheta
                    iRec = 10 + ( isx    - 1 ) * Pos%Nsy * Pos%Ntheta * Pos%NSz * Pos%NRz + &
                                ( isy    - 1 ) *           Pos%Ntheta * Pos%NSz * Pos%NRz + &
                                ( itheta - 1 ) *                        Pos%NSz * Pos%NRz + &
                                ( isz    - 1 ) *                                  Pos%NRz + irz
                    WRITE( SHDFile, REC = IRec ) P( itheta, irz, 1 : Pos%NRr )

                 END DO RcvrBearing
              END DO
              CASE ( 'A' )             ! arrivals calculation, ascii
                 CALL WriteArrivalsASCII3D(  Pos%Rr, Pos%Ntheta, NRz_per_range, Pos%NRr )
              CASE ( 'a' )             ! arrivals calculation, binary
                 CALL WriteArrivalsBinary3D( Pos%Rr, Pos%Ntheta, NRz_per_range, Pos%NRr )
           END SELECT
        END DO Source_y
     END DO Source_x
  END DO Source_z

  ! close all files
  SELECT CASE ( Beam%RunType( 1 : 1 ) )
  CASE ( 'C', 'S', 'I' )      ! TL calculation
     CLOSE( SHDFile )
  CASE ( 'A', 'a' )           ! arrivals calculation
     CLOSE( ARRFile )
  CASE ( 'R' )                ! ray trace
     CLOSE( RAYFile )
  END SELECT

  ! Display run time
  CALL CPU_TIME( Tstop )
  WRITE( PRTFile, "( /, ' CPU Time = ', G15.3, 's' )" ) Tstop - Tstart

END SUBROUTINE BellhopCore

! **********************************************************************!

SUBROUTINE PickEpsilon( BeamType, omega, c, Dalpha, Dbeta, rLoop, EpsMultiplier, epsilon )

  ! Picks the optimum value for epsilon

  REAL      (KIND=8), INTENT( IN  ) :: omega, c             ! angular frequency, sound speed
  REAL      (KIND=8), INTENT( IN  ) :: Dalpha, Dbeta        ! angular spacing for ray fan
  REAL      (KIND=8), INTENT( IN  ) :: epsMultiplier, Rloop ! multiplier, loop range
  COMPLEX   (KIND=8), INTENT( OUT ) :: epsilon( 2 )         ! beam initial conditions
  CHARACTER (LEN= 2), INTENT( IN  ) :: BeamType
  LOGICAL, SAVE      :: INIFlag = .TRUE.
  REAL      (KIND=8) :: HalfWidth(  2 ) = [ 0.0, 0.0 ]
  COMPLEX   (KIND=8) :: epsilonOpt( 2 ) = [ 0.0, 0.0 ]
  CHARACTER (LEN=80) :: TAG

  SELECT CASE ( BeamType( 1 : 1 ) )
  CASE ( 'C', 'R' )   ! Cerveny beams
     TAG    = 'Cerveny style beam'

     SELECT CASE ( BeamType( 2 : 2 ) )
     CASE ( 'F' )
        TAG            = 'Space filling beams'
        HalfWidth( 1 ) = 2.0 / ( ( omega / c ) * Dalpha )
        HalfWidth( 2 ) = 0.0
        IF ( Dbeta /= 0.0 ) HalfWidth( 2 ) = 2.0 / ( ( omega / c ) * Dbeta )
        epsilonOpt     = i * 0.5 * omega * HalfWidth ** 2
     CASE ( 'M' )
        TAG             = 'Minimum width beams'
        HalfWidth(  1 ) = SQRT( 2.0 * c * 1000.0 * rLoop / omega )
        HalfWidth(  2 ) = HalfWidth( 1 )
        epsilonOPT      = i * 0.5 * omega * HalfWidth **2
     CASE ( 'C' )
        TAG    = 'Cerveny style beam'
     END SELECT

  CASE ( 'g' )
     TAG             = 'Geometric beam, hat-shaped, Ray coord.'
     epsilonOPT      = 0.0
  CASE ( 'G', '^' )
     TAG             = 'Geometric beam, hat-shaped, Cart. coord.'
     epsilonOPT      = 0.0
  CASE ( 'b' )
     TAG             = 'Geometric beam, Gaussian-shaped, Ray coord.'
     epsilonOPT      = 0.0
  CASE ( 'B' )
     TAG             = 'Geometric beam, Gaussian-shaped, Cart. coord.'
     epsilonOPT      = 0.0
  CASE ( 'S' )
     TAG        = 'Simple Gaussian beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2
  END SELECT

  epsilon = epsMultiplier * epsilonOPT

  ! On first call write info to prt file
  IF ( INIFlag ) THEN
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) TAG
     WRITE( PRTFile, * ) 'HalfWidth1  = ', HalfWidth( 1 )
     WRITE( PRTFile, * ) 'HalfWidth2  = ', HalfWidth( 2 )
     WRITE( PRTFile, * ) 'epsilonOPT1 = ', epsilonOPT( 1 )
     WRITE( PRTFile, * ) 'epsilonOPT2 = ', epsilonOPT( 2 )
     WRITE( PRTFile, * ) 'EpsMult     = ', EpsMultiplier
     WRITE( PRTFile, * )
     INIFlag = .FALSE.
  END IF

END SUBROUTINE PickEpsilon


!**********************************************************************!

SUBROUTINE TraceRay2D( alpha, beta, Amp0 )

  ! Traces the beam corresponding to a particular take-off angle

  USE ReflectMod

  REAL     (KIND=8), INTENT( IN ) :: alpha, beta, Amp0 ! initial angles, amplitude
  INTEGER           :: is, is1                   ! index for a step along the ray
  REAL     (KIND=8) :: x( 3 )                    ! ray coordinate
  REAL     (KIND=8) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
  REAL     (KIND=8) :: DistBegTop, DistEndTop, DistBegBot, DistEndBot ! Distances from ray beginning, end to top and bottom
  REAL     (KIND=8) :: tradial( 2 ), BotnInt( 3 ), TopnInt( 3 ), s1, s2
  REAL     (KIND=8) :: z_xx, z_xy, z_yy, kappa_xx, kappa_xy, kappa_yy

  ! *** Initial conditions ***

  iSegr = 1   ! this is not really necessary, but ensures the segment search process begins in the same segment for each ray
  iSegz = 1

  iSmallStepCtr = 0
  tradial = [ COS( beta ), SIN( beta ) ]
  ray2D( 1 )%x = [ 0.0D0, xs_3D( 3 ) ]

  CALL EvaluateSSP2D( ray2D( 1 )%x, c, cimag, gradc, crr, crz, czz, rho, xs_3D, tradial, freq )
  ray2D( 1 )%t         = [ COS( alpha ), SIN( alpha ) ] / c
  ray2D( 1 )%p         = [ 1.0, 0.0 ]
  ray2D( 1 )%q         = [ 0.0, 1.0 ]
  ray2D( 1 )%tau       = 0.0
  ray2D( 1 )%Amp       = Amp0
  ray2D( 1 )%Phase     = 0.0
  ray2D( 1 )%NumTopBnc = 0
  ray2D( 1 )%NumBotBnc = 0

  !IsegTopx = 1   ! this is not really necessary, but ensures the segment search process begins in the same segment for each ray
  !IsegTopy = 1
  !IsegBotx = 1
  !IsegBoty = 1
  iSegx = 1   ! this is not really necessary, but ensures the segment search process begins in the same segment for each ray
  iSegy = 1
  iSegz = 1

  CALL GetTopSeg3D( xs_3D )   ! identify the top    segment above the source
  CALL GetBotSeg3D( xs_3D )   ! identify the bottom segment below the source

  ! Trace the beam (note that Reflect alters the step index is)
  is = 0
  CALL Distances3D( xs_3D, Topx, Botx, Topn, Botn, DistBegTop, DistBegBot )

  IF ( DistBegTop <= 0 .OR. DistBegBot <= 0 ) THEN
     Beam%Nsteps = 1
     RETURN       ! source must be within the medium
  END IF

  Stepping: DO istep = 1, MaxN - 1
     is  = is + 1
     is1 = is + 1

     CALL Step2D( ray2D( is ), ray2D( is1 ), tradial )

     ! convert polar coordinate of ray to x-y coordinate
     x( 1 ) = xs_3D( 1 ) + ray2D( is1 )%x( 1 ) * tradial( 1 )
     x( 2 ) = xs_3D( 2 ) + ray2D( is1 )%x( 1 ) * tradial( 2 )
     x( 3 ) = ray2D( is1 )%x( 2 )

     CALL GetTopSeg3D( x )    ! identify the top    segment above the source
     CALL GetBotSeg3D( x )    ! identify the bottom segment below the source

     IF ( IsegTopx == 0 .OR. IsegTopy == 0 .OR. IsegBotx == 0 .OR. IsegBoty == 0 ) THEN ! we escaped the box
        Beam%Nsteps = is
        EXIT Stepping
     END IF

     ! Reflections?
     ! Tests that ray at step is is inside, and ray at step is+1 is outside
     ! to detect only a crossing from inside to outside
     ! DistBeg is the distance at step is,   which is saved
     ! DistEnd is the distance at step is+1, which needs to be calculated
  
     CALL Distances3D( x, Topx, Botx, Topn, Botn, DistEndTop, DistEndBot )

     IF      ( DistBegTop > 0.0d0 .AND. DistEndTop <= 0.0d0 ) THEN  ! test top reflection
        IF ( atiType == 'C' ) THEN

           x = [ xs_3D( 1 ) + ray2D( is + 1 )%x( 1 ) * tradial( 1 ),   &
                 xs_3D( 2 ) + ray2D( is + 1 )%x( 1 ) * tradial( 2 ),   &
                           ray2D( is + 1 )%x( 2 ) ]

           s1     = ( x( 1 ) - Topx( 1 ) ) / Top_deltax   ! proportional distance along segment
           s2     = ( x( 2 ) - Topx( 2 ) ) / Top_deltay   ! proportional distance along segment

           TopnInt = Top( IsegTopx,     IsegTopy     )%Noden * ( 1 - s1 ) * ( 1 - s2 ) +  &
                     Top( IsegTopx + 1, IsegTopy     )%Noden * ( s1     ) * ( 1 - s2 ) +  &
                     Top( IsegTopx + 1, IsegTopy + 1 )%Noden * ( s1     ) * ( s2     ) +  &
                     Top( IsegTopx,     IsegTopy + 1 )%Noden * ( 1 - s1 ) * ( s2     )

           z_xx = Top( IsegTopx, IsegTopy )%z_xx
           z_xy = Top( IsegTopx, IsegTopy )%z_xy
           z_yy = Top( IsegTopx, IsegTopy )%z_yy

           kappa_xx = Top( IsegTopx, IsegTopy )%kappa_xx
           kappa_xy = Top( IsegTopx, IsegTopy )%kappa_xy
           kappa_yy = Top( IsegTopx, IsegTopy )%kappa_yy
        ELSE
           TopnInt = Topn   ! normal is constant in a segment
           z_xx = 0
           z_xy = 0
           z_yy = 0
           kappa_xx = 0
           kappa_xy = 0
           kappa_yy = 0
        END IF

        CALL Reflect2D( is, Bdry%Top%HS, 'TOP', TopnInt, z_xx, z_xy, z_yy, kappa_xx, kappa_xy, kappa_yy, RTop, NTopPTS, tradial)
        ray2D( is + 1 )%NumTopBnc = ray2D( is )%NumTopBnc + 1

        x = [ xs_3D( 1 ) + ray2D( is + 1 )%x( 1 ) * tradial( 1 ),   &
              xs_3D( 2 ) + ray2D( is + 1 )%x( 1 ) * tradial( 2 ),   &
                        ray2D( is + 1 )%x( 2 ) ]

        CALL Distances3D( x, Topx, Botx, Topn, Botn, DistEndTop, DistEndBot )

     ELSE IF ( DistBegBot > 0.0d0 .AND. DistEndBot <= 0.0d0 ) THEN  ! test bottom reflection
        ! write( *, * ) 'Reflecting', x, Botx
        ! write( *, * ) 'Botn', Botn
        ! write( *, * ) 'Distances', DistEndTop, DistEndBot
        IF ( btyType == 'C' ) THEN

           x = [ xs_3D( 1 ) + ray2D( is + 1 )%x( 1 ) * tradial( 1 ),   &
                 xs_3D( 2 ) + ray2D( is + 1 )%x( 1 ) * tradial( 2 ),   &
                           ray2D( is + 1 )%x( 2 ) ]

           s1     = ( x( 1 ) - Botx( 1 ) ) / Bot_deltax   ! proportional distance along segment
           s2     = ( x( 2 ) - Botx( 2 ) ) / Bot_deltay   ! proportional distance along segment

           BotnInt = Bot( IsegBotx,     IsegBoty     )%Noden * ( 1 - s1 ) * ( 1 - s2 ) +  &
                     Bot( IsegBotx + 1, IsegBoty     )%Noden * ( s1     ) * ( 1 - s2 ) +  &
                     Bot( IsegBotx + 1, IsegBoty + 1 )%Noden * ( s1     ) * ( s2     ) +  &
                     Bot( IsegBotx,     IsegBoty + 1 )%Noden * ( 1 - s1 ) * ( s2     )

           z_xx = Bot( IsegBotx, IsegBoty )%z_xx
           z_xy = Bot( IsegBotx, IsegBoty )%z_xy
           z_yy = Bot( IsegBotx, IsegBoty )%z_yy

           kappa_xx = Bot( IsegBotx, IsegBoty )%kappa_xx
           kappa_xy = Bot( IsegBotx, IsegBoty )%kappa_xy
           kappa_yy = Bot( IsegBotx, IsegBoty )%kappa_yy
        ELSE
           BotnInt = Botn   ! normal is constant in a segment
           z_xx = 0
           z_xy = 0
           z_yy = 0
           kappa_xx = 0
           kappa_xy = 0
           kappa_yy = 0
        END IF

        CALL Reflect2D( is, Bdry%Bot%HS, 'BOT', BotnInt, z_xx, z_xy, z_yy, kappa_xx, kappa_xy, kappa_yy, RBot, NBotPTS, tradial)
        ray2D( is + 1 )%NumBotBnc = ray2D( is )%NumBotBnc + 1

        x = [ xs_3D( 1 ) + ray2D( is + 1 )%x( 1 ) * tradial( 1 ),   &
              xs_3D( 2 ) + ray2D( is + 1 )%x( 1 ) * tradial( 2 ),   &
                        ray2D( is + 1 )%x( 2 ) ]

        CALL Distances3D( x, Topx, Botx, Topn, Botn, DistEndTop, DistEndBot )

     END IF

     ! Has the ray left the box, lost its energy, escaped the boundaries, or exceeded storage limit?
     !!!! this should be modified to have a single box
     !!!! no need to test x( 1 ), for instance, against several limits; calculate one limit in advance
     !!!IF ( ABS( x( 1 ) - xs_3D( 1 ) ) > Beam%Box%x .OR. &
     !!!     ABS( x( 2 ) - xs_3D( 2 ) ) > Beam%Box%y .OR. &
     !!!     ABS( x( 3 ) - xs_3D( 3 ) ) > Beam%Box%z .OR. &
     IF( &
          x( 1 ) < MAX( BotGlobalx( 1            ), TopGlobalx( 1            ) ) .OR. &
          x( 2 ) < MAX( BotGlobaly( 1            ), TopGlobaly( 1            ) ) .OR. &
          x( 1 ) > MIN( BotGlobalx( NBTYPts( 1 ) ), TopGlobalx( NATIPts( 1 ) ) ) .OR. &
          x( 2 ) > MIN( BotGlobaly( NBTYPts( 2 ) ), TopGlobaly( NATIPts( 2 ) ) ) .OR. &
          ray2D( is + 1 )%Amp < 0.005 .OR. &
          ! ray2D( is + 1 )%t( 1 ) < 0  .OR. & ! kills off a backward traveling ray
          iSmallStepCtr > 50 ) THEN
        Beam%Nsteps = is + 1
        EXIT Stepping
     ELSE IF ( is >= MaxN - 3 ) THEN
        WRITE( PRTFile, * ) 'Warning in TraceRay2D : Insufficient storage for ray trajectory'
        WRITE( PRTFile, * ) 'Angles are  alpha = ', alpha * RadDeg, '    beta = ', beta * RadDeg
        Beam%Nsteps = is
        EXIT Stepping
     END IF

     DistBegTop = DistEndTop
     DistBegBot = DistEndBot

  END DO Stepping

END SUBROUTINE TraceRay2D

! **********************************************************************!

SUBROUTINE Step2D( ray0, ray2, tradial )

  ! Does a single step along the ray
  ! x denotes the ray coordinate, (r,z)
  ! t denotes the scaled tangent to the ray (previously (rho, zeta))
  ! c * t would be the unit tangent

  USE Step3DMod

  REAL (KIND=8), INTENT( IN ) :: tradial( 2 )   ! coordinate of source and ray bearing angle
  TYPE( ray2DPt )    :: ray0, ray1, ray2
  INTEGER            :: iSegx0, iSegy0, iSegz0
  REAL     (KIND=8 ) :: gradc0( 2 ), gradc1( 2 ), gradc2( 2 ), rho, &
                        c0, cimag0, crr0, crz0, czz0, csq0, cnn0_csq0, &
                        c1, cimag1, crr1, crz1, czz1, csq1, cnn1_csq1, &
                        c2, cimag2, crr2, crz2, czz2, urayt0( 2 ), urayt1( 2 ), &
                        h, halfh, hw0, hw1, ray2n( 2 ), RM, RN, gradcjump( 2 ), cnjump, csjump, w0, w1, &
                        rayx3D( 3 ), rayt3D( 3 ) 

  ! The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
  ! to the Heun (second order Runge-Kutta method).
  ! However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

  ! *** Phase 1 (an Euler step)

  CALL EvaluateSSP2D( ray0%x, c0, cimag0, gradc0, crr0, crz0, czz0, rho, xs_3D, tradial, freq )

  csq0      = c0 * c0
  cnn0_csq0 = crr0 * ray0%t( 2 )**2 - 2.0 * crz0 * ray0%t( 1 ) * ray0%t( 2 ) + czz0 * ray0%t( 1 )**2

  iSegx0 = iSegx     ! make note of current layer
  iSegy0 = iSegy
  iSegz0 = iSegz

  h = Beam%deltas            ! initially set the step h, to the basic one, deltas

  urayt0 = c0 * ray0%t  ! unit tangent
  rayx3D = [ xs_3D( 1 ) + ray0%x( 1 ) * tradial( 1 ), xs_3D( 2 ) + ray0%x( 1 ) * tradial( 2 ), ray0%x( 2 ) ]
  rayt3D = [           urayt0( 1 ) * tradial( 1 ),           urayt0( 1 ) * tradial( 2 ), urayt0( 2 ) ]

  CALL ReduceStep3D( rayx3D, rayt3D, iSegx0, iSegy0, iSegz0, h ) ! reduce h to land on boundary

  halfh = 0.5 * h   ! first step of the modified polygon method is a half step

  ray1%x = ray0%x + halfh * urayt0
  ray1%t = ray0%t - halfh * gradc0 / csq0
  ray1%p = ray0%p - halfh * cnn0_csq0 * ray0%q
  ray1%q = ray0%q + halfh * c0        * ray0%p

  ! *** Phase 2

  CALL EvaluateSSP2D( ray1%x, c1, cimag1, gradc1, crr1, crz1, czz1, rho, xs_3D, tradial, freq )
  csq1      = c1 * c1
  cnn1_csq1 = crr1 * ray1%t( 2 )**2 - 2.0 * crz1 * ray1%t( 1 ) * ray1%t( 2 ) + czz1 * ray1%t( 1 )**2

  ! The Munk test case with a horizontally launched ray caused problems.
  ! The ray vertexes on an interface and can ping-pong around that interface.
  ! Have to be careful in that case about big changes to the stepsize (that invalidate the leap-frog scheme) in phase II.
  ! A modified Heun or Box method could also work.

  urayt1 = c1 * ray1%t   ! unit tangent
  rayt3D = [           urayt1( 1 ) * tradial( 1 ),           urayt1( 1 ) * tradial( 2 ), urayt1( 2 ) ]

  CALL ReduceStep3D( rayx3D, rayt3D, iSegx0, iSegy0, iSegz0, h ) ! reduce h to land on boundary

  ! use blend of f' based on proportion of a full step used.
  w1  = h / ( 2.0d0 * halfh )
  w0  = 1.0d0 - w1
  hw0 = h * w0
  hw1 = h * w1

  ray2%x   = ray0%x   + hw0 * c0 * ray0%t        + hw1 * c1 * ray1%t
  ray2%t   = ray0%t   - hw0 * gradc0 / csq0      - hw1 * gradc1 / csq1
  ray2%p   = ray0%p   - hw0 * cnn0_csq0 * ray0%q - hw1 * cnn1_csq1 * ray1%q
  ray2%q   = ray0%q   + hw0 * c0        * ray0%p + hw1 * c1        * ray1%p
  ray2%tau = ray0%tau + hw0 / CMPLX( c0, cimag0, KIND=8 ) + hw1 / CMPLX( c1, cimag1, KIND=8 )

  ray2%Amp       = ray0%Amp
  ray2%Phase     = ray0%Phase
  ray2%NumTopBnc = ray0%NumTopBnc
  ray2%NumBotBnc = ray0%NumBotBnc

  ! If we crossed an interface, apply jump condition

  CALL EvaluateSSP2D( ray2%x, c2, cimag2, gradc2, crr2, crz2, czz2, rho, xs_3D, tradial, freq )
  ray2%c = c2

  !!! this needs modifying like the full 3D version to handle jumps in the x-y direction
  IF ( iSegz /= iSegz0 ) THEN
     gradcjump =  gradc2 - gradc0  ! this is precise only for c-linear layers
     ray2n     = [ -ray2%t( 2 ), ray2%t( 1 ) ]

     cnjump    = DOT_PRODUCT( gradcjump, ray2n  )
     csjump    = DOT_PRODUCT( gradcjump, ray2%t )

     RM        = ray2%t( 1 ) / ray2%t( 2 )
     RN        = RM * ( 2 * cnjump - RM * csjump ) / c2
     RN        = -RN
     ray2%p    = ray2%p + ray2%q * RN

  END IF

END SUBROUTINE Step2D

!**********************************************************************!

SUBROUTINE TraceRay3D( alpha, beta, epsilon, Amp0 )

  ! Traces the beam corresponding to a particular take off angle

  USE Step3DMod
  USE Reflect3DMod

  REAL     ( KIND=8 ), INTENT( IN ) :: Amp0  ! source coordinate, initial amplitude
  REAL     ( KIND=8 ), INTENT( IN ) :: alpha, beta    ! take-off angles of the ray
  COMPLEX  ( KIND=8 ), INTENT( IN ) :: epsilon( 2 )   ! beam initial conditions
  INTEGER             :: is, is1
  REAL     ( KIND=8 ) :: DistBegTop, DistEndTop, DistBegBot, DistEndBot, &
                         c, cimag, gradc( 3 ), cxx, cyy, czz, cxy, cxz, cyz, rho   ! soundspeed derivatives
  REAL     ( KIND=8 ) :: TopnInt( 3 ), BotnInt( 3 )
  REAL     ( KIND=8 ) :: s1, s2
  REAL     ( KIND=8 ) :: z_xx, z_xy, z_yy, kappa_xx, kappa_xy, kappa_yy

  ! *** Initial conditions ***

  iSegx = 1   ! this is not really necessary, but ensures the segment search process begins in the same segment for each ray
  iSegy = 1
  iSegz = 1

  iSmallStepCtr = 0
  ray3D( 1 )%x    = xs_3D
  CALL EvaluateSSP3D(  ray3D( 1 )%x, c, cimag, gradc, cxx, cyy, czz, cxy, cxz, cyz, rho, freq, 'TAB' )

  ray3D( 1 )%t    = [ COS( alpha ) * COS( beta ) / c, COS( alpha ) * SIN( beta ) / c, SIN( alpha ) / c ]
  !ray3D( 1 )%f    = epsilon( 2 )
  !ray3D( 1 )%g    = epsilon( 1 )
  !ray3D( 1 )%h    = 0.0
  !ray3D( 1 )%DetP = 1.0
  !ray3D( 1 )%DetQ = epsilon( 1 ) * epsilon( 2 )
  ray3D( 1 )%c    = c
  ray3D( 1 )%phi  = 0.0

  ray3D( 1 )%tau       = 0.0
  ray3D( 1 )%Amp       = Amp0
  ray3D( 1 )%Phase     = 0.0
  ray3D( 1 )%NumTopBnc = 0
  ray3D( 1 )%NumBotBnc = 0

  !ray3D( 1 )%p_tilde = [ 1.0, 0.0 ]   use these for complex Cerveny beams
  !ray3D( 1 )%q_tilde = [ epsilon( 1 ), CMPLX( 0.0, KIND=8 ) ]
  !ray3D( 1 )%p_hat   = [ 0.0, 1.0 ]
  !ray3D( 1 )%q_hat   = [ CMPLX( 0.0, KIND=8 ), epsilon( 2 ) ]

  ray3D( 1 )%p_tilde = [ 1.0, 0.0 ]
  ray3D( 1 )%q_tilde = [ 0.0, 0.0 ]
  ray3D( 1 )%p_hat   = [ 0.0, 1.0 ]
  ray3D( 1 )%q_hat   = [ 0.0, 0.0 ]

  ! dummy BotSeg info to force GetBotSeg to search for the active segment on first call
  xTopSeg = [ +big, -big ]
  yTopSeg = [ +big, -big ]
  xBotSeg = [ +big, -big ]
  yBotSeg = [ +big, -big ]
    
  CALL GetTopSeg3D( xs_3D )   ! identify the top    segment above the source
  CALL GetBotSeg3D( xs_3D )   ! identify the bottom segment below the source
  
  ! Trace the beam (note that Reflect alters the step index is)
  is = 0

  CALL Distances3D( ray3D( 1 )%x, Topx,  Botx, Topn, Botn, DistBegTop, DistBegBot )

  IF ( DistBegTop <= 0 .OR. DistBegBot <= 0 ) THEN
     Beam%Nsteps = 1
     RETURN       ! source must be within the medium
  END IF

  Stepping: DO istep = 1, MaxN - 1
     is  = is + 1
     is1 = is + 1

     CALL Step3D( ray3D( is ), ray3D( is1 ) )
     CALL GetTopSeg3D( ray3D( is1 )%x )   ! identify the top    segment above the source
     CALL GetBotSeg3D( ray3D( is1 )%x )   ! identify the bottom segment below the source

     IF ( IsegTopx == 0 .OR. IsegTopy == 0 .OR. IsegBotx == 0 .OR. IsegBoty == 0 ) THEN ! we escaped the box
        Beam%Nsteps = is
        EXIT Stepping
     END IF

     ! Reflections?
     ! Tests that ray at step is is inside, and ray at step is+1 is outside
     ! to detect only a crossing from inside to outside
     ! DistBeg is the distance at step is, which is saved
     ! DistEnd is the distance at step is+1, which needs to be calculated
  
     CALL Distances3D( ray3D( is1 )%x, Topx, Botx, Topn, Botn, DistEndTop, DistEndBot )

     IF      ( DistBegTop > 0.0d0 .AND. DistEndTop <= 0.0d0 ) THEN  ! test top reflection
        IF ( atiType == 'C' ) THEN
           s1 = ( ray3D( is1 )%x( 1 ) - Topx( 1 ) ) / Top_deltax   ! proportional distance along segment
           s2 = ( ray3D( is1 )%x( 2 ) - Topx( 2 ) ) / Top_deltay   ! proportional distance along segment

           TopnInt = Top( IsegTopx,     IsegTopy     )%Noden * ( 1 - s1 ) * ( 1 - s2 ) +  &
                     Top( IsegTopx + 1, IsegTopy     )%Noden * ( s1     ) * ( 1 - s2 ) +  &
                     Top( IsegTopx + 1, IsegTopy + 1 )%Noden * ( s1     ) * ( s2     ) +  &
                     Top( IsegTopx,     IsegTopy + 1 )%Noden * ( 1 - s1 ) * ( s2     )
           z_xx = Top( IsegTopx, IsegTopy )%z_xx
           z_xy = Top( IsegTopx, IsegTopy )%z_xy
           z_yy = Top( IsegTopx, IsegTopy )%z_yy

           kappa_xx = Top( IsegBotx, IsegBoty )%kappa_xx
           kappa_xy = Top( IsegBotx, IsegBoty )%kappa_xy
           kappa_yy = Top( IsegBotx, IsegBoty )%kappa_yy
        ELSE
           TopnInt = Topn   ! normal is constant in a segment
           z_xx = 0
           z_xy = 0
           z_yy = 0

           kappa_xx = 0
           kappa_xy = 0
           kappa_yy = 0
        END IF

        CALL Reflect3D( is, Bdry%Top%HS, 'TOP', TopnInt, z_xx, z_xy, z_yy, kappa_xx, kappa_xy, kappa_yy, RTop, NTopPTS )
        ray3D( is + 1 )%NumTopBnc = ray3D( is )%NumTopBnc + 1
        CALL Distances3D( ray3D( is + 1 )%x, Topx,  Botx, Topn, Botn, DistEndTop, DistEndBot )

     ELSE IF ( DistBegBot > 0.0d0 .AND. DistEndBot <= 0.0d0 ) THEN  ! test bottom reflection

        IF ( btyType == 'C' ) THEN
           s1 = ( ray3D( is1 )%x( 1 ) - Botx( 1 ) ) / Bot_deltax   ! proportional distance along segment
           s2 = ( ray3D( is1 )%x( 2 ) - Botx( 2 ) ) / Bot_deltay   ! proportional distance along segment

           BotnInt = Bot( IsegBotx,     IsegBoty     )%Noden * ( 1 - s1 ) * ( 1 - s2 ) +  &
                     Bot( IsegBotx + 1, IsegBoty     )%Noden * ( s1     ) * ( 1 - s2 ) +  &
                     Bot( IsegBotx + 1, IsegBoty + 1 )%Noden * ( s1     ) * ( s2     ) +  &
                     Bot( IsegBotx,     IsegBoty + 1 )%Noden * ( 1 - s1 ) * ( s2     )
           z_xx = Bot( IsegBotx, IsegBoty )%z_xx
           z_xy = Bot( IsegBotx, IsegBoty )%z_xy
           z_yy = Bot( IsegBotx, IsegBoty )%z_yy

           kappa_xx = Bot( IsegBotx, IsegBoty )%kappa_xx
           kappa_xy = Bot( IsegBotx, IsegBoty )%kappa_xy
           kappa_yy = Bot( IsegBotx, IsegBoty )%kappa_yy
        ELSE
           BotnInt = Botn   ! normal is constant in a segment
           z_xx = 0
           z_xy = 0
           z_yy = 0

           kappa_xx = 0
           kappa_xy = 0
           kappa_yy = 0
        END IF

        CALL Reflect3D( is, Bdry%Bot%HS, 'BOT', BotnInt, z_xx, z_xy, z_yy, kappa_xx, kappa_xy, kappa_yy, RBot, NBotPTS )
        ray3D( is + 1 )%NumBotBnc = ray3D( is )%NumBotBnc + 1
        CALL Distances3D( ray3D( is + 1 )%x, Topx, Botx, Topn, Botn, DistEndTop, DistEndBot )

     END IF

     ! Has the ray exited the beam box, lost its energy, escaped the BTY, ATI boundaries, or exceeded storage limit?
     IF ( &
          ABS( ray3D( is + 1 )%x( 1 ) - xs_3D( 1 ) ) > Beam%Box%x .OR. &
          ABS( ray3D( is + 1 )%x( 2 ) - xs_3D( 2 ) ) > Beam%Box%y .OR. &
          ABS( ray3D( is + 1 )%x( 3 )              ) > Beam%Box%z .OR. & ! box is centered at z=0
          ray3D( is + 1 )%x( 1 ) < MAX( BotGlobalx( 1            ), TopGlobalx( 1            ) ) .OR. &
          ray3D( is + 1 )%x( 2 ) < MAX( BotGlobaly( 1            ), TopGlobaly( 1            ) ) .OR. &
          ray3D( is + 1 )%x( 1 ) > MIN( BotGlobalx( NBTYPts( 1 ) ), TopGlobalx( NATIPts( 1 ) ) ) .OR. &
          ray3D( is + 1 )%x( 2 ) > MIN( BotGlobaly( NBTYPts( 2 ) ), TopGlobaly( NATIPts( 2 ) ) ) .OR. &
          ray3D( is + 1 )%Amp < 0.005 .OR. &
          iSmallStepCtr > 50 ) THEN
        Beam%Nsteps = is + 1
        EXIT Stepping
     ELSE IF ( is >= MaxN - 3 ) THEN
        Beam%Nsteps = is
        WRITE( PRTFile, * ) 'Warning in TraceRay3D : Insufficient storage for ray trajectory'
        WRITE( PRTFile, * ) 'Angles are  alpha = ', alpha * RadDeg, '    beta = ', beta * RadDeg
        EXIT Stepping
     END IF

     DistBegTop = DistEndTop
     DistBegBot = DistEndBot

  END DO Stepping
END SUBROUTINE TraceRay3D

! **********************************************************************!

SUBROUTINE Distances3D( rayx, Topx, Botx, Topn, Botn, DistTop, DistBot )

  ! Computes distances from ray to boundaries
  ! Formula differs from JKPS because code uses outward pointing normals
  ! Note that Topx, Botx, just need to be any node on the diagonal that divides each square into triangles
  ! In bdry3DMod, the node is selected as the one at the lowest x, y, index and that defines the triangles

  REAL (KIND=8), INTENT( IN  ) :: rayx( 3 )             ! ray coordinate
  REAL (KIND=8), INTENT( IN  ) :: Topx( 3 ), Botx( 3 )  ! top, bottom boundary coordinate for node
  REAL (KIND=8), INTENT( IN  ) :: Topn( 3 ), Botn( 3 )  ! top, bottom boundary normal
  REAL (KIND=8), INTENT( OUT ) :: DistTop, DistBot      ! distance from the ray to top, bottom boundaries 
  REAL (KIND=8)                :: dTop( 3 ), dBot( 3 )

  dTop    = rayx - Topx  ! vector pointing from top    to ray
  dBot    = rayx - Botx  ! vector pointing from bottom to ray
  DistTop = -DOT_PRODUCT( Topn, dTop )
  DistBot = -DOT_PRODUCT( Botn, dBot )

!!$  write( *, * )
!!$  write( *, * ) 'Distances3D', DistBot
!!$  write( *, * ) 'rayx', rayx
!!$  write( *, * ) 'Botx', Botx
!!$  write( *, * ) 'dBot', dBot
!!$  write( *, * ) 'Normal', Botn
!!$  write( *, * )
END SUBROUTINE Distances3D

END PROGRAM BELLHOP3D
