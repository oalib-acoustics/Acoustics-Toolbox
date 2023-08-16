PROGRAM SCOOTER

  ! Finite-element wavenumber integration program

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

  ! Initial version developed at the Naval Research Laboratory in 1985.

  USE MathConstants
  USE SourceReceiverPositions
  USE sspMod
  USE AttenMod
  USE RWSHDFile
  USE FatalError

  IMPLICIT NONE

  SAVE
  INTEGER, PARAMETER :: ENVFile = 5, PRTFile = 6, GRNFile = 25, IRCFile = 12, MaxMedium = 500
  INTEGER            :: FirstAcoustic, LastAcoustic, Nk, ii, j
  REAL      (KIND=8) :: cMin, cLow, cHigh, omega, omega2, RMax, Atten
  CHARACTER  (LEN=8) :: TopOpt, BotOpt

  ! finite-difference grid
  INTEGER              :: Loc( MaxMedium ), N( MaxMedium ), NG( MaxMedium )
  REAL      (KIND=8)   :: h(   MaxMedium )

  ! storage for finite-difference equations
  REAL    (KIND=8),    ALLOCATABLE :: rho( : )
  COMPLEX (KIND=8),    ALLOCATABLE :: B1( : ), B2( : ), B3( : ), B4( : )

  INTEGER              :: ik, NPoints, NTotal1, IAllocStat, ifreq
  REAL (KIND=8)        :: kMin, kMax, deltak, TStart, TEnd
  REAL (KIND=8), ALLOCATABLE :: k( : ), cVec( : )   ! wavenumbers and their phase speeds
  REAL     (KIND=8)    :: freq, freq0
  COMPLEX, ALLOCATABLE :: Green( :, :, : )     ! G( NSz, NRz, Nk )
  CHARACTER (LEN=80)   :: Title, FileRoot
  CHARACTER (LEN=10)   :: PlotType = 'Green'

  CALL CPU_TIME( Tstart )

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  CALL GetPar

  IF ( SSP%NMedia > 1 ) THEN
     IF ( ANY( SSP%sigma( 2 : SSP%NMedia ) /= 0.0 ) ) CALL ERROUT( 'SCOOTER', 'Rough interfaces not allowed' )
  ENDIF

  ! Set up vector of wavenumber samples based on highest frequency
  kMin = REAL( 2 * pi * freqVec( Nfreq ) / cHigh )
  kMAX = REAL( 2 * pi * freqVec( Nfreq ) / cLow  )
  Nk   = INT( 2000.0 * RMax * ( kMax - kMin ) / REAL( pi ) )

  WRITE( PRTFile, * ) 'Nk = ', Nk
  IF ( ALLOCATED( k     ) ) DEALLOCATE( k, cVec, Pos%Rr )
  ALLOCATE( k( Nk ), cVec( Nk ), Pos%Rr( Nk ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( 'SCOOTER', 'Insufficient memory to allocate k( Nk ) vector' )

  ! Set-up the vector of k-space points
  Deltak = ( kMax - kMin ) / ( Nk - 1 )
  k(    1 : Nk )   = kMin + [ ( ik, ik = 0, Nk - 1 ) ] * Deltak
  cVec( 1 : Nk )   = 2 * pi * freqVec( Nfreq ) / k   ! phase speeds associated with wavenumbers

  IF ( TopOpt( 7 : 7 ) == '0' ) THEN   ! Option to zero out the stabilizing attenuation
     WRITE( PRTFile, * ) 'Option selected to zero out stabilizing attenuation'
  END IF

  IF ( ALLOCATED( Green ) ) DEALLOCATE( Green )
  ALLOCATE ( Green( Pos%NSz, Pos%NRz, Nk ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) THEN
     CALL ERROUT( 'SCOOTER', 'Insufficient memory to allocate Green''s function matrix; reduce Rmax, NSz, or NRz' )
  END IF

  FreqLoop: DO ifreq = 1, Nfreq
     freq   = freqVec( ifreq )
     omega  = 2.0D0 * pi * freq
     omega2 = omega ** 2

     IF ( Nfreq > 1 ) THEN
        WRITE( PRTFile, * ) '__________________________________________________________________________'
        WRITE( PRTFile, * )
        WRITE( PRTFile, "( ' Frequency = ', G11.4, ' Hz' )" ) freq
        FLUSH( PRTFile )
     END IF

     CALL UpdateSSPLoss( freq, freq0 )
     CALL UpdateHSLoss(  freq, freq0 )

     N( 1 : SSP%NMedia ) = INT( freq / freq0 * NG( 1 : SSP%NMedia ) )   ! scale 'global' grid in proportion to frequency

     ! following done as a loop to avoid mixing vector and scalar in the MAX function
     DO iSSP = 1, SSP%NMedia
        N( iSSP ) = MAX( N( iSSP ), 100 )   ! make sure the grid has at least 100 points
     END DO

     h( 1 : SSP%NMedia ) = ( SSP%Depth( 2 : SSP%NMedia + 1 ) - SSP%Depth( 1 : SSP%NMedia ) ) / N( 1 : SSP%NMedia )  ! vector of mesh widths
     NPoints             = SUM( N( 1 : SSP%NMedia ) ) + SSP%NMedia                                  ! number of solution points

     IF ( ALLOCATED( B1 ) ) DEALLOCATE( B1, B2, B3, B4, rho )
     ALLOCATE ( B1( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints ), rho( NPoints ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) CALL ERROUT( 'SCOOTER', 'Insufficient memory to allocate B1, B2, B3, B4 vectors' )

     CALL Initialize( NPoints )    ! Initialize matrices

     ! Set up vector of wavenumber samples, k
     kMin   = omega / cHigh
     kMAX   = omega / cLow
     Deltak = ( kMax - kMin ) / ( Nk - 1 )

     k( 1 : Nk ) = kMin + [ ( ik, ik = 0, Nk - 1 ) ] * Deltak

     Atten  = Deltak   ! stabilizing attenuation
     IF ( TopOpt( 7 : 7 ) == '0' ) Atten = 0.0   ! Option to zero out the stabilizing attenuation

     ! Write header for Green's function file, if this is the first frequency
     IF ( ifreq == 1 ) THEN
        Pos%Rr  = cVec   ! phase speed vector goes where r (its conjugate variable) is normally stored in the file
        Pos%NRr = Nk
        CALL WriteHeader( TRIM( FileRoot ) // '.grn', Title, freq0, Atten, PlotType )
     END IF

     NTotal1 = SUM( N( FirstAcoustic : LastAcoustic ) ) + 1       ! size of matrix for acoustic part
     CALL CalculateKernel( NPoints, k, Green, NTotal1, ifreq )

     CALL CPU_TIME( Tend )
     WRITE( PRTFile, "(' CPU Time: ', G15.5, 's   ---cumulative' )" ) Tend - Tstart
  END DO FreqLoop

  CLOSE( ENVFile )
  CLOSE( PRTFile )
  CLOSE( GRNFile )

CONTAINS

  !**********************************************************************!

  SUBROUTINE GetPar

    ! Read in the ENVFile data

    USE ReadEnvironmentMod
    USE RefCoef

    INTEGER, PARAMETER :: iProf = 1
    REAL               :: zMin, zMAX

    Title = 'SCOOTER- '
    CALL ReadEnvironment( FileRoot, Title, freq0, MaxMedium, TopOpt, NG, BotOpt, cLow, cHigh, RMax, ENVFile, PRTFile )

    IF ( RMax <= 0.0 ) CALL ERROUT( 'GetPar', 'RMax must be positive'  )

    zMin = REAL( SSP%Depth( 1              ) )
    zMAX = REAL( SSP%Depth( SSP%NMedia + 1 ) )
    CALL ReadSzRz(    zMin, zMAX )   ! Read source/receiver depths
    CALL ReadfreqVec( freq0, TopOpt( 6 : 6 ) )

    CLOSE ( ENVFile )

    CALL ReadReflectionCoefficient( FileRoot, BotOpt( 1 : 1 ), TopOpt( 2 : 2 ), PRTFile )   ! Bot, Top refl. coef.

  END SUBROUTINE GetPar
  
  !**********************************************************************!

  SUBROUTINE Initialize( NPoints )

    ! Initializes arrays defining difference equations

    INTEGER, INTENT( IN ) :: NPoints
    LOGICAL               :: ElasticFlag = .FALSE.
    INTEGER               :: Medium, N1
    REAL     (KIND=8)     :: cMinV, two_h, freq
    COMPLEX  (KIND=8)     :: cP( NPoints ), cS( NPoints ), cP2, cS2
    CHARACTER (LEN=8)     :: Task

    cMin          = HUGE( cMin )
    FirstAcoustic = 0
    Loc( 1 )      = 0

    MediumLoop: DO Medium = 1, SSP%NMedia
       IF ( Medium /= 1 ) Loc( Medium ) = Loc( Medium - 1 ) + N( Medium - 1 ) + 1
       N1   = N(   Medium ) + 1
       ii   = Loc( Medium ) + 1
       Task = 'TAB'
       CALL EvaluateSSP( cP( ii ), cS( ii ), rho( ii ), Medium, N1, freq, Task )

       IF ( cS( ii ) == ( 0.0, 0.0 ) ) THEN   ! Acoustic medium
          SSP%Material( Medium )  =  'ACOUSTIC'
          IF ( FirstAcoustic == 0 ) FirstAcoustic = Medium
          LastAcoustic = Medium

          cMinV = MINVAL( DBLE( cP( ii:ii + N( Medium ) ) ) )
          cMin  = MIN( cMin, cMinV )
          B1( ii : ii + N( Medium ) ) = omega2 / cP( ii : ii + N( Medium ) ) ** 2
       ELSE                                   ! Elastic medium
          SSP%Material( Medium ) = 'ELASTIC'
          two_h         = 2.0 * h( Medium )

          DO j = ii, ii + N( Medium )
             cMin = MIN( DBLE( cS( j ) ), cMin )

             cP2 = cP( j ) ** 2
             cS2 = cS( j ) ** 2

             B1(  j ) = two_h / ( rho( j ) * cS2 )
             B2(  j ) = two_h / ( rho( j ) * cP2 )
             B3(  j ) = 4.0 * two_h * rho( j ) * cS2 * ( cP2 - cS2 ) / cP2
             B4(  j ) = two_h * ( cP2 - 2.0 * cS2 ) / cP2
             rho( j ) = two_h * omega2 * rho( j )
          END DO

       ENDIF
    END DO MediumLoop

    ! (cLow, cHigh) = phase speed interval for the wavenumber integration
    ! user specified interval is reduced if it exceeds domain of possible phase speeds

    ! Bottom properties
    IF ( HSBot%BC == 'A' ) THEN
       IF ( HSBot%cS /= ( 0.0, 0.0 ) ) THEN    ! Elastic  bottom half-space
          ElasticFlag = .TRUE.
          cMin   = MIN( cMin, DBLE( HSBot%cS ) )
       ELSE                                    ! Acoustic bottom half-space
          cMin   = MIN( cMin, DBLE( HSBot%cP ) )
       END IF
    END IF

    ! Top properties
    IF ( HSTop%BC == 'A' ) THEN
       IF ( HSTop%cS /= ( 0.0, 0.0 ) ) THEN    ! Elastic  top half-space
          ElasticFlag = .TRUE.
          cMin   = MIN( cMin, DBLE( HSTop%cS ) )
       ELSE                                    ! Acoustic top half-space
          cMin   = MIN( cMin, DBLE( HSTop%cP ) )
       END IF
    END IF
    
    IF ( ElasticFlag ) THEN
       cMin = 0.85 * cMin   ! phase speed a bit less than Scholte/Stoneley wave speed
    ELSE
       cMin = 0.90 * cMin   ! phase speed must be some 10% below slowest acoustic speed
    END IF

    ! The code above calculates phase speed limits, but I'm not using them at the moment
    ! because the k-vector has already been allocated at this stage
    !!!cLow = MAX( cLow, cMin )
    
  END SUBROUTINE Initialize
  
  !**********************************************************************!

  SUBROUTINE BCImpedance( x, BotTop, HS, f, g, iPower )

    !     Compute Boundary Condition Impedance
    !     Same subroutine as in KRAKENC except
    !        PEKRT    is replaced by SQRT
    !        cInside  is related to B1 differently

    USE RefCoef

    COMPLEX  (KIND=8), PARAMETER :: zero = (0.0D0, 0.0D0 )
    COMPLEX  (KIND=8), INTENT( IN  ) :: x           ! wavenumber squared
    CHARACTER (LEN=3), INTENT( IN  ) :: BotTop      ! flag indicating bottom or top boundary
    TYPE( HSInfo ),    INTENT( IN  ) :: HS          ! structure containing halfspace parameters
    INTEGER,           INTENT( OUT ) :: iPower      ! powers of 10 for large dynamic range of the impedance
    COMPLEX  (KIND=8), INTENT( OUT ) :: f, g        ! functions defining the impedance for that value of x
    TYPE( ReflectionCoef ) :: RInt
    INTEGER                :: Ibot, Itop, Medium
    REAL     (KIND=8)      :: rhoInside = 1.
    COMPLEX  (KIND=8)      :: kx, kz, gammaS2, gammaP2, gammaS, gammaP, mu, yV( 5 ), RCmplx, cInside = 1500.
    ! complex  (kind=8)      :: k0, znorm   ! air acoustics tests for Richard Evans

    iPower = 0

    ! Get rho, C just inside the boundary
    SELECT CASE ( BotTop )
    CASE ( 'TOP' )
       Itop      = 1
       rhoInside = rho( Itop )
       cInside   = SQRT( omega2 / B1( Itop ) )
    CASE ( 'BOT' )
       Ibot      = Loc( LastAcoustic ) + N( LastAcoustic ) + 1
       rhoInside = rho( Ibot )
       cInside   = SQRT( omega2 / B1( Ibot ) )
    END SELECT

    ! impedance for different boundary types

    SELECT CASE ( HS%BC )
    CASE ( 'V' )                   ! Vacuum
       f       = 1.0
       g       = -i * SQRT( omega2 / cInside ** 2 - x ) * SSP%sigma( 1 ) ** 2
       yV      = CMPLX( [ f, g, zero, zero, zero ] )
    CASE ( 'R' )                    ! Rigid
       f       = 0.0
       g       = 1.0
       yV      = CMPLX( [ f, g, zero, zero, zero ] )
    CASE ( 'A' )                    ! Acousto-elastic half-space
       IF ( REAL( HS%cS ) > 0.0 ) THEN
          gammaS2 = x - omega2 / HS%cS ** 2
          gammaP2 = x - omega2 / HS%cP ** 2
          gammaS  = SQRT( gammaS2 )
          gammaP  = SQRT( gammaP2 )
          mu      = HS%rho * HS%cS ** 2

          yV( 1 ) = ( gammaS * gammaP - x ) / mu
          yV( 2 ) = ( ( gammaS2 + x ) ** 2 - 4.0 * gammaS * gammaP * x ) * mu
          yV( 3 ) = 2.0 * gammaS * gammaP - gammaS2 - x
          yV( 4 ) = gammaP * ( x - gammaS2 )
          yV( 5 ) = gammaS * ( gammaS2 - x )

          f = omega2 * yV( 4 )
          g = yV( 2 )
       ELSE
          gammaP = SQRT( x - omega2 / HS%cP ** 2 )
          f = 1.0
          g = HS%rho / gammaP
       ENDIF
    CASE ( 'F' )                    ! Tabulated reflection coefficient
       ! Compute the grazing angle Theta
       kx         = SQRT( x )
       kz         = SQRT( omega2 / cInside ** 2 - kx ** 2 )
       RInt%theta = RadDeg * ATAN2( REAL( kz ), REAL( kx ) )
       ! RInt%theta = RadDeg * ATAN( kz / kx )

       ! Evaluate R( TheInt )
       SELECT CASE ( BotTop )
       CASE ( 'TOP' )
          CALL InterpolateReflectionCoefficient( RInt, RTop, NTopPts, PRTFile )
       CASE ( 'BOT' )
          CALL InterpolateReflectionCoefficient( RInt, RBot, NBotPts, PRTFile )
       END SELECT

       ! Convert R( Theta ) to ( f, g ) in Robin BC
       RCmplx = RInt%R * EXP( i * RInt%phi )

       f      = 1.0
       g      = rhoInside * ( 1.0 + RCmplx ) / ( i * kz * ( 1.0 - RCmplx ) )

    CASE ( 'P' )                    ! Precalculated reflection coef
       CALL InterpolateIRC( x, f, g, iPower, xTab, fTab, gTab, ITab, NkTab )
    END SELECT

    IF ( BotTop == 'TOP' ) g = -g    ! A top BC has the sign flipped relative to a bottom BC

    ! Shoot through elastic layers
    SELECT CASE ( BotTop )
    CASE ( 'TOP' )
       IF ( FirstAcoustic > 1 ) THEN
          DO Medium = 1, FirstAcoustic - 1                ! Shooting down from top
             CALL ElasticDn( x, yV, iPower, Medium )
          END DO
          f = omega2 * yV( 4 )
          g = yV( 2 )
       END IF
    CASE ( 'BOT' )
       IF ( LastAcoustic < SSP%NMedia ) THEN
          DO Medium = SSP%NMedia, LastAcoustic + 1, -1    ! Shooting up from bottom
             CALL ElasticUp( x, yV, iPower, Medium )
          END DO
          f = omega2 * yV( 4 )
          g = yV( 2 )
       END IF
    END SELECT

  END SUBROUTINE BCImpedance

  !**********************************************************************!

  SUBROUTINE ElasticUp( x, yV, iPower, Medium )

    ! Propagates through an elastic layer using compound matrix formulation

    INTEGER,          PARAMETER :: iPowerR = 5, iPowerF = -5
    REAL    (KIND=8), PARAMETER :: Roof = 1.0E5, Floor = 1.0E-5
    INTEGER,          INTENT( IN    ) :: Medium
    INTEGER,          INTENT( INOUT ) :: iPower
    COMPLEX (KIND=8), INTENT( IN    ) :: x                ! wavenumber, k^2
    REAL    (KIND=8)                  :: two_h
    COMPLEX (KIND=8)                  :: xV( 5 ), yV( 5 ), zV( 5 )   ! solution of differential equation at 3 successive steps
    COMPLEX (KIND=8)                  :: two_x, xB3, four_h_x

    ! Euler's method for first step
    two_x    = 2.0 * x
    two_h    = 2.0 * h( Medium )
    four_h_x = 4.0 * h( Medium ) * x
    j        = Loc( Medium ) + N( Medium ) + 1
    xB3      = x * B3( j ) - rho( j )

    zV( 1 ) = yV( 1 ) - 0.5 * (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
    zV( 2 ) = yV( 2 ) - 0.5 * ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
    zV( 3 ) = yV( 3 ) - 0.5 * (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
    zV( 4 ) = yV( 4 ) - 0.5 * (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
    zV( 5 ) = yV( 5 ) - 0.5 * (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

    ! Modified midpoint method
    Step: DO ii = N( Medium ), 1, -1
       j = j - 1

       xV = yV
       yV = zV

       xB3 = x * B3( j ) - rho( j )

       zV( 1 ) = xV( 1 ) - (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
       zV( 2 ) = xV( 2 ) - ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
       zV( 3 ) = xV( 3 ) - (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
       zV( 4 ) = xV( 4 ) - (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
       zV( 5 ) = xV( 5 ) - (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

       ! Scale if necessary
       IF ( ii /= 1 ) THEN
          IF      ( ABS( DBLE( zV( 2 ) ) ) < Floor ) THEN
             zV     = Roof * zV
             yV     = Roof * yV
             iPower = iPower - iPowerR
          ELSE IF ( ABS( DBLE( zV( 2 ) ) ) > Roof ) THEN
             zV     = Floor * zV
             yV     = Floor * yV
             iPower = iPower - iPowerF
          END IF

       END IF
    END DO Step

    yV = ( xV + 2.0 * yV + zV ) / 4.0   ! Apply the standard filter at the terminal point

  END SUBROUTINE ElasticUp

  !**********************************************************************!

  SUBROUTINE ElasticDn( x, yV, iPower, Medium )

    ! Propagates through an elastic layer using compound matrix formulation

    INTEGER,          PARAMETER :: iPowerR = 5, iPowerF = -5
    REAL    (KIND=8), PARAMETER :: Roof = 1.0E5, Floor = 1.0E-5
    INTEGER,          INTENT( IN    ) :: Medium
    INTEGER,          INTENT( INOUT ) :: iPower
    COMPLEX (KIND=8), INTENT( IN    ) :: x                ! wavenumber, k^2
    REAL    (KIND=8) :: two_h
    COMPLEX (KIND=8) :: xV( 5 ), yV( 5 ), zV( 5 )   ! solution of differential equation at 3 successive steps
    COMPLEX (KIND=8) :: two_x, xB3, four_h_x

    ! Euler's method for first step

    two_x    = 2.0 * x
    two_h    = 2.0 * h( Medium )
    four_h_x = 4.0 * h( Medium ) * x
    j        = Loc( Medium ) + 1
    xB3      = x * B3( j ) - rho( j )

    zV( 1 ) = yV( 1 ) + 0.5 * (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
    zV( 2 ) = yV( 2 ) + 0.5 * ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
    zV( 3 ) = yV( 3 ) + 0.5 * (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
    zV( 4 ) = yV( 4 ) + 0.5 * (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
    zV( 5 ) = yV( 5 ) + 0.5 * (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

    ! Modified midpoint method
    Step: DO ii = 1, N( Medium )
       j = j + 1

       xV = yV
       yV = zV

       xB3 = x * B3( j ) - rho( j )

       zV( 1 ) = xV( 1 ) + (   B1( j ) * yV( 4 ) - B2( j ) * yV( 5 ) )
       zV( 2 ) = xV( 2 ) + ( -rho( j ) * yV( 4 ) -     xB3 * yV( 5 ) )
       zV( 3 ) = xV( 3 ) + (     two_h * yV( 4 ) + B4( j ) * yV( 5 ) )
       zV( 4 ) = xV( 4 ) + (       xB3 * yV( 1 ) + B2( j ) * yV( 2 ) - two_x * B4( j ) * yV( 3 ) )
       zV( 5 ) = xV( 5 ) + (  rho( j ) * yV( 1 ) - B1( j ) * yV( 2 ) -        four_h_x * yV( 3 ) )

       ! Scale if necessary
       IF ( ii /= N( Medium ) ) THEN
          IF     ( ABS( DBLE( zV( 2 ) ) ) < Floor ) THEN
             zV     = Roof * zV
             yV     = Roof * yV
             iPower = iPower - iPowerR
          ELSE IF ( ABS( DBLE( zV( 2 ) ) ) > Roof ) THEN
             zV     = Floor * zV
             yV     = Floor * yV
             iPower = iPower - iPowerF
          ENDIF
       ENDIF
    END DO Step

    yV = ( xV + 2.0 * yV + zV ) / 4.0   ! Apply the standard filter at the terminal point

  END SUBROUTINE ElasticDn

  !**********************************************************************!

  SUBROUTINE CalculateKernel( NPoints, k, Green, NTotal1, ifreq )

    ! Solve system for a sequence of k-values

    USE calculateweights

    INTEGER,          INTENT( IN  ) :: NPoints, NTotal1, ifreq
    REAL    (KIND=8), INTENT( IN  ) :: k( Nk )
    COMPLEX,          INTENT( OUT ) :: Green( Pos%NSz, Pos%NRz, Nk )
    INTEGER                         :: ik, iS, iR, l, Medium
    REAL                            :: z( NTotal1 ), rhoElement( NPoints )
    REAL    (KIND=8)                :: rhoh, rhoSz( Pos%NSz )
    COMPLEX (KIND=8)                :: BElement, DF( NTotal1 ), EF( NTotal1 ), x

    ! Tabulate z coordinates
    z( 1 ) = REAL( SSP%Depth( FirstAcoustic ) )
    j      = 2

    DO Medium = FirstAcoustic, LastAcoustic
       ! z( j : j + N( Medium ) - 1 ) = REAL( SSP%Depth( Medium ) + [ ( ii * h( Medium ), ii = 1, N( Medium ) ) ] )
       ! gfortran generated erroneous code for the above so we rewrite as:
       z( j : j + N( Medium ) - 1 ) = REAL( [ ( ii * h( Medium ), ii = 1, N( Medium ) ) ] )
       z( j : j + N( Medium ) - 1 ) = REAL( SSP%Depth( Medium ) ) + z( j : j + N( Medium ) - 1 )
       j = j + N( Medium )
    END DO

    ! Compute weights for source/rcvr depth interpolation
    CALL WEIGHT( z, NTotal1, Pos%Sz, Pos%NSz, Pos%ws, Pos%iSz )
    CALL WEIGHT( z, NTotal1, Pos%Rz, Pos%NRz, Pos%wr, Pos%iRz )

    ! Tabulate the density at each source depth, using a piecewise constant approx. to density
    ! A medium with N f.d. points has N+1 SSP, rho samples
    ! and N elements, so you need to offset by 1 for each medium
    ! The density at an interface between two media is set to that of the shallower medium

    Medium = FirstAcoustic
    DO iS = 1, Pos%NSz
       ! search for medium for this source depth
       DO WHILE ( Pos%Sz( iS ) >= SSP%Depth( Medium + 1 ) .AND. Medium < LastAcoustic )
          Medium = Medium + 1
          !write( *, * ) 'Medium ', Medium,  Pos%Sz( is ), SSP%Depth( Medium )
       END DO

       ii          = Loc( FirstAcoustic ) + 1 + Pos%isz( is ) + ( Medium - 1 )
       rhoSz( iS ) = rho( ii )
       !write( *, * ) iS, Pos%Sz( iS ), ii, rhoSz( iS )
    END DO

    ! Assemble matrix
    j       = 1
    l       = Loc( FirstAcoustic ) + 1
    DF( 1 ) = 0.0

    MediumLoop: DO Medium = FirstAcoustic, LastAcoustic
       DO ii = 1, N( Medium )
          rhoElement( l ) = REAL( rho( l ) + rho( l + 1 ) ) / 2.0
          rhoH            = rhoElement( l ) * h( Medium )
          BElement        = h( Medium ) * ( ( B1( l ) + B1( l + 1 ) ) / 2.0 ) / ( 12.0 * rhoElement( l ) )

          DF( j     ) = DF( j ) - 1.0 / rhoH + 5.0 * BElement
          DF( j + 1 ) =         - 1.0 / rhoH + 5.0 * BElement
          EF( j + 1 ) =           1.0 / rhoH +       BElement

          j = j + 1
          l = l + 1
       END DO
       l = l + 1
    END DO MediumLoop

    DO ik = 1, Nk   ! Step through each point in k-space
       x = ( k( ik ) + i * Atten ) ** 2
       CALL Solve( NTotal1, x, Green, ik, DF, EF, rhoElement, rhoSz )  ! Solve for G( k )
    END DO

    ! Write Green's function to file
    SourceDepth: DO iS = 1, Pos%NSz
       RcvrDepth: DO iR = 1, Pos%NRz
          WRITE( GRNFile, REC = 10 + ( ifreq - 1 )  * Pos%NSz * Pos%NRz + ( iS - 1 ) * Pos%NRz + iR ) Green( iS, iR, : )
       END DO RcvrDepth
    END DO SourceDepth

  END SUBROUTINE CalculateKernel

  !**********************************************************************!

  SUBROUTINE Solve( NTotal1, x, Green, ik, DF, EF, rhoElement, rhoSz )

    ! Set up the linear system and solve

    USE factor_Mod
    USE backsub_Mod
    INTEGER,           INTENT( IN  ) :: NTotal1, ik
    REAL,              INTENT( IN  ) :: rhoElement( * )
    REAL     (KIND=8), INTENT( IN  ) :: rhoSz( Pos%NSz )
    COMPLEX  (KIND=8), INTENT( IN  ) :: DF( * ), EF( * ), x
    COMPLEX,           INTENT( OUT ) :: Green( Pos%NSz, Pos%NRz, Nk )
    INTEGER                          :: elt, Medium, is, iPower
    COMPLEX  (KIND=8)                :: d( NTotal1 ), e( NTotal1 ), RV1( NTotal1 ), RV2( NTotal1 ), RV4( NTotal1 )
    COMPLEX  (KIND=8)                :: BElement, xT, f, g

    ! Complete assembly of matrix by adding in x
    j   = 1
    elt = Loc( FirstAcoustic ) + 1
    d( 1 ) = DF( 1 )

    MediumLoop: DO Medium = FirstAcoustic, LastAcoustic
       xT = -h( Medium ) * x / 12.0

       DO ii = 1, N( Medium )
          BElement = xT / rhoElement( elt )

          d( j     ) = d(  j     ) + 5.0 * BElement
          d( j + 1 ) = DF( j + 1 ) + 5.0 * BElement
          e( j + 1 ) = EF( j + 1 ) +       BElement

          j   = j + 1
          elt = elt + 1
       END DO

       elt = elt + 1
    END DO MediumLoop

    ! Corner elt requires top impedance

    CALL BCImpedance( x, 'TOP', HSTop, f, g, iPower )
    IF ( g == 0.0 ) THEN
       d( 1 ) = 1.0D0
       e( 2 ) = 0.0D0
    ELSE
       d( 1 ) = d( 1 ) +  f / g
    ENDIF

    ! Corner elt requires bottom impedance

    CALL BCImpedance( x, 'BOT', HSBot, f, g, iPower )
    IF ( g == 0.0 ) THEN
       d( NTotal1 ) = 1.0D0
       e( NTotal1 ) = 0.0D0
    ELSE
       d( NTotal1 ) =  d( NTotal1 ) - f / g
    ENDIF

    CALL Factor( NTotal1, d, e, RV1, RV2, RV4 )     ! Do LU decomposition

    SourceDepth: DO iS = 1, Pos%NSz

       ! Set up RHS in d (previously used for diagonal)
       d           = 0.0
       ii          = Pos%isz( iS )

       ! we need to get the density at the source
       d( ii     ) = 2.0 * ( 1.0 - Pos%ws( iS ) ) / rhosz( iS )
       d( ii + 1 ) = 2.0 *     Pos%ws( iS )       / rhosz( iS )

       CALL BackSub( NTotal1, RV1, RV2, RV4, d )    ! Solve the system

       Green( iS, :, ik ) = CMPLX( d( Pos%iRz ) + Pos%wr * ( d( Pos%iRz + 1 ) - d( Pos%iRz ) ) )   ! extract the solution at the rcvr depths
    END DO SourceDepth

  END SUBROUTINE Solve
END PROGRAM SCOOTER
