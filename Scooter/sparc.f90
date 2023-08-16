PROGRAM SPARC

  ! Finite-element, time-domain, wavenumber-integration program

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

  ! Initial version developed at the SACLANT Undersea Research Center in 1987.

  USE MathConstants
  USE SourceReceiverPositions
  USE sspMod
  USE RWSHDFile
  USE FatalError

  IMPLICIT NONE

  SAVE
  INTEGER,   PARAMETER :: ENVFile = 5, PRTFile = 6, GRNFile = 25, RTSFile = 35, MaxN = 17000, MaxMedium = 500, MaxIt = 1000000
  INTEGER              :: Loc( MaxMedium), N( MaxMedium ), Nk, NTot1, NTout, ii, ir
  REAL                 :: C2R( MaxN ), C2I( MaxN ), Z( MaxN ), &
       H( MaxMedium ), CMin, CLow, CHigh, omega2, DeltaK, &
       Deltat, CrossT, CMax, TStart, V, TMult, alpha, beta, FMin, FMax
  REAL        (KIND=8) :: rho( MaxN ), freq
  CHARACTER (LEN=8 )   :: TopOpt, BotOpt
  CHARACTER (LEN=4 )   :: Pulse
  CHARACTER (LEN=80)   :: Title
  REAL,    ALLOCATABLE :: Tout( : ), RTSrz( :, : ), RTSrr( :, : )
  REAL (KIND=8), ALLOCATABLE :: k( : )
  COMPLEX, ALLOCATABLE :: Green( :, :, : )
  REAL                 :: Tend
  CHARACTER (LEN=80)   :: FileRoot

  CALL CPU_TIME( Tstart )
  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  CALL GetPar( FileRoot )

  h( 1 : SSP%NMedia ) = REAL( SSP%Depth( 2 : SSP%NMedia + 1 ) - SSP%Depth( 1 : SSP%NMedia ) ) / REAL( N( 1 : SSP%NMedia ) )

  CALL Initialize
  CALL WriteHeaderSparc( FileRoot )
  CALL KERNEL

  CALL CPU_TIME( Tend )
  WRITE( PRTFile, "( ' CPU TIME: ', G15.5, 's' )" ) Tend - Tstart

  CLOSE( ENVFile )
  CLOSE( PRTFile )
  STOP

CONTAINS

  !**********************************************************************

  SUBROUTINE GetPar( FileRoot )

    ! Read in the ENVFile data

    USE ReadEnvironmentMod
    USE SubTabulate
    USE AttenMod

    INTEGER            :: IAllocStat, ik
    REAL               :: zMin, zMax
    REAL      (KIND=8) :: kMin, kMax, RMax, cLowD, cHighD
    CHARACTER (LEN=80) :: FileRoot

    Title = 'SPARC-   '
    CALL ReadEnvironment( FileRoot, Title, freq, MaxMedium, TopOpt, N, BotOpt, cLowD, cHighD, RMax, ENVFile, PRTFile )
    IF ( HSBot%BC == 'F' .OR. HSTop%BC == 'P' ) &
       CALL ERROUT( 'SPARC', 'The option to read a file for the reflection loss is not implemented in SPARC' )

    cLow  = REAL( cLowD  )   ! converting double to single precision
    cHigh = REAL( cHighD )

    zMin = SNGL( SSP%Depth( 1              ) )
    zMax = SNGL( SSP%Depth( SSP%NMedia + 1 ) )
    CALL ReadSzRz( zMin, zMax ) !  Read source/receiver depths, Sz, Rz

    omega2 = SNGL( ( 2.0 * pi * freq ) ** 2 )
    CALL UpdateSSPLoss( freq, freq )
    ! update loss in the halfspaces based on the frequency
    ! depth of 1e20 is a large number to ensure bio loss is not included

    IF ( ( HSTop%BC /= 'V' .AND. HSTop%BC /= 'R' ) .OR. &
         ( HSBot%BC /= 'V' .AND. HSBot%BC /= 'R' ) ) THEN
           CALL ERROUT( 'GETPAR', 'SPARC only allows Vacuum or Rigid boundary conditions'  )
    END IF

    READ(  ENVFile, * ) Pulse
    READ(  ENVFile, * ) fMin, fMax   ! Upper and lower frequency limits
    WRITE( PRTFile, * ) 'fMin, fMax = ', fMin, fMax

    ! integration parameters
    kMin = 2.0 * pi * fMin / cHigh
    kMax = 2.0 * pi * fMax / cLow

    kMin = MAX( 1D-20, kMin )   ! avoid a zero that would produce a divide check when the phase speed is written

    Nk = INT( 1000.0 * RMax * ( kMax - kMin ) / ( 2.0 * pi ) )

    WRITE( PRTFile, * ) 'Nk = ', Nk

    ALLOCATE( k( Nk ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( 'GETPAR', 'Too many pts in k-space'  )

    deltak      = SNGL( ( kMax - kMin ) / ( Nk - 1 ) )
    k( 1 : Nk ) = REAL( kMin ) + [ ( ik, ik = 0, Nk - 1 ) ] * deltak

    SELECT CASE ( Pulse( 1 : 1 ) )
    CASE ( 'P' )
       WRITE( PRTFile, * ) 'Pseudo-gaussian pulse'
    CASE ( 'R' )
       WRITE( PRTFile, * ) 'Ricker wavelet'
    CASE ( 'A' )
       WRITE( PRTFile, * ) 'Approximate Ricker wavelet'
    CASE ( 'S' )
       WRITE( PRTFile, * ) 'Single sine source'
    CASE ( 'H' )
       WRITE( PRTFile, * ) 'Hanning weighted four sine pulse'
    CASE ( 'N' )
       WRITE( PRTFile, * ) 'N-wave pulse'
    CASE ( 'M' )
       WRITE( PRTFile, * ) 'Miracle-wave pulse'
    CASE ( 'G' )
       WRITE( PRTFile, * ) 'Gaussian pulse'
    CASE ( 'F' )
       WRITE( PRTFile, * ) 'Source time series from File'
    CASE ( 'B' )
       WRITE( PRTFile, * ) 'Source time series reversed from file'
    CASE DEFAULT
       CALL ERROUT( 'GETPAR', 'Unknown source type'  )
    END SELECT

    CALL ReadRcvrRanges ! Read receiver ranges, Rr

    IF ( Pos%Rr( Pos%NRr ) > 1000.0 * Rmax ) THEN
       WRITE( PRTFile, * ) ' Pos%Rr( Pos%NRr ) = ', Pos%Rr( Pos%NRr ) / 1000.0, '   RMax = ', RMax
       CALL ERROUT( 'GETPAR', 'Largest receiver range exceeds RMax'  )
    END IF

    ! Read in the output times
    CALL ReadVector( Ntout, tout, 'Output times', 's' )

    ! Integration parameters
    ! alpha = 0.0  lumped     mass matrix,
    !         1.0  consistent mass matrix
    ! beta  = 0.0  standard explicit,
    !         0.25 standard implicit

    WRITE( PRTFile, * )
    READ(  ENVFile, * ) tStart, tMult, alpha, beta, V
    WRITE( PRTFile, * ) 'tStart = ', tStart
    WRITE( PRTFile, * ) 'tMult  = ', tMult
    WRITE( PRTFile, * ) 'alpha  = ', alpha
    WRITE( PRTFile, * ) 'beta   = ', beta
    WRITE( PRTFile, * ) 'V      = ', V

    CLOSE ( ENVFile )

    IF ( ANY( SSP%sigma( 1 : SSP%NMedia ) /= 0.0 ) ) CALL ERROUT( 'GetPar', 'Rough interfaces not allowed'  )

  END SUBROUTINE GetPar

  !**********************************************************************

  SUBROUTINE Initialize

    ! Initializes arrays defining difference equations

    USE calculateweights

    INTEGER           :: J, medium, N1
    REAL              :: cpR, crosst1
    COMPLEX           :: cpT
    REAL     (KIND=8) :: freq
    COMPLEX  (KIND=8) :: cp( MaxN ), cs( MaxN )
    CHARACTER (LEN=8) :: Task

    cMin     =  1.0E6
    cMax     = -1.0E6
    crosst   =  1.0E6
    Loc( 1 ) = 0
    J        = 0

    MediumLoop: DO medium = 1, SSP%NMedia
       IF ( medium /= 1 ) Loc( medium ) = Loc( medium - 1 ) + N ( medium - 1 ) + 1
       N1 = N( medium ) + 1
       IF ( Loc( medium ) + N1 > MaxN ) THEN
          WRITE( PRTFile, * ) 'FATAL ERROR: Insufficient storage for mesh'
          STOP                'FATAL ERROR: Insufficient storage for mesh'
       ENDIF

       Task = 'TAB'
       CALL EvaluateSSP( cp, cs, rho( Loc( medium ) + 1 ),  medium, N1, freq, Task )

       DO ii = 1, N1
          cpR     = REAL( cp( ii ), 4 )
          cMin    = MIN( cpR, cMin )
          cMax    = MAX( cpR, cMax )
          crosst1 = REAL( h( medium ) / cp( ii ), 4 )
          crosst  = AMIN1( crosst1, crosst )

          J        = J + 1
          cpT      = REAL(  cp( ii ), 4 ) ! convert to single precision
          c2R( J ) = REAL(  cpT**2 )  ! CHECK THIS
          c2I( J ) = AIMAG( cpT**2 ) / SQRT( omega2 ) / c2R( J )
       END DO

    END DO MediumLoop

    ! Tabulate z-coordinates
    z( 1 ) = SNGL( SSP%Depth( 1 ) )
    J      = 2
    NTot1  = SUM( N( 1:SSP%NMedia ) ) + 1

    DO medium = 1, SSP%NMedia
       z( J : J + N( medium ) - 1 ) = SNGL( SSP%Depth( medium ) ) + [ ( ii * h( medium ), ii = 1, N( medium ) ) ]
       J = J + N( medium )
    END DO

    CALL WEIGHT( z, NTot1, Pos%sz, Pos%NSz, Pos%ws, Pos%isz ) ! Compute weights for source depth interpolation
    CALL WEIGHT( z, NTot1, Pos%rz, Pos%NRz, Pos%wr, Pos%irz ) ! Compute weights for rcvr   depth interpolation

  END SUBROUTINE Initialize

  !**********************************************************************

  SUBROUTINE KERNEL

    ! Solve system for a sequence of k-values

    INTEGER       :: ik, Itout, iG
    REAL          :: x
    REAL (KIND=8) :: Scale

    ! Allocate and clear matrices for storing time-series
    SELECT CASE ( TopOpt( 5 : 5 ) )
    CASE ( 'D')
       ALLOCATE( RTSrz( Pos%NRz, Ntout ) )
       RTSrz = 0.0
    CASE ( 'R')
       ALLOCATE( RTSrr( Pos%NRr,  Ntout ) )
       RTSRR = 0.0
    CASE DEFAULT
       ALLOCATE( Green( Ntout, Pos%NRz, Nk ) )
    END SELECT

    deltat = REAL( tMult / SQRT( 1.0 / crosst**2 + ( 0.5 * cMax * k( Nk ) ) **2 ) ) ! Courant condition to set time step
    WRITE( PRTFile, * ) 'Time step = ', deltat
    WRITE( PRTFile, * ) 'Estimated fl. pt. ops (billions) = ', ( tout( Ntout ) / deltat ) * Nk * NTot1 / 25000000

    Wavenumber: DO ik = 1, Nk
       x      = REAL( k( ik ) ** 2 )
       deltat = REAL( tMult / SQRT( 1.0 / crosst**2 + ( 0.5 * cMax * k( ik ) ) **2 ) )
       ! IF ( 10*(ik/10) == ik ) WRITE( 6, * ) 'ik, Nk', ik, Nk
       IF ( MOD( ik - 1, max( Nk / 50, 1 ) ) == 0 ) THEN
          WRITE( PRTFile, FMT = "( 'Wavenumber ', I7, ' of ', I7 )" ) ik, Nk
          FLUSH( PRTFile )
       END IF
       CALL MARCH( x, ik )  ! March that component for all time
    END DO Wavenumber

    ! write out the field
    SELECT CASE ( TopOpt( 5 : 5 ) )
    CASE ( 'S' )   ! snapshot
       ! time is treated like a source depth in terms of the file format
       DO Itout = 1, Ntout
          RcvrDepth: DO ir = 1, Pos%NRz
             iG = ( Itout - 1 ) * Pos%NRz + ir   ! Index of source/rcvr combo
             WRITE( GRNFile, REC = 10 + iG ) Green( Itout, ir, 1 : Nk )
          END DO RcvrDepth
       END DO
       CLOSE( GRNFile )
    CASE ( 'D' )   ! RTS (vertical array)
       Scale = 1.0 / SQRT( pi * Pos%Rr( 1 ) )
       DO Itout = 1, Ntout
          WRITE( RTSFile, '( 12G15.6 )' ) tout( Itout ), Scale * RTSrz( 1 : Pos%NRz, Itout )
       END DO
       CLOSE( RTSFile )
    CASE ( 'R' )   ! RTS (horizontal array)
       DO Itout = 1, Ntout
          WRITE( RTSFile, '( 12G15.6 )' ) tout( Itout ), RTSrr( 1 : Pos%NRr, Itout )
       END DO
       CLOSE( RTSFile )
    END SELECT

  END SUBROUTINE KERNEL

  !**********************************************************************

  SUBROUTINE WriteHeaderSparc( FileRoot )

    ! Writes headers for disk files

    CHARACTER (LEN=80), INTENT( IN ) :: FileRoot
    REAL      (KIND=8), PARAMETER    :: Atten = 0   ! no stabilizing attenuation in a Sparc run
    CHARACTER (LEN=10)               :: PlotType = 'Green'

    SELECT CASE ( TopOpt( 5 : 5 ) )
    CASE ( 'S' )          ! snapshot
       Nfreq   = Ntout
       ALLOCATE( FreqVec( Nfreq ) )
       freqVec = tout     ! time vector goes in place of frequency vector

       DEALLOCATE( Pos%Rr )   ! This vector is not used in a snapshot run
       ALLOCATE( Pos%Rr( Nk ) )
       Pos%Rr  = SQRT( omega2 ) / k   ! phase speed vector goes where r (its conjugate variable) is normally stored in the file
       IF ( k( 1 ) == 0.0 ) Pos%Rr( 1 ) = 1e20   ! give infinite phase speed where wavenumber vanishes
       Pos%NRr = Nk

       CALL WriteHeader( TRIM( FileRoot ) // '.grn', Title, freq, Atten, PlotType )
    CASE ( 'R' )          ! Horizontal array
       OPEN ( FILE = TRIM( FileRoot ) // '.rts', UNIT = RTSFile, STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
       WRITE( RTSFile, * ) '''' // Title( 1 : 75 ) // ''''
       WRITE( RTSFile, * ) Pos%NRr, Pos%Rr( 1 : Pos%NRr )
    CASE ( 'D' )          ! Vertical array
       OPEN ( FILE = TRIM( FileRoot ) // '.rts', UNIT = RTSFile, STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
       WRITE( RTSFile, * ) '''' // Title( 1 : 75 ) // ''''
       WRITE( RTSFile, * ) Pos%NRz, Pos%rz( 1 : Pos%NRz )
    END SELECT

  END SUBROUTINE WriteHeaderSparc

  !**********************************************************************

  SUBROUTINE MARCH( x, ik )

    ! March a single spectral component forward in time

    USE factor_Mod

    INTEGER, INTENT( IN ) :: ik    ! index of wavenumber
    REAL,    INTENT( IN ) :: x     ! wavenumber squared
    LOGICAL :: IniFlag
    INTEGER :: Itime, medium, Node, Itout, L
    REAL    :: c2RElt, c2IElt, fLoCut, fHiCut, hElt, rhoElt, rkT, time
    COMPLEX :: AD2( NTot1 ), AE2( NTot1 ), AD1( NTot1 ), AE1( NTot1 ), AD0( NTot1 ), AE0( NTot1 ), &
               RV1( NTot1 ), RV2( NTot1 ), RV4( NTot1 ), U0(  NTot1 ), U1(  NTot1 ), U2(  NTot1 )

    ! Assemble mass and stiffness matrices
    rkT  = SQRT( x )
    Node = 1
    L    = 1

    AD2( 1 ) = 0.0
    AD1( 1 ) = 0.0
    AD0( 1 ) = 0.0

    MediumLoop: DO medium = 1, SSP%NMedia
       hElt = h( medium )

       DO ii = 1, N(medium)
          rhoElt = SNGL( ( rho( L ) + rho( L + 1 ) ) / 2.0 )
          c2RElt =       ( c2R( L ) + c2R( L + 1 ) ) / 2.0
          c2IElt =       ( c2I( L ) + c2I( L + 1 ) ) / 2.0

          CALL CONTRIB( Node, x, rkT, V, rhoElt, c2RElt, c2IElt, hElt, deltat, AD2, AE2, AD1, AE1, AD0, AE0, alpha, beta )

          Node = Node + 1
          L    = L    + 1
       END DO

       L = L + 1
    END DO MediumLoop

    CALL Factor( NTot1, AD2, AE2, RV1, RV2, RV4 )  ! * Factor A2 *

    ! Initialize pressure vectors
    U0 = 0.0
    U1 = 0.0
    U2 = 0.0

    ! Initializate parameters for bandpass filtering of the source time series
    IF ( Pulse( 4 : 4 ) == 'L' .OR. Pulse( 4 : 4 ) == 'B' ) THEN
       fLoCut = rkT * cLow / ( 2.0 * REAL( pi ) )
    ELSE
       fLoCut = 0.0
    ENDIF

    IF ( Pulse( 4 : 4 ) == 'H' .OR. Pulse( 4 : 4 ) == 'B' ) THEN
       fHiCut = rkT * cHigh / ( 2.0 * REAL( pi ) )
    ELSE
       fHiCut = 10.0 * fMax
    ENDIF

    time   = 0.0
    IniFlag = .TRUE.   ! need to tell source routine to refilter for each new value of k

    ! Begin forward march
    Itout = 1
    DO Itime = 1, MaxIT
       time = tStart + ( Itime - 1 ) * deltat
       ! Take a step forward in time
       ! AE2 does not need to be passed: if AE2 exists, then STEP uses an implicit solver using RV1, 2, 4
       CALL STEP( AD0, AE0, AD1, AE1, AD2, RV1, RV2, RV4, U0, U1, U2, time, rkT, fLoCut, fHiCut, IniFlag )
       CALL EXTRACT( U0, U1, ik, time, Itout, rkT ) ! Extract soln for desired receivers
       IF ( Itout > Ntout ) RETURN
    END DO

  END SUBROUTINE MARCH

  !**********************************************************************

  SUBROUTINE CONTRIB( Node, x, rkT, V, rhoElt, c2RElt, c2IElt, hElt, deltat, AD2, AE2, AD1, AE1, AD0, AE0, alpha, beta )

    ! Computes the contribution for the given element

    COMPLEX, PARAMETER     :: i = ( 0.0, 1.0 )
    INTEGER, INTENT( IN  ) :: Node
    REAL,    INTENT( IN  ) :: x, rkT, V, rhoElt, c2Relt, c2IElt, hElt, deltat, alpha, beta
    COMPLEX, INTENT( OUT ) :: AD2( * ), AE2( * ), AD1( * ), AE1( * ), AD0( * ), AE0( * )
    REAL                   :: deltat2, rhoh
    COMPLEX                :: Ke( 2, 2 ), Ce( 2, 2 ), Me( 2, 2 ), y, z

    !  alpha controls lumping and beta controls implicitness
    deltat2 = deltat ** 2
    rhoh    = rhoElt * hElt

    ! * Elemental stiffness matrix ( +k2 term ) *
    y = 1.0 + i * rkT * V * c2IElt
    z = x * ( y - V * V / c2RElt )

    Ke( 1, 1 ) =  y / rhoh + ( 3.0 - alpha ) * hElt * z / rhoElt / 6.0
    Ke( 1, 2 ) = -y / rhoh +         alpha   * hElt * z / rhoElt / 6.0
    Ke( 2, 2 ) =  y / rhoh + ( 3.0 - alpha ) * hElt * z / rhoElt / 6.0

    ! * Elemental damping matrix *
    z = 2.0 * i * rkT * V / c2RElt
    Ce( 1, 1 ) = c2IElt * ( 1.0/rhoh + ( 3.0 - alpha ) * hElt * x / rhoElt / 6.0) &
         &                           + ( 3.0 - alpha ) * hElt * z / rhoElt / 6.0
    Ce( 1, 2 ) = c2IElt * (-1.0/rhoh +         alpha   * hElt * x / rhoElt / 6.0) &
         &                           +         alpha   * hElt * z / rhoElt / 6.0
    Ce( 2, 2 ) = c2IElt * ( 1.0/rhoh + ( 3.0 - alpha ) * hElt * x / rhoElt / 6.0) &
         &                           + ( 3.0 - alpha ) * hElt * z / rhoElt / 6.0

    ! * Elemental mass matrix *
    Me( 1, 1 ) = ( 3.0 - alpha ) * hElt / ( rhoElt * c2RElt ) / 6.0
    Me( 1, 2 ) =         alpha   * hElt / ( rhoElt * c2RElt ) / 6.0
    Me( 2, 2 ) = ( 3.0 - alpha ) * hElt / ( rhoElt * c2RElt ) / 6.0

    ! * A2 matrix *
    AD2( Node   ) = AD2( Node ) + &
         &          Me( 1, 1 ) + 0.5 * deltat * Ce( 1, 1 ) + beta * deltat2 * Ke( 1, 1 )
    AE2( Node+1 ) = Me( 1, 2 ) + 0.5 * deltat * Ce( 1, 2 ) + beta * deltat2 * Ke( 1, 2 )
    AD2( Node+1 ) = Me( 2, 2 ) + 0.5 * deltat * Ce( 2, 2 ) + beta * deltat2 * Ke( 2, 2 )

    ! * A1 matrix *
    AD1( Node   ) = AD1( Node ) + &
         &          2.0*Me( 1, 1 ) - ( 1.0 - 2.0*beta ) *deltat2 * Ke( 1, 1 )
    AE1( Node+1 ) = 2.0*Me( 1, 2 ) - ( 1.0 - 2.0*beta ) *deltat2 * Ke( 1, 2 )
    AD1( Node+1 ) = 2.0*Me( 2, 2 ) - ( 1.0 - 2.0*beta ) *deltat2 * Ke( 2, 2 )

    ! * A0 matrix *
    AD0( Node   ) = AD0( Node ) &
         &          -Me( 1, 1 ) + 0.5 * deltat * Ce( 1, 1 ) - beta * deltat2 * Ke( 1, 1 )
    AE0( Node+1 ) = -Me( 1, 2 ) + 0.5 * deltat * Ce( 1, 2 ) - beta * deltat2 * Ke( 1, 2 )
    AD0( Node+1 ) = -Me( 2, 2 ) + 0.5 * deltat * Ce( 2, 2 ) - beta * deltat2 * Ke( 2, 2 )

  END SUBROUTINE CONTRIB

  !**********************************************************************

  SUBROUTINE STEP( AD0, AE0, AD1, AE1, AD2, RV1, RV2, RV4, U0, U1, U2, time, rkT, fLoCut, fHiCut, IniFlag )

    !     Take a time step

    !     This is the inner loop.
    !     A significant speed-up can be obtained by using real arithmetic and eliminating the
    !     tridiagonal solver if the Hilbert transform is bypassed and the
    !     explicit solver is used.

    USE factor_Mod
    USE backsub_Mod
    USE sourceMod

    LOGICAL :: IniFlag
    REAL,    INTENT( IN   ) :: time, rkT, fLoCut, fHiCut
    COMPLEX, INTENT( IN   ) :: AD2( NTot1 ), AD1( NTot1 ), AE1( NTot1 ), AD0( NTot1 ), &
                               AE0( NTot1 ), RV1( NTot1 ), RV2( NTot1 ), RV4( NTot1 )
    COMPLEX, INTENT( INOUT) :: U0( NTot1 ), U1( NTot1 )
    COMPLEX, INTENT( OUT  ) :: U2( NTot1 )
    INTEGER                 :: is, J, js, L, medium, Ntpts
    REAL                    :: timeV( 1 ), deltat2, omega
    COMPLEX                 :: ST( Pos%NSz )
    CHARACTER (LEN=60)      :: PulseTitle

    omega   = SQRT( omega2 )
    deltat2 = deltat ** 2

    ! Form U2TEMP = A1 * U1 + A0 * U0 + S

    ! Surface point
    J = 1
    U2( J ) = AD1( J ) * U1( J ) + AE1( J + 1 ) * U1( J + 1 ) + AD0( J ) * U0( J ) + AE0( J + 1 ) * U0( J + 1 )

    ! Interior points
    L = NTot1 - 1
    U2( 2 : L ) = AD1( 2 : L ) * U1( 2 : L ) + AE1( 2 : L ) * U1( 1 : L - 1 ) + AE1( 3 : L + 1 ) * U1( 3 : L + 1 ) &
                + AD0( 2 : L ) * U0( 2 : L ) + AE0( 2 : L ) * U0( 1 : L - 1 ) + AE0( 3 : L + 1 ) * U0 (3 : L + 1 )

    ! Bottom point
    J = NTot1
    U2( J ) = AD1( J ) * U1( J ) + AE1( J ) * U1( J - 1 ) + AD0( J ) * U0( J ) + AE0( J ) * U0( J - 1 )

    ! Source terms
    Ntpts      = 1
    timeV( 1 ) = time   ! source expects a vector, not a scalar
    CALL SOURCE( timeV, ST, Pos%Sz, Pos%NSz, Ntpts, omega, fLoCut, fHiCut, Pulse, PulseTitle, IniFlag )

    medium = 1          ! Assumes source in first medium

    ST = CMPLX( deltat2 * ST * EXP( -i * rkT * V * time ) )

    SourceDepth: DO iS = 1, Pos%NSz
       js = Pos%iSz( iS ) + ( medium - 1 )

       U2( js     ) = U2( js     ) + ( 1.0 - Pos%ws( iS ) ) * ST( iS )
       U2( js + 1 ) = U2( js + 1 ) +         Pos%ws( iS )   * ST( iS )
    END DO SourceDepth

    ! Solve A2 * U2 = U2Temp

    ! Implicit or explicit solver?
    IF ( alpha == 0.0 .AND. beta == 0.0 ) THEN
       U2 = U2 / AD2
    ELSE
       CALL BackSub( NTot1, RV1, RV2, RV4, U2 )
    ENDIF

    ! Do a roll
    U0 = U1
    U1 = U2

    ! * Boundary conditions (natural is natural) *
    IF ( TopOpt( 2 : 2 ) == 'V' ) U1( 1 )     = 0.0
    IF ( BotOpt( 1 : 1 ) == 'V' ) U1( NTot1 ) = 0.0

  END SUBROUTINE STEP

  !**********************************************************************

  SUBROUTINE EXTRACT( U0, U1, ik, time, Itout, rkT )

    ! Extract solution (snapshot, vertical or horizontal array)

    INTEGER, INTENT( IN    ) :: ik                 ! index of wavenumber
    REAL,    INTENT( IN    ) :: time, rkT          ! Extract solution a time, wavenumber
    COMPLEX, INTENT( IN    ) :: U0( * ), U1( * )   ! Solution at two successive time steps
    INTEGER, INTENT( INOUT ) :: Itout
    REAL    :: T, wt
    COMPLEX :: UT1, UT2, U, const

    ! Three cases: snapshot, vertical time series, horizontal time series

    DO WHILE ( Itout <= Ntout .AND. time + deltat >= tout( Itout ) )
       ! note above could access tout( ntout + 1 ) which is outside array dimension
       ! but result is irrelevant
       ! exit statement below added to prevent that

       wt = ( tout( Itout ) - time ) / deltat  ! Weight for temporal interpolation:

       SELECT CASE ( TopOpt( 5 : 5 ) )
       CASE ( 'S' ) ! Case of a snapshot
          const = CMPLX( EXP( i * rkT * V * tout( Itout ) ) )

          DO ir = 1, Pos%NRz
             ! Linear interpolation in depth
             ii  = Pos%irz( ir )
             UT1 = U0( ii ) + Pos%wr( ir ) * ( U0( ii+1 ) - U0( ii ) )
             UT2 = U1( ii ) + Pos%wr( ir ) * ( U1( ii+1 ) - U1( ii ) )

             ! Linear interpolation in time
             Green( Itout, ir, ik ) = const * ( UT1 + wt * ( UT2 - UT1 ) )
          END DO

       CASE ( 'D' ) ! Case of an RTS (vertical array)
          ! const = -SQRT( 2.0 ) * deltak * SQRT( rkT ) *  COS(      rkT * r( 1 ) - pi / 4.0 )
          const = CMPLX( SQRT( 2.0 ) * deltak * SQRT( rkT ) * &
                  EXP( i * ( rkT * ( V * tout( Itout ) - Pos%Rr( 1 ) ) + SNGL( pi ) / 4.0 ) ) )

          DO ir = 1, Pos%NRz
             ! Linear interpolation in depth
             ii  = Pos%irz( ir )
             UT1 = U0( ii ) + Pos%wr( ir ) * ( U0( ii+1 ) - U0( ii ) )
             UT2 = U1( ii ) + Pos%wr( ir ) * ( U1( ii+1 ) - U1( ii ) )

             ! Linear interpolation in time
             U = UT1 + wt * ( UT2 - UT1 )
             RTSrz( ir, Itout ) = RTSrz( ir, Itout ) + REAL( const * U )
          END DO

       CASE ( 'R' ) ! Case of an RTS (horizontal array)
          IF ( time >= tout( 1 ) ) THEN
             ! Linear interpolation in depth
             ii  = Pos%irz( 1 )
             T   = ( Pos%rz( 1 ) - z( ii ) ) / ( z( ii + 1 ) - z( ii ) )

             UT1 = U0( ii ) + T * ( U0( ii + 1 ) - U0( ii ) )
             UT2 = U1( ii ) + T * ( U1( ii + 1 ) - U1( ii ) )

             ! Linear interpolation in time
             U = UT1 + wt * ( UT2 - UT1 )
             const = CMPLX( SQRT( 2.0 ) * EXP( i * rkT * V * time ) )

             RTSrr( 1 : Pos%NRr, Itout ) = SNGL( RTSrr( 1 : Pos%NRr, Itout ) + REAL( deltak * U *  &
                  const * EXP( i * ( -rkT * Pos%Rr( 1 : Pos%NRr ) + SNGL( pi ) / 4.0 ) ) ) * SQRT( rkT / Pos%Rr( 1 : Pos%NRr ) ) )

          ENDIF
       END SELECT
       Itout = Itout + 1
       IF ( Itout == Ntout + 1 ) EXIT   ! did last time point, so exit the do loop
    END DO

  END SUBROUTINE EXTRACT

END PROGRAM SPARC
