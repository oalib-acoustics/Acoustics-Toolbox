PROGRAM KRAKEN

  ! Ocean acoustic normal modes.

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

  ! Originally developed as part of the author's dissertation under the supervision
  ! of Prof. Edward L. Reiss, Northwestern University

  USE KrakenMod
  USE ReadEnvironmentMod
  USE SourceReceiverPositions
  USE RefCoef
  USE sspMod
  USE AttenMod
  USE FatalError

  IMPLICIT NONE
  INTEGER            :: Min_Loc( 1 ), iFirst, iLast, iRec, ii, j
  REAL               :: zMin, zMax
  REAL      (KIND=8) :: Error, freq0
  CHARACTER (LEN=80) :: FileRoot

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  Profile: DO iProf = 1, 9999   ! Loop over a sequence of profiles
     NV( 1 : 5 ) = [ 1, 2, 4, 8, 16 ]

     Title = 'KRAKEN- '
     CALL ReadEnvironment( FileRoot, Title, freq0, MaxMedium, TopOpt, NG, BotOpt, cLow, cHigh, RMax, ENVFile, PRTFile )
     IF ( HSBot%BC == 'F' .OR. HSTop%BC == 'P' ) &
        CALL ERROUT( 'KRAKEN', 'The option to read a file for the reflection loss is not implemented in KRAKEN' )
     zMin = SNGL( SSP%Depth( 1              ) )
     zMax = SNGL( SSP%Depth( SSP%NMedia + 1 ) )
     CALL ReadSzRz(    zMin, zMax )   ! Read source/receiver depths
     CALL ReadfreqVec( freq0, TopOpt( 6 : 6 ) )

     IF ( iProf == 1 ) CALL ReadReflectionCoefficient( FileRoot, HSBot%BC, HSTop%BC, PRTFile )

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

        WRITE( PRTFile, * )
        WRITE( PRTFile, * ) 'Mesh multiplier   CPU seconds'

        DO iSet = 1, NSets   ! Main loop: solve the problem for a sequence of meshes
           N( 1 : SSP%NMedia ) = INT( NG( 1 : SSP%NMedia ) * NV( iSet ) * freq / freq0 )   ! scaled by frequency
           h( 1 : SSP%NMedia ) = ( SSP%Depth( 2 : SSP%NMedia + 1 ) - SSP%Depth( 1 : SSP%NMedia ) ) / N( 1 : SSP%NMedia )
           hV( iSet )      = h( 1 )
           CALL Solve

           IF ( Error * 1000.0 * RMax < 1.0 ) THEN   ! Convergence achieved
              EXIT
           ELSE
              IF ( iSet == NSets ) WRITE( PRTFile, * ) 'Warning in KRAKEN : Too many meshes needed: check convergence'
           END IF
        END DO

        ! Solution complete: discard modes with phase velocity above cHigh

        Min_Loc = MINLOC( Extrap( 1, 1 : M ), Extrap( 1, 1 : M ) > omega2 / cHigh ** 2 )
        M       = Min_Loc( 1 )

        ! Write eigenvalues to PRTFile and MODFile

        WRITE( PRTFile, * )
        WRITE( PRTFile, * ) '   I       k             alpha     Phase Speed     Group Speed'
        WRITE( PRTFile, * ) '         (1/m)           (1/m)        (m/s)           (m/s)'

        ! k() contains scatter losses; Extrap contains extrapolated values of wavenumbers, so combine ...
        k( 1 : M ) = SQRT( Extrap( 1, 1 : M ) + k( 1 : M ) )

        DO mode = 1, M, MAX( 1, M / 30 )  ! print every mode unless there are an awful lot
           WRITE( PRTFile, '( I5, G18.10, G10.2, G18.10, G14.6  )' ) mode, k( mode ), omega / DBLE( k( mode ) ), VG( mode )
        END DO
        FLUSH( PRTFile )

        WRITE( MODFile, REC = iRecProfile ) M

        iFirst = 1
        DO IREC = 1, 1 + ( 2 * M - 1 ) / LRecordLength
           iLast  = MIN( M, iFirst + LRecordLength / 2 - 1 )
           WRITE( MODFile, REC = iRecProfile + 1 + M + IREC ) CMPLX( k( iFirst : iLast ) )
           iFirst = iLast + 1
        END DO
        FLUSH( MODFile )

        ! set record pointer to beginning of next mode set
        iRecProfile = iRecProfile + 3 + M + ( 2 * M - 1 ) / LRecordLength
     END DO FreqLoop
  END DO Profile

  CLOSE( ENVFile )
  CLOSE( PRTFile )
  CLOSE( MODFile )

CONTAINS

  !**********************************************************************

  SUBROUTINE Initialize

    ! Initializes arrays defining difference equations

    LOGICAL           :: ElasticFlag = .FALSE.
    INTEGER           :: IAllocStat, Medium, NPoints, N1
    REAL     (KIND=8) :: Two_h
    REAL     (KIND=8) :: cP2, cS2
    COMPLEX  (KIND=8), ALLOCATABLE :: cP( : ), cS( : )
    CHARACTER (LEN=8) :: Task

    cMin          = HUGE( cMin )
    FirstAcoustic = 0
    Loc( 1 )      = 0

    ! Allocate storage for finite-difference coefficients

    NPoints = SUM( N( 1 : SSP%NMedia ) ) + SSP%NMedia

    IF ( ALLOCATED( B1 ) ) DEALLOCATE( B1, B1C, B2, B3, B4, rho )
    ALLOCATE ( B1( NPoints ), B1C( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints ), rho( NPoints ), &
         cP( NPoints ), cS( NPoints ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) &
         CALL ERROUT( 'KRAKEN - Initialize', 'Insufficient memory to allocate B1, B2, B3, B4 vectors: Reduce mesh.' )

    Media: DO Medium = 1, SSP%NMedia

       IF ( Medium /= 1 ) Loc( Medium ) = Loc( Medium - 1 ) + N( Medium - 1 ) + 1
       N1 = N(   Medium ) + 1
       ii = Loc( Medium ) + 1

       ! EvaluateSSP reads in the data for a medium
       Task = 'TAB'
       CALL EvaluateSSP( cP( ii ), cS( ii ), rho( ii ), Medium, N1, freq, Task )

       ! Load diagonals of the finite-difference equations

       IF ( REAL( cS( ii ) ) == 0.0 ) THEN   ! Case of an acoustic medium

          SSP%Material( Medium ) = 'ACOUSTIC'
          IF ( FirstAcoustic == 0 ) FirstAcoustic = Medium
          LastAcoustic = Medium

          cMin = MIN( cMin, MINVAL( DBLE( cP( ii : ii + N( Medium ) ) ) ) )

          B1(  ii : ii + N( Medium ) ) = -2.0D0 + h( Medium ) ** 2 * DBLE( omega2 / cP( ii : ii + N( Medium ) ) ** 2 )
          B1C( ii : ii + N( Medium ) ) =                            AIMAG( omega2 / cP( ii : ii + N( Medium ) ) ** 2 )

       ELSE                                  ! Case of an elastic medium
          IF ( SSP%sigma( Medium ) /= 0.0D0 ) CALL ERROUT( 'KRAKEN', 'Rough elastic interfaces are not allowed' )

          SSP%Material( Medium ) = 'ELASTIC'
          ElasticFlag        = .TRUE.
          Two_h              = 2.0D0 * h( Medium )

          DO j = ii, ii + N( Medium )
             cMin = MIN( DBLE( cS( j ) ), cMin )

             cP2 = DBLE( cP( j ) ** 2 )
             cS2 = DBLE( cS( j ) ** 2 )

             B1(  j ) = Two_h / ( rho( j ) * cS2 )
             B2(  j ) = Two_h / ( rho( j ) * cP2 )
             B3(  j ) = 4.0D0 * Two_h * rho( j ) * cS2 * ( cP2 - cS2 ) / cP2
             B4(  j ) = Two_h * ( cP2 - 2.0D0 * cS2 ) / cP2
             rho( j ) = Two_h * omega2 * rho( j )
          END DO

       END IF

    END DO Media

    ! (cLow, cHigh) = phase speed interval for the mode search
    ! user specified interval is reduced if it exceeds domain of possible mode phase speeds

    ! Bottom properties
    IF ( HSBot%BC == 'A' ) THEN
       IF ( REAL( HSBot%cS ) > 0.0 ) THEN      ! Elastic  bottom half-space
          ElasticFlag = .TRUE.
          cMin  = MIN( cMin,  DBLE( HSBot%cS ) )
          cHigh = MIN( cHigh, DBLE( HSBot%cS ) )
       ELSE                                    ! Acoustic bottom half-space
          cMin  = MIN( cMin,  DBLE( HSBot%cP ) )
          ! cHigh = MIN( cHigh, DBLE( HSBot%cP ) )
       END IF
    END IF

    ! Top properties
    IF ( HSTop%BC == 'A' ) THEN
       IF ( REAL( HSTop%cS ) > 0.0 ) THEN       ! Elastic  top half-space
          ElasticFlag = .TRUE.
          cMin   = MIN( cMin,  DBLE( HSTop%cS ) )
          cHigh  = MIN( cHigh, DBLE( HSTop%cS ) )
       ELSE                                     ! Acoustic top half-space
          cMin   = MIN( cMin,  DBLE( HSTop%cP ) )
          ! cHigh  = MIN( cHigh, DBLE( HSTop%cP ) )
       END IF
    END IF

    ! If elastic medium then reduce cMin for Scholte wave
    IF ( ElasticFlag ) cMin = 0.85 * cMin
    cLow = MAX( cLow, cMin )

  END SUBROUTINE Initialize

  !**********************************************************************

  SUBROUTINE Solve

    ! Solves the eigenvalue problem at the current mesh and produces a new extrapolation

    INTEGER          :: Min_Loc( 1 ), Key
    INTEGER          :: NTotal, NTotal1   ! number of mesh points (where eigenvector is sampled)
    REAL    (KIND=8) :: x1, x2, TStart, TEnd
    REAL    (KIND=8) :: T1, T2, F1, F2

    CALL CPU_TIME( Tstart )
    CALL Initialize           ! set up the finite-difference mesh

    ! Choose a solver ...
    IF ( iProf > 1 .AND. iSet <= 2 .AND. TopOpt( 4 : 4 ) == 'C' ) THEN
       CALL Solve3 ! continuation from last profile if option selected and doing the first or second mesh
    ELSE IF ( ( iSet <= 2 ) .AND. ( SSP%NMedia <= LastAcoustic - FirstAcoustic + 1 ) ) THEN
       CALL Solve1 ! bisection for first two sets (applicable if elasticity is limited to halfspaces)
    ELSE
       CALL Solve2 ! extrapolation from first two meshes
    END IF

    Extrap( iSet, 1 : M ) = EVMat( iSet, 1 : M )

    ! Remove any eigenvalues outside the spectral limits
    ! Typically the last 'eigenvalue' results from forcing a zero in funct( x ) for x outside the limits
    ! Inverse iteration would usually fail for that mode

    Min_Loc = MINLOC( Extrap( 1, 1 : M ), Extrap( 1, 1 : M ) > omega2 / cHigh ** 2 )
    M       = Min_Loc( 1 )

    NTotal  = SUM( N( FirstAcoustic : LastAcoustic ) )
    NTotal1 = NTotal + 1

    IF ( iSet == 1 ) CALL Vector( NTotal, NTotal1 )   ! If this is the first mesh, compute the eigenvectors

    ! Richardson extrapolation to improve the accuracy

    Error = 1.0D10          ! initialize error to a large number
    KEY   = 2 * M / 3 + 1   ! index of element used to check convergence

    IF ( iSet > 1 ) THEN
       T1 = Extrap( 1, KEY )

       DO j = iSet - 1, 1, -1
          ModeLoop: DO mode = 1, M
             x1 = NV( j    ) ** 2
             x2 = NV( iSet ) ** 2
             F1 = Extrap( j,     mode )
             F2 = Extrap( j + 1, mode )
             Extrap( j, mode ) = F2 - ( F1 - F2 ) / ( x2 / x1 - 1.0D0 )
          END DO ModeLoop
       END DO

       T2    = Extrap( 1, KEY )
       Error = ABS( T2 - T1 )
    END IF

    CALL CPU_TIME( Tend )   ! check elapsed time
    ET( iSet ) = Tend - Tstart
    WRITE( PRTFile, '( 1X, I8, 6X, G15.3, ''s'' )' ) NV( iSet ), ET( iSet )

  END SUBROUTINE Solve

  !**********************************************************************

  SUBROUTINE Solve2

    ! Provides initial guess to root finder for each EVMat( I )

    USE RootFinderSecantMod
    INTEGER            :: IAllocStat
    INTEGER            :: Iteration, MaxIteration = 2000   ! iterations in root finder
    REAL      (KIND=8) :: P( 10 ), x, x1, x2, Tolerance
    CHARACTER (LEN=80) :: ErrorMessage

    CountModes = .FALSE.

    x = omega2 / cLow ** 2

    ! solve1 has already allocated space for the following unless the problem has shear
    IF ( .NOT. ALLOCATED( k ) ) THEN
       M = 3000   ! this sets the upper limit on how many modes can be calculated
       ALLOCATE( EVMat( NSets, M ), Extrap( NSets, M ), k( M ), VG( M ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
            CALL ERROUT( 'KRAKEN - Solve2', 'Insufficient memory (too many modes).' )
    END IF

    ! start looking for each root in succession
    ModeLoop: DO mode = 1, M

       ! For first or second meshes, use a high guess
       ! Otherwise use extrapolation to produce an initial guess

       x = 1.00001 * x

       IF ( iSet >= 2 ) THEN

          P( 1 : iSet - 1 ) = EVMat( 1 : iSet - 1, mode )

          IF ( iSet >= 3 ) THEN
             DO ii = 1, iSet - 2
                DO j = 1, iSet - ii - 1
                   x1 = hV( j      ) ** 2
                   x2 = hV( j + ii ) ** 2

                   P( j ) = ( ( hV( iSet ) ** 2 - x2 ) * P( j     ) - &
                              ( hV( iSet ) ** 2 - x1 ) * P( j + 1 ) ) &
                        / ( x1 - x2 )
                END DO
             END DO
             x = P( 1 )
          END IF

       END IF

       ! Use the secant method to refine the eigenvalue
       ! With large matrices one needs to allow a larger tolerance 
       Tolerance = ABS( x ) * SIZE( B1 ) * 10.0D0 ** ( 1.0 - PRECISION( x ) )

       CALL RootFinderSecant( x, Tolerance, Iteration, MaxIteration, ErrorMessage, Funct )

       IF ( ErrorMessage /= ' ' ) THEN   ! any problems in root finder?
          WRITE( PRTFile, * ) 'iSet, mode = ', iSet, mode
          WRITE( PRTFile, * ) 'Warning in KRAKEN - RootFinderSecant', ErrorMessage
          x = TINY( x )   ! make sure value is discarded
       END IF

       EVMat( iSet, mode ) = x

       ! Toss out modes outside user specified spectrum
       IF ( omega2 / cHigh ** 2 > x ) THEN
          M = mode - 1
          RETURN
       END IF

    END DO ModeLoop

  END SUBROUTINE Solve2

  !**********************************************************************

  SUBROUTINE Solve3

    ! Provides initial guess to root finder for each EVMat(I)
    ! This solver tries to use eigenvalues from a previous profile as initial guesses

    USE RootFinderSecantMod
    INTEGER            :: It, MaxIt, iPower
    REAL      (KIND=8) :: x, xMin, Tolerance, Delta
    CHARACTER (LEN=80) :: ErrorMessage

    CountModes = .FALSE.

    MaxIT = 500

    ! Determine number of modes

    xMin = 1.00001D0 * omega2 / cHigh ** 2

    CALL FUNCT( xMin, Delta, iPower )
    M = modeCount

    ModeLoop: DO mode = 1, M

       x         = EVMat( iSet, mode )
       Tolerance = ABS( x ) * 10.0D0 ** ( 2.0 - PRECISION( x ) )

       CALL RootFinderSecant( x, Tolerance, IT, MaxIT, ErrorMessage, Funct )  ! Use the secant method to refine the eigenvalue

       IF ( ErrorMessage /= ' ' ) THEN
          WRITE( PRTFile, * ) 'iSet, mode = ', iSet, mode
          WRITE( PRTFile, * ) 'Warning in KRAKEN - RootFinderSecant', ErrorMessage
          x = TINY( x )   ! make sure value is discarded
       END IF

       EVMat( iSet, mode ) = x

       IF ( omega2 / cHigh ** 2 > x ) THEN  ! Toss out modes outside user specified spectrum
          M = mode - 1
          RETURN
       END IF

    END DO ModeLoop

  END SUBROUTINE Solve3

  !**********************************************************************

  SUBROUTINE FUNCT( x, Delta, iPower )

    ! FUNCT( x ) = 0 is the dispersion relation

    USE BCImpedanceMod

    INTEGER,          PARAMETER  :: iPowerR = 50, iPowerF = -50
    REAL    (KIND=8), PARAMETER  :: Roof = 1.0D+50, Floor = 1.0D-50
    INTEGER,          INTENT( INOUT ) :: iPower
    REAL    (KIND=8), INTENT( IN    ) :: x
    REAL    (KIND=8), INTENT( OUT   ) :: Delta
    INTEGER                           :: iPowerBot
    REAL    (KIND=8)                  :: f, g
    COMPLEX (KIND=8)                  :: fTop, gTop, fBot, gBot

!!$  IF ( x <= omega2 / cHigh ** 2 ) THEN    ! For a k below the cts spectrum limit, force a zero
!!$     Delta  = 0.0D0
!!$     iPower = 0
!!$     RETURN
!!$  END IF

    modeCount = 0

    CALL BCImpedance(    x, 'BOT', HSBot, fTop, gTop, iPower,    .FALSE. ) ! Bottom impedance
    f = DBLE( fTop )
    g = DBLE( gTop )
    CALL AcousticLayers( x,               f,    g,    iPower             ) ! Shoot through acoustic layers
    CALL BCImpedance(    x, 'TOP', HSTop, fBot, gBot, iPowerBot, .FALSE. ) ! Top impedance

    Delta  = REAL( f * gBot - g * fBot )
    iPower = iPower + iPowerBot
    IF ( g * Delta > 0.0D0 ) modeCount = modeCount + 1

    ! Deflate previous roots

    IF ( ( mode > 1 ) .AND. ( SSP%NMedia > LastAcoustic - FirstAcoustic + 1 ) ) THEN
       ModeLoop: DO j = 1, mode - 1
          Delta = Delta / ( x - EVMat( iSet, j ) )

          ! Scale if necessary
          DO WHILE ( ABS( Delta ) < Floor .AND. ABS( Delta ) > 0.0D0 )
             Delta  = Roof * Delta
             iPower = iPower - iPowerR
          END DO

          DO WHILE ( ABS( Delta ) > Roof )
             Delta  = Floor * Delta
             iPower = iPower - iPowerF
          END DO

       END DO ModeLoop
    END IF

  END SUBROUTINE FUNCT

  !**********************************************************************

  SUBROUTINE AcousticLayers( x, f, g, iPower )

    ! Shoot through acoustic layers

    INTEGER,          PARAMETER       :: iPowerF = -50
    REAL    (KIND=8), PARAMETER       :: Roof = 1.0D+50, Floor = 1.0D-50
    INTEGER,          INTENT( INOUT ) :: iPower
    REAL    (KIND=8), INTENT( IN    ) :: x
    REAL    (KIND=8), INTENT( INOUT ) :: f, g
    INTEGER                           :: Medium
    REAL    (KIND=8)                  :: rhoMedium
    REAL    (KIND=8)                  :: p0 = 0, p1, p2, h2k2

    IF ( FirstAcoustic == 0 ) RETURN

    Media: DO Medium = LastAcoustic, FirstAcoustic, -1   ! Loop over successive acoustic media

       h2k2      = h( Medium ) ** 2 * x
       ii        = Loc( Medium ) + N( Medium ) + 1
       rhoMedium = rho(  Loc( Medium ) + 1  )   ! density is made homogeneous using value at top of each medium

       p1 = -2.0D0 * g
       p2 = ( B1( ii ) - h2k2 ) * g - 2.0D0 * h( Medium ) * f * rhoMedium

       ! Shoot (towards surface) through a single medium
       DO ii = Loc( Medium ) + N( Medium ), Loc( Medium ) + 1, -1

          p0 = p1
          p1 = p2
          p2 = ( h2k2 - B1( ii ) ) * p1 - p0

          IF ( CountModes ) THEN
             IF ( p0 * p1 <= 0.0D0 ) modeCount = modeCount + 1
          END IF

          IF ( ABS( p2 ) > Roof ) THEN   ! Scale if necessary
             p0     = Floor * p0
             p1     = Floor * p1
             p2     = Floor * p2
             iPower = iPower - iPowerF
          END IF

       END DO

       ! f = P' / rho and g = -P since f P + g P' / rho = 0
       rhoMedium = rho( Loc( Medium ) + 1 )   ! density at top of layer
       f         = -( p2 - p0 ) / ( 2.0D0 * h( Medium ) ) / rhoMedium
       g         = -p1
    END DO Media

  END SUBROUTINE AcousticLayers

  !**********************************************************************

  SUBROUTINE Vector( NTotal, NTotal1 )

    ! Do inverse iteration to compute each of the eigenvectors and write these to the disk file

    USE BCImpedanceMod
    USE InverseIterationMod
    USE calculateweights
    USE MergeVectorsMod

    INTEGER, INTENT( IN ) :: NTotal, NTotal1   ! number of mesh points (where eigenvector is sampled)
    INTEGER              :: iError, iPower, iTurningPoint, jj, L, Medium, NzTab
    INTEGER, ALLOCATABLE :: IzTab( : )
    REAL                 :: zTab( Pos%NSz + Pos%NRz ), z( NTotal1 )
    REAL,    ALLOCATABLE :: WTS( : )
    REAL    (KIND=8)     :: h_rho = 0
    REAL    (KIND=8)     :: x, xh2
    REAL    (KIND=8)     :: Phi( NTotal1 ), d( NTotal1 ), e( NTotal1 + 1 )
    COMPLEX, ALLOCATABLE :: PhiTab( : )
    COMPLEX (KIND=8)     :: fTop, gTop, fBot, gBot

    ! Tabulate z-coordinates and off-diagonals of matrix

    j      = 1
    z( 1 ) = SNGL( SSP%Depth( FirstAcoustic ) )

    DO Medium = FirstAcoustic, LastAcoustic

       h_rho = h( Medium ) * rho( Loc( Medium ) + 1 )        ! density at the top of each layer

       e( j + 1 : j + N( Medium ) ) = 1.0D0 / h_rho
       z( j + 1 : j + N( Medium ) ) = z( j ) + SNGL( h( Medium ) * [ (jj, jj = 1, N( Medium ) ) ] )

       j = j + N( Medium )
    END DO

    e( NTotal1 + 1 ) = 1.0D0 / h_rho       ! Dummy value; never used

    ! Calculate the indices, weights, ... for mode interpolation
    CALL MergeVectors( Pos%Sz, Pos%NSz, Pos%Rz, Pos%NRz, zTab, NzTab )

    ! following lines would merge vectors without sorting and removing duplicates
    ! NzTab = Pos%NSz + Pos%NRz
    ! zTab( 1           : Pos%NSz ) = Pos%sz( 1 : Pos%NSz )
    ! zTab( Pos%NSz + 1 : NzTab   ) = Pos%rz( 1 : Pos%NRz )

    ALLOCATE( WTS( NzTab ), IzTab( NzTab ), PhiTab( NzTab ) )
    CALL Weight( z, NTotal1, zTab, NzTab, WTS, IzTab )

    ! Open MODFile and write header

    IF ( ifreq == 1 .AND. iProf == 1 ) THEN
       ! LRecordLength must not increase between profiles !!!
       LRecordLength = MAX( 2 * Nfreq, 2 * NzTab, 32, 3 * ( LastAcoustic - FirstAcoustic + 1 ) )   ! Logical record length in `longwords' (4 bytes)
       OPEN ( FILE = TRIM( FileRoot) //'.mod', UNIT = MODFile, ACCESS = 'DIRECT', STATUS = 'REPLACE', &
              RECL = 4 * LRecordLength, FORM = 'UNFORMATTED' )
    END IF

    IF ( ifreq == 1 ) THEN
       WRITE( MODFile, REC = iRecProfile     ) LRecordLength, Title, Nfreq, LastAcoustic - FirstAcoustic + 1, NzTab, NzTab
       WRITE( MODFile, REC = iRecProfile + 1 ) ( N( Medium ), SSP%Material( Medium ), Medium = FirstAcoustic, LastAcoustic )
       WRITE( MODFile, REC = iRecProfile + 2 ) ( REAL( SSP%Depth( Medium ) ), REAL( rho( Loc( Medium ) + 1 ) ), &
            Medium = FirstAcoustic, LastAcoustic )
       WRITE( MODFile, REC = iRecProfile + 3 ) freqVec( 1 : Nfreq )
       WRITE( MODFile, REC = iRecProfile + 4 ) zTab( 1 : NzTab )
       iRecProfile = iRecProfile + 5
    END IF

    ! top and bottom halfspace info (changes with frequency and possibly profile)
    WRITE( MODFile, REC = iRecProfile + 1 ) &
         HSTop%BC, CMPLX( HSTop%cP ), CMPLX( HSTop%cS ), REAL( HSTop%rho ), REAL( SSP%Depth( 1              ) ), &
         HSBot%BC, CMPLX( HSBot%cP ), CMPLX( HSBot%cS ), REAL( HSBot%rho ), REAL( SSP%Depth( SSP%NMedia + 1 ) )

    ! Main loop: for each eigenvalue call InverseIteration to get eigenvector

    ModeLoop: DO mode = 1, M
       x = EVMat( 1, mode )

       ! Corner elt requires top impedance
       CALL BCImpedance( x, 'TOP', HSTop, fTop, gTop, iPower, .FALSE. )

       IF ( gTop == 0.0D0 ) THEN
          d( 1 ) = 1.0D0
          e( 2 ) = 0.0D0
       ELSE
          L      = Loc( FirstAcoustic ) + 1
          xh2    = x * h( FirstAcoustic ) * h( FirstAcoustic )
          h_rho  = h( FirstAcoustic ) * rho( L )
          d( 1 ) = ( B1( L ) - xh2 ) / h_rho / 2.0D0 + REAL( fTop / gTop )
       END IF

       ! Set up the diagonal
       ! iTurningPoint was intended to be an index of the near surface turning point for bowl shaped profiles
       ! It is actually the index closest to the surface where the mode is still oscillatory
       iTurningPoint = NTotal
       j   = 1
       L   = Loc( FirstAcoustic ) + 1

       Media: DO Medium = FirstAcoustic, LastAcoustic
          xh2   = x * h( Medium ) ** 2
          h_rho = h( Medium ) * rho( Loc( Medium ) + 1 )

          IF ( Medium >= FirstAcoustic + 1 ) THEN
             L      = L + 1
             d( j ) = ( d( j ) + ( B1( L ) - xh2 ) / h_rho ) / 2.0D0
          END IF

          DO ii = 1, N( Medium )
             j      = j + 1
             L      = L + 1
             d( j ) = ( B1( L ) - xh2 ) / h_rho

             IF ( B1( L ) - xh2 + 2.0D0 > 0.0D0 ) THEN   ! Find index of turning point nearest top
                iTurningPoint = MIN( j, iTurningPoint )
             END IF
          END DO

       END DO Media

       ! Corner elt requires bottom impedance
       CALL BCImpedance( x, 'BOT', HSBot, fBot, gBot, iPower, .FALSE. )

       IF ( gBot == 0.0D0 ) THEN
          d( NTotal1 ) = 1.0D0
          e( NTotal1 ) = 0.0D0
       ELSE
          d( NTotal1 ) = d( NTotal1 ) / 2.0D0 - REAL( fBot / gBot )
       END IF

       CALL InverseIteration( NTotal1, d, e, iError, Phi )   ! Inverse iteration to compute eigenvector

       IF ( iError /= 0 ) THEN
          WRITE( PRTFile, * ) 'mode = ', mode
          WRITE( PRTFile, * ) 'Warning in KRAKEN - InverseIteration: Inverse iteration failed to converge'
          Phi = 0.0   ! zero out the errant eigenvector
       ELSE
          CALL Normalize( Phi, iTurningPoint, NTotal1, x )   ! Normalize the eigenvector
       END IF

       ! Tabulate the modes at the source/rcvr depths and write to disk
       PhiTab = CMPLX( Phi( IzTab ), KIND=4 ) + WTS * CMPLX( Phi( IzTab + 1 ) - Phi( IzTab ), KIND=4 )
       WRITE( MODFile, REC = iRecProfile + 1 + mode ) PhiTab

    END DO ModeLoop

    DEALLOCATE( WTS, IzTab, PhiTab )

  END SUBROUTINE Vector

  !**********************************************************************

  SUBROUTINE Normalize( Phi, iTurningPoint, NTotal1, x )

    ! Normalize the eigenvector:
    ! SqNorm = Integral(Phi ** 2) by the trapezoidal rule:
    ! Integral( F ) = H * ( F(1) + ... + F(N-1) + 0.5 * ( F(0) + F(N) ) )

    ! Compute perturbation due to material absorption
    ! Compute the group velocity
    ! Call ScatterLoss to figure interfacial scatter loss

    USE BCImpedanceMod
    INTEGER,          INTENT( IN    ) :: iTurningPoint, NTotal1     ! index of turning point, number of points in mesh
    REAL    (KIND=8), INTENT( IN    ) :: x                ! eigenvalue
    REAL    (KIND=8), INTENT( INOUT ) :: Phi( NTotal1 )   ! eigenvector
    INTEGER                           :: iPower, j1, L, L1, Medium
    REAL    (KIND=8)                  :: rhoMedium, rho_omega_h2
    REAL    (KIND=8)                  :: x1, x2, DrhoDx, DetaDx, SqNorm, RN, ScaleFactor, Slow
    COMPLEX (KIND=8)                  :: Perturbation_k, Del
    COMPLEX (KIND=8)                  :: fTop1, gTop1, fTop2, gTop2, fBot1, gBot1, fBot2, gBot2

    SqNorm = 0.0D0
    Perturbation_k   = 0.0D0
    Slow   = 0.0D0

    ! Compute contribution from the top half-space
    SELECT CASE( TopOpt( 2 : 2 ) )
    CASE ( 'A' )
       !Del  = -0.5D0 * ( omega2 / HSTop%cP ** 2 - DBLE( omega2 / HSTop%cP ** 2 ) ) / &
       !            SQRT( x                      - DBLE( omega2 / HSTop%cP ** 2 ) )
       Del  = i * AIMAG( SQRT( CMPLX( x - omega2 / HSTop%cP ** 2 ) ) )
       Perturbation_k = Perturbation_k - Del * Phi( 1 ) ** 2 / HSTop%rho
       Slow = Slow + Phi( 1 ) ** 2 / ( 2 * SQRT( x - DBLE( omega2 / HSTop%cP ** 2 ) ) ) / ( HSTop%rho * DBLE( HSTop%cP ) ** 2 )
    CASE ( 'F', 'P' )
       CALL BCImpedance( x, 'TOP', HSTop, fTop1, gTop1, iPower, .FALSE. )
       CALL BCImpedance( x, 'TOP', HSTop, fTop2, gTop2, iPower, .TRUE.  )
       Del  = fTop2 / gTop2 - fTop1 / gTop1
       Perturbation_k = Perturbation_k - Del * Phi( 1 ) ** 2
    END SELECT

    ! Compute contribution from the volume
    L = Loc( FirstAcoustic )
    j = 1

    Media: DO Medium = FirstAcoustic, LastAcoustic
       L            = L + 1
       rhoMedium    = rho( L )
       rho_omega_h2 = rhoMedium * omega2 * h( Medium ) ** 2

       ! top interface
       SqNorm = SqNorm + 0.5D0 * h( Medium ) *                      Phi( j ) ** 2 / rhoMedium
       Slow   = Slow   + 0.5D0 * h( Medium ) * ( B1( L ) + 2.D0 ) * Phi( j ) ** 2 / rho_omega_h2
       Perturbation_k   = Perturbation_k   + 0.5D0 * h( Medium ) * i * B1C( L )       * Phi( j ) ** 2 / rhoMedium

       ! medium
       L1 = L + 1
       L  = L + N( Medium ) - 1
       j1 = j + 1
       j  = j + N( Medium ) - 1

       SqNorm = SqNorm + h( Medium ) * SUM(                           Phi( j1 : j ) ** 2 ) / rhoMedium
       Slow   = Slow   + h( Medium ) * SUM( ( B1( L1 : L ) + 2.D0 ) * Phi( j1 : j ) ** 2 ) / rho_omega_h2
       Perturbation_k   = Perturbation_k   + h( Medium ) * i * SUM(  B1C( L1 : L )          * Phi( j1 : j ) ** 2 ) / rhoMedium

       ! bottom interface
       L = L + 1
       j = j + 1

       SqNorm = SqNorm + 0.5D0 * h( Medium ) *                      Phi( j ) ** 2 / rhoMedium
       Slow   = Slow   + 0.5D0 * h( Medium ) * ( B1( L ) + 2.D0 ) * Phi( j ) ** 2 / rho_omega_h2
       Perturbation_k   = Perturbation_k   + 0.5D0 * h( Medium ) * i * B1C( L )        * Phi( j ) ** 2 / rhoMedium

    END DO Media

    ! Compute contribution from the bottom half-space
    SELECT CASE( BotOpt( 1 : 1 ) )
    CASE ( 'A' )
       ! This is the original perturbation method
       !   Del  = i * AIMAG( SQRT( CMPLX( x - omega2 / HSBot%cP ** 2 ) ) )
       !   Perturbation_k = Perturbation_k - Del * Phi( j ) ** 2 / HSBot%rho
       !   Slow = Slow + Phi( j ) ** 2 / ( 2 * SQRT( x - DBLE( omega2 / HSBot%cP ** 2 ) ) ) / ( HSBot%rho * DBLE( HSBot%cP ) ** 2 )

       ! This formula uses a loop length and works for leaky modes
       CALL BCImpedance( x, 'BOT', HSBot, fBot1, gBot1, iPower, .FALSE. )
       CALL BCImpedance( x, 'BOT', HSBot, fBot2, gBot2, iPower, .TRUE.  )
       Del  = fBot2 / gBot2 - fBot1 / gBot1
       Perturbation_k = Perturbation_k - Del * Phi( j ) ** 2

       ! group velocity still uses the perturbation method so group speeds will be wrong for leaky modes
       Slow = Slow + Phi( j ) ** 2 / ( 2 * SQRT( x - DBLE( omega2 / HSBot%cP ** 2 ) ) ) / ( HSBot%rho * DBLE( HSBot%cP ) ** 2 )
    CASE ( 'F', 'P' )
       CALL BCImpedance( x, 'BOT', HSBot, fBot1, gBot1, iPower, .FALSE. )
       CALL BCImpedance( x, 'BOT', HSBot, fBot2, gBot2, iPower, .TRUE.  )
       Del  = fBot2 / gBot2 - fBot1 / gBot1
       Perturbation_k = Perturbation_k - Del * Phi( j ) ** 2
    END SELECT

    ! Compute derivative of top admitance
    x1 = 0.9999999D0 * x
    x2 = 1.0000001D0 * x

    CALL BCImpedance( x1, 'TOP', HSTop, fTop1, gTop1, iPower, .FALSE. )
    CALL BCImpedance( x2, 'TOP', HSTop, fTop2, gTop2, iPower, .FALSE. )

    DrhoDx = 0.0D0
    IF ( gTop1 /= 0.0D0 ) DrhoDx = REAL( ( fTop2 / gTop2 - fTop1 / gTop1 ) ) / ( x2 - x1 )

    ! Compute derivative of bottom admitance

    CALL BCImpedance( x1, 'BOT', HSBot, fBot1, gBot1, iPower, .FALSE. )
    CALL BCImpedance( x2, 'BOT', HSBot, fBot2, gBot2, iPower, .FALSE. )

    DetaDx = 0.0D0
    IF ( gBot1 /= 0.0D0 ) DetaDx = REAL( ( fBot2 / gBot2 - fBot1 / gBot1 ) ) / ( x2 - x1 )

    ! Scale the mode
    RN          = SqNorm - DrhoDx * Phi( 1 ) ** 2 + DetaDx * Phi( NTotal1 ) ** 2

    IF ( RN <= 0.0D0 ) THEN
       RN = -RN
       WRITE( PRTFile, * ) 'mode = ', mode
       WRITE( PRTFile, * ) 'Warning in KRAKEN - Normalize : Normalization constant non-positive; suggests grid too coarse'
    END IF

    ScaleFactor = 1.0D0 / SQRT( RN )
    IF ( Phi( iTurningPoint ) < 0.0D0 ) ScaleFactor = -ScaleFactor

    Phi            = ScaleFactor * Phi
    Slow           = ScaleFactor ** 2 * Slow * omega / SQRT( x )
    Perturbation_k = ScaleFactor ** 2 * Perturbation_k
    VG( mode )     = 1 / Slow

    CALL ScatterLoss( Perturbation_k, Phi, x )   ! Compute interfacial scatter loss

  END SUBROUTINE Normalize

  !**********************************************************************

  SUBROUTINE ScatterLoss( Perturbation_k, Phi, x )

    ! Figure scatter loss

    USE Scattering

    REAL     (KIND=8), INTENT( IN    ) :: x              ! eigenvalue
    REAL     (KIND=8), INTENT( IN    ) :: Phi( * )       ! eigenvector
    COMPLEX  (KIND=8), INTENT( INOUT ) :: Perturbation_k ! perturbation to wavenumber due to loss
    INTEGER           :: L, Medium
    REAL     (KIND=8) :: rho1, rho2, rhoInside, h2
    COMPLEX  (KIND=8) :: kx, eta1Sq, eta2Sq, U, PhiC

    kx    = SQRT( x )

    j = 1
    L = Loc( FirstAcoustic )

    Media: DO Medium = FirstAcoustic - 1, LastAcoustic   ! Loop over media

       ! Calculate rho1, eta1Sq, Phi, U

       IF ( Medium == FirstAcoustic - 1 ) THEN   ! Top properties
          SELECT CASE ( HSTop%BC )
          CASE ( 'A' )          ! Acousto-elastic
             rho1      = HSTop%rho
             eta1Sq    = x - omega2 / HSTop%cP ** 2
             U         = SQRT( eta1Sq ) * Phi( 1 ) / HSTop%rho
          CASE ( 'V' )          ! Vacuum
             rho1      = 1.0D-9
             eta1Sq    = 1.0D0
             rhoInside = rho( Loc( FirstAcoustic ) + 1 )
             U         = Phi( 2 ) / h( FirstAcoustic ) / rhoInside
          CASE ( 'R' )          ! Rigid
             rho1      = 1.0D+9
             eta1Sq    = 1.0D0
             U         = 0.0D0
          CASE DEFAULT          ! Tabulated
             rho1      = 0.0D0
             eta1Sq    = 0.0D0
             U         = 0.0D0
          END SELECT
       ELSE
          h2 = h( Medium ) ** 2
          j  = j + N( Medium )
          L  = Loc( Medium ) + N( Medium ) + 1

          rho1   = rho( L )
          eta1Sq = ( 2.0D0 + B1( L ) ) / h2 - x
          U      = ( -Phi( j - 1 ) - 0.5D0 * ( B1( L ) - h2 * x ) * Phi( j ) ) / ( h( Medium ) * rho1 )
       END IF

       ! Calculate rho2, eta2

       IF ( Medium == LastAcoustic ) THEN   ! Bottom properties
          SELECT CASE ( HSBot%BC )
          CASE ( 'A' )          ! Acousto-elastic
             rho2   = HSBot%rho
             eta2Sq = omega2 / HSBot%cP ** 2 - x
          CASE ( 'V' )          ! Vacuum
             rho2   = 1.0D-9
             eta2Sq = 1.0D0
          CASE ( 'R' )          ! Rigid
             rho2   = 1.0D+9
             eta2Sq = 1.0D0
          CASE DEFAULT          ! Tabulated
             rho2   = 0.0D0
             eta2Sq = 0.0D0
          END SELECT
       ELSE
          rho2   = rho( L + 1 )
          eta2Sq = ( 2.0D0 + B1( L + 1 ) ) / h( Medium + 1 ) ** 2 - x
       END IF

       PhiC = Phi( j )   ! convert to complex*16
       Perturbation_k = Perturbation_k + KupIng( SSP%sigma( Medium + 1 ), eta1Sq, rho1, eta2Sq, rho2, PhiC, U )

    END DO Media

    k( mode ) = Perturbation_k

  END SUBROUTINE ScatterLoss

  !**********************************************************************

  SUBROUTINE Solve1

    ! Uses Sturm sequences to isolate the eigenvalues
    ! and Brent's method to refine them

    USE SourceReceiverPositions   ! just to get Nfreq
    USE RootFinderBrent

    INTEGER                          :: iPower, NTotal, NzTab, IAllocStat
    REAL      (KIND=8)               :: x, x1, x2, xMin, xMax, Eps, Delta
    REAL      (KIND=8), ALLOCATABLE  :: xL( : ), xR( : )
    CHARACTER (LEN=80)               :: ErrorMessage

    CountModes = .TRUE.

    ! Determine number of modes

    xMin = 1.00001D0 * omega2 / cHigh ** 2
    
    CALL FUNCT( xMin, Delta, iPower )
    M = modeCount
    WRITE( PRTFile, * ) '               --- Number of modes = ', M

    IF ( ALLOCATED( xL ) ) DEALLOCATE( xL, xR )
    ALLOCATE( xL( M + 1 ), xR( M + 1 ) )

    IF ( iSet == 1 ) THEN
       IF ( ALLOCATED ( EVMat ) ) DEALLOCATE( EVMat, Extrap, k, VG )
       ALLOCATE( EVMat( NSets, M ), Extrap( NSets, M ), k( M ), VG( M ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( 'KRAKEN - Solve1', 'Insufficient memory (too many modes).' )
    END IF

    xMax = omega2 / cLow ** 2
    CALL FUNCT( xMax, Delta, iPower )
    M = M - modeCount

    IF ( M == 0 ) THEN   ! Open a dummy MODFile for possible use by FIELD3D

       LRecordLength = 32
       NzTab = 0

       ! Open MODFile and write header
       IF ( ifreq == 1 .AND. iProf == 1 ) THEN
          OPEN ( FILE = TRIM( FileRoot) //'.mod', UNIT = MODFile, ACCESS = 'DIRECT', RECL = 4 * LRecordLength, FORM = 'UNFORMATTED')
       END IF

       WRITE( MODFile, REC = iRecProfile     ) LRecordLength, Title, Nfreq, 1, NzTab, NzTab

!!! is this still correct with the new modefile format ???
       WRITE( MODFile, REC = iRecProfile + 6 ) M
       CALL ERROUT( 'KRAKEN', 'No modes for given phase speed interval' )
       iRecProfile = iRecProfile + 6
    END IF

    NTotal = SUM( N( FirstAcoustic : LastAcoustic ) )

    IF ( M > NTotal / 5 ) THEN
       WRITE( PRTFile, * ) 'Approximate number of modes =', M
       WRITE( PRTFile, * ) 'Warning in KRAKEN - Solve1 : Mesh too coarse to sample the modes adequately'
    END IF

    CALL Bisection( xMin, xMax, xL, xR )   ! Initialize upper and lower bounds

    ! Call ZBRENT to refine each eigenvalue in turn
    CountModes = .FALSE.

    DO mode = 1, M
       x1  = xL( mode )
       x2  = xR( mode )
       EPS = ABS( x2 ) * 10.0D0 ** ( 2.0 - PRECISION( x2 ) )
       CALL ZBRENTX( x, x1, x2, EPS, ErrorMessage, Funct )

       IF ( ErrorMessage /= ' ' ) THEN
          WRITE( PRTFile, * ) 'iSet, mode = ', iSet, mode
          WRITE( PRTFile, * ) 'Warning in KRAKEN-ZBRENTX : ', ErrorMessage
       END IF

       EVMat( iSet, mode ) = x
    END DO

    DEALLOCATE( xL, xR )

  END SUBROUTINE Solve1

  !**********************************************************************

  SUBROUTINE Bisection( xMin, xMax, xL, xR )

    ! Returns an isolating interval (xL, xR) for each eigenvalue
    ! in the given interval [ xMin, xMax ]

    INTEGER, PARAMETER           :: MaxBisections = 50
    INTEGER                      :: j, NZeros, NZer1, iPower
    REAL (KIND=8), INTENT( IN  ) :: xMin, xMax
    REAL (KIND=8), INTENT( OUT ) :: xL( M + 1 ), xR( M + 1 )
    REAL (KIND=8)                :: x, x1, x2, Delta

    CountModes = .TRUE.

    xL = xMin   ! initial left  boundary
    xR = xMax   ! initial right boundary

    CALL FUNCT( xMax, Delta, iPower )
    NZer1 = modeCount

    IF ( M == 1 ) RETURN   ! quick exit if only one mode is sought

    ModeLoop: DO mode = 1, M - 1  ! loop over eigenvalue

       ! Obtain initial guesses for x1 and x2
       IF ( xL( mode ) == xMin ) THEN
          x2 = xR( mode )
          x1 = MAX( MAXVAL( xL( mode + 1 : M ) ), xMin )

          ! Begin bisection (allowing no more than MaxBisections per mode)
          Bisect: DO j = 1, MaxBisections
             x = x1 + ( x2 - x1 ) / 2
             CALL FUNCT( x, Delta, iPower )
             NZeros = modeCount - NZer1

             IF ( NZeros < mode ) THEN   ! not too many zeros, this is a new right bdry
                x2         = x
                xR( mode ) = x
             ELSE                        ! this is a new left bdry
                x1 = x
                IF ( xR( NZeros + 1 ) >= x ) xR( NZeros + 1 ) = x
                IF ( xL( NZeros     ) <= x ) xL( NZeros     ) = x
             END IF

             ! when we have replaced the default, initial value for xL, we are done
             IF ( xL( mode ) /= xMin ) CYCLE ModeLoop
          END DO Bisect
       END IF
    END DO ModeLoop

  END SUBROUTINE Bisection

END PROGRAM KRAKEN
