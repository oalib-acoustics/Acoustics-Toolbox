PROGRAM KRAKEL
  
  ! Program for solving for ocean acoustic normal modes
  
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

  ! Sept. 24, 2020, MBP
  ! Someone asked about compiling this with the latest version of the Acoustics Toolbox.
  ! I made a huge number of changes to bring it (mostly) up to date
  ! It now compiles and runs and I think it's 99% of where it needs to be.
  ! However, further updates are needed to write the mode files in a form consistent with the rest of the package.
  ! In particular, KRAKENC writes the modes at depths representing the union of source and receiver depths.
  ! KRAKEL uses the finite difference grid and is writing to the wrong record number
  ! I'm stoppping here for now. A user interested in elastic normal modes is likely better off looking at the
  ! Matlab version of KRAKEL. The Matlab version has access to more recent solvers for Algebraic Eigenvalue Problems.
  ! However, it also needs some work ...

  USE krakelmod
  USE sspMod
  USE ReadEnvironmentMod
  USE SourceReceiverPositions
  USE RefCoef
  USE AttenMod
  IMPLICIT NONE
  INTEGER            :: IFirst, ILast, IREC
  REAL      (KIND=8) :: Error, freq0
  CHARACTER (LEN=80) :: FileRoot

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  NV( 1 : 5 ) = [ 1, 2, 4, 8, 16 ]
  CALL GetPar

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

     M = MaxM

     CALL UpdateSSPLoss( freq0, freq0 )
     CALL UpdateHSLoss(  freq0, freq0 )

     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) 'Mesh multiplier   CPU seconds'

     ! Main loop: solve the problem for a sequence of meshes

     DO iSet = 1, NSets

        N( 1 : SSP%NMedia ) = NG( 1 : SSP%NMedia ) * NV( iSet )
        h( 1 : SSP%NMedia ) = ( SSP%Depth( 2 : SSP%NMedia + 1 ) - SSP%Depth( 1 : SSP%NMedia ) ) / N( 1 : SSP%NMedia )
        hV( iSet ) = h( 1 )

        CALL SOLVE

        IF ( Error * 1000.0 * RMax < 1.0 ) THEN   ! Convergence achieved
           EXIT
        ELSE
           IF ( ISet == NSets ) WRITE( PRTFile, * ) 'Warning in KRAKEL : Too many meshes needed: check convergence'
        END IF

     END DO

  ! Solution complete: discard modes with phase velocity above cHigh

  DO Mode = M, 1, -1
     IF ( omega / DBLE( SQRT( Extrap( 1, Mode ) ) ) < cHigh ) THEN
        M = Mode
        EXIT
     END IF
  END DO

  ! Write eigenvalues to PRTFile and MODFile
  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) '   I          k             alpha          Phase Speed'

  k( 1 : M ) = SQRT( Extrap( 1, 1 : M ) )
  DO Mode = 1, M
     WRITE( PRTFile, "( I5, 3G18.10 )" ) Mode, k( Mode ), omega / DBLE( k( Mode ) )
  END DO

  WRITE( MODFile, REC = IRecProfile ) M

  IFirst = 1
  DO IREC = 1, 1 + ( 2 * M - 1 ) / LRecordLength
     ILast  = MIN( M, IFirst + LRecordLength / 2 - 1 )
     WRITE( MODFile, REC = IRecProfile + 1 + M + IREC ) CMPLX( k( IFirst : ILast ) )
     IFirst = ILast + 1
  END DO

  ! set record pointer to beginning of next mode set
  IRecProfile = IRecProfile + 3 + M + ( 2 * M - 1 ) / LRecordLength
END DO FreqLoop

CONTAINS

  !**********************************************************************

  SUBROUTINE GetPar

    ! Read in the ENVFile data

    REAL               :: zMin, zMax

    TITLE = 'KRAKEL-  '
    CALL ReadEnvironment( FileRoot, Title, freq0, MaxMedium, TopOpt, NG, BotOpt, cLow, cHigh, RMax, ENVFile, PRTFile )

    IF ( cLow <= 0.0 .OR. cHigh <= 0.0 .OR. cLow >= cHigh ) &
         CALL ERROUT( PRTFile, 'F', 'KRAKEL - GetPar', 'Need phase speeds cLow, cHigh > 0 and cLow < cHigh'  )

    zMin = SNGL( SSP%Depth( 1              ) )
    zMax = SNGL( SSP%Depth( SSP%NMedia + 1 ) )
    CALL ReadSzRz(    zMin, zMax )   ! Read source/receiver depths
    CALL ReadfreqVec( freq0, TopOpt( 6 : 6 ) )

    CLOSE( ENVFile )

    CALL ReadReflectionCoefficient( FileRoot, HSBot%BC, HSTop%BC, PRTFile )   ! Optionally read in Bot, Top reflection coefficients
    M = MaxM

  END SUBROUTINE GetPar

  !**********************************************************************

  SUBROUTINE INIT

    ! Initializes arrays defining difference equations
    ! Identical to KRAKENC apart for B1, B2, B3, B4 and rho

    LOGICAL           :: ElasticFlag = .FALSE.
    INTEGER           :: IAllocStat, ii, j, Medium, N1, NPoints
    REAL     (KIND=8) :: cMinV
    COMPLEX (KIND=8), ALLOCATABLE :: cP( : ), cS( : )
    COMPLEX (KIND=8)  :: lambda, mu
    CHARACTER (LEN=8) :: Task

    cMin     = HUGE( cMin )
    FirstAcoustic    = 0
    Loc( 1 ) = 0

    ! Allocate storage for finite-difference coefficients

    NPoints = SUM( N( 1 : SSP%NMedia ) ) + SSP%NMedia

    IF ( ALLOCATED( B1 ) ) DEALLOCATE( B1, B1C, B2, B3, B4, rho )
    ALLOCATE ( B1( NPoints ), B1C( NPoints ), B2( NPoints ), B3( NPoints ), B4( NPoints ), rho( NPoints ), &
         cP( NPoints ), cS( NPoints ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) &
         CALL ERROUT( PRTFile, 'F', 'KRAKEN - INIT', 'Insufficient memory: Reduce mesh.' )

    j = 0

    Media: DO Medium = 1, SSP%NMedia   ! Loop over media
       IF ( Medium /= 1 ) Loc( Medium ) = Loc( Medium - 1 ) + N( Medium - 1 ) + 1
       N1 = N( Medium ) + 1
       ii = Loc( Medium ) + 1

       IF ( Loc( Medium ) + N1 > MaxN ) THEN
          WRITE( PRTFile, * ) 'Insufficient storage for next mesh'
          WRITE( PRTFile, * ) 'PROGRAM ABORTING'
          STOP 'ERROR IN KRAKENC: Insufficient storage for next mesh'
       END IF

       ! EvaluateSSP reads in the data for a medium
       Task = 'TAB'
       CALL EvaluateSSP( cP( ii ), cS( ii ), rho( ii ), Medium, N1, freq, Task, ENVFile, PRTFile )

       IF ( cS( ii ) == ( 0.0, 0.0 ) ) THEN   ! Case of an acoustic medium

          SSP%Material( Medium ) = 'ACOUSTIC'
          IF ( FirstAcoustic == 0 ) FirstAcoustic = Medium
          LastAcoustic = Medium

          cMinV = MINVAL( DBLE( cp( ii:ii + N( Medium ) ) ) )
          cMin  = MIN( cMin, cMinV )
          B1( ii : ii + N( Medium ) ) = -2.0 +  h( Medium ) ** 2 * omega2 / cp( ii : ii + N( Medium ) ) ** 2

       ELSE                                   ! Case of an elastic medium

          IF ( SSP%sigma( Medium ) /= 0.0 ) THEN
             WRITE( PRTFile, * ) 'Rough elastic interface not allowed'
             WRITE( PRTFile, * ) 'PROGRAM ABORTING'
             STOP 'ERROR IN KRAKENC: Rough elastic interface not allowed'
          END IF

          SSP%Material( Medium ) = 'ELASTIC'
          ElasticFlag        = .TRUE.

          DO j = ii, ii + N( Medium )
             cMin = MIN( DBLE( cs( j ) ), cMin )

             lambda = rho( j ) * ( cP( j ) ** 2 - 2.0 * cS( j ) ** 2 )
             mu     = rho( j ) * cS( j ) ** 2

             B1( j  ) = 1.0 / mu
             B2( j  ) = 1.0 / ( lambda + 2.0 * mu )
             B3( j  ) = 4.0 * mu * ( lambda + mu ) / ( lambda + 2.0 * mu )
             B4( j  ) = lambda / ( lambda + 2.0 * mu )
             rho( j ) = omega2 * rho( j )

          END DO

       END IF
    END DO Media

    ! Bottom properties

    SELECT CASE ( HSBot%BC )
    CASE ( 'A' )
       IF ( HSBot%cS /= ( 0.0, 0.0 ) ) THEN    ! Elastic  bottom half-space
          SSP%Material( SSP%NMedia + 1 ) = 'ELASTIC'
          ElasticFlag = .TRUE.
          cMin = MIN( cMin, DBLE( HSBot%cS ) )
       ELSE                                    ! Acoustic bottom half-space
          SSP%Material( SSP%NMedia + 1 ) = 'ACOUSTIC'
          cMin = MIN( cMin, DBLE( HSBot%cP ) )
       END IF
    END SELECT

    ! Top properties

    SELECT CASE ( HSTop%BC )
    CASE ( 'A' )
       IF ( HSTop%cS /= ( 0.0, 0.0 ) ) THEN    ! Elastic  top half-space
          ElasticFlag = .TRUE.
          cMin = MIN( cMin, DBLE( HSTop%cS ) )
       ELSE                                    ! Acoustic top half-space
          cMin = MIN( cMin, DBLE( HSTop%cP ) )
       END IF
    END SELECT

    ! If elastic medium then reduce cMin for Scholte wave
    IF ( ElasticFlag ) cMin = 0.85 * cMin
    cLow = MAX( cLow, 0.99 * cMin )

  END SUBROUTINE INIT

  !**********************************************************************

  SUBROUTINE SOLVE

    ! Solves the eigenvalue problem at the current mesh
    ! and produces a new extrapolation

    INTEGER       :: j, Key
    REAL      (KIND=8)                :: T1, T2, F1, F2, TStart, TEnd, x1, x2

    CALL CPU_TIME( Tstart )
    CALL INIT     ! set up the finite difference mesh                      
    CALL SOLVE2   ! solve for the eigenvalues

    Extrap( iSet, 1 : M ) = EVMat( iSet, 1 : M ) 
    IF ( iSet == 1 ) CALL VECTOR   ! If this is the first mesh, compute the eigenvectors

    ! Richardson extrapolation to improve the accuracy          

    Error = 1.0D10          ! initialize error to a large number
    KEY   = 2 * M / 3 + 1   ! index of element used to check convergence

    IF ( iSet > 1 ) THEN
       T1 = Extrap( 1, KEY )

       DO j = iSet, 2, -1
          ModeLoop: DO Mode = 1, M
             F2 = Extrap( j,     Mode )
             F1 = Extrap( j - 1, Mode )
             x2 = NV( iSet  ) ** 2
             x1 = NV( j - 1 ) ** 2
             Extrap( j - 1, Mode ) = F2 - ( F1 - F2 ) / ( x2 / x1 - 1.0 )
          END DO ModeLoop
       END DO

       T2    = Extrap( 1, KEY )
       Error = ABS( T2 - T1 )
    END IF

    CALL CPU_TIME( Tend)   ! check elapsed time
    ET( iSet ) = Tend - Tstart
    WRITE( PRTFile, '( 1X, I8, 6X, G15.3 )' ) NV( iSet ), ET( iSet )

  END SUBROUTINE SOLVE

  !**********************************************************************

  SUBROUTINE SOLVE2

    ! Provides initial guess to root finder for each EVMat(I)

    USE RootFinderSecantMod

    INTEGER            :: MaxIteration, NTotal, ii, j, IAllocStat, Iteration
    REAL      (KIND=8) :: Tolerance
    REAL      (KIND=8) :: x, P( 10 )
    CHARACTER (LEN=80) :: ErrorMessage

    x     = omega2 / cLow ** 2
    MaxIteration = 500

    ! solve1 has already allocated space for the following unless the problem has shear
    IF ( .NOT. ALLOCATED( k ) ) THEN
       M = 3000   ! this sets the upper limit on how many modes can be calculated
       ALLOCATE( EVMat( NSets, M ), Extrap( NSets, M ), k( M ), VG( M ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
            CALL ERROUT( PRTFile, 'F', 'KRAKEN - SOLVE2', 'Insufficient memory (too many modes).' )
    END IF

    NTotal = SUM( N( FirstAcoustic : LastAcoustic ) )

    DO Mode = 1, M

       ! For first or second meshes, use a high guess
       ! Otherwise use extrapolation to produce an initial guess

       x = 1.00001 * x
       IF ( TopOpt( 4 : 4 ) == 'R' ) x = k( Mode ) ** 2
       IF ( iSet >= 2 ) THEN

          DO ii = 1, iSet-1
             P( ii ) = EVMat( ii, Mode)
          END DO

          IF ( iSet >= 3 ) THEN
             DO ii = 1, iSet-2
                DO j = 1, iSet - ii - 1
                   P( j ) = ( ( hV( iSet ) ** 2 - hV( j + ii ) ** 2 ) * P( j ) - ( hV( iSet ) ** 2 &
                        -hV( j ) ** 2 ) * P( j+1 ) ) / ( hV( j ) ** 2 - hV( j + ii ) ** 2 )
                END DO
             END DO
             x = P( 1 )
          END IF

       END IF

       ! Use the secant method to refine the eigenvalue
       Tolerance = ABS( x ) * 10.0 ** ( 3.0 - PRECISION( x ) )
       CALL RootFinderSecant( x, Tolerance, Iteration, MaxIteration, ErrorMessage, Funct )

       IF ( ErrorMessage /= ' ' ) THEN
          WRITE( PRTFile, * ) ErrorMessage
          WRITE( PRTFile, * ) 'for iSet, Mode = ',iSet, Mode
       END IF

       EVMat( iSet, Mode ) = x

       ! Toss out modes outside user specified spectral limits
       IF ( SQRT( omega2 ) / cHigh > REAL( SQRT( x ) ) ) THEN
          M = Mode - 1
          RETURN
       END IF

       IF ( Mode > NTotal / 5 ) THEN
          WRITE( PRTFile, * ) 'Caution: mesh is too coarse to compute all modes in spectrum'
          M = Mode
          RETURN
       END IF
    END DO   ! next mode

    ! Check whether we have hit the ceiling on max # modes

    IF ( M == MaxM ) THEN
       WRITE( PRTFile, * ) 'Too many modes in phase speed interval--'
       WRITE( PRTFile, * ) 'program will compute as many as it can'
    END IF

  END SUBROUTINE SOLVE2

  !**********************************************************************

  SUBROUTINE FUNCT( x, Delta, IPOW )

    ! FUNCT( x ) = 0 is the dispersion relation

    INTEGER,       INTENT( OUT ) :: IPOW
    INTEGER                      :: IPVT( MaxN ), NMat, lda, Ml, Mu, j, Info
    REAL (KIND=8)                :: A( 13, MaxN ), VDET( 2 )
    REAL (KIND=8), INTENT( IN  ) :: x
    REAL (KIND=8), INTENT( OUT ) :: Delta

    ! If we get to a k-value below the user's spectrum limit
    ! then force a zero by returning a new function value
    ! that is ten orders of magnitude less than the last value

    ! IF ( REAL( x ) <= omega2 / cHigh ** 2 ) THEN
    !    IPOW = IPOW - 10
    !    RETURN
    ! END IF

    CALL SETA( A, NMat, x ) ! Set up the matrix

    lda = 13
    ML  = 4
    MU  = 4

    CALL DGBTRF( NMat, NMat, ML, MU, A, lda, IPVT, INFO )  ! Factor the matrix

    IF ( INFO /= 0 ) THEN
       WRITE( PRTFile, * ) 'INFO = ', INFO
       STOP
    END IF

    CALL DGBDI( A, lda, NMat, ML, MU, IPVT, VDET )  ! Calculate the determinant

    Delta = VDET( 1 )
    IPOW  = INT( VDET( 2 ) )

    ! Deflate previous roots

    IF ( Mode > 1 ) THEN
       Delta = Delta / ( X - EVMat( iSet, Mode - 1 ) )
       IF ( Mode > 2 ) THEN
          ModeLoop: DO j = 1, Mode-2
             Delta = ( EVMat( iSet, Mode - 1 ) - EVMat( iSet, j ) ) * Delta / &
                  ( X - EVMat( iSet, j ) )
          END DO ModeLoop
       END IF
    END IF

  END SUBROUTINE FUNCT

  !**********************************************************************

  SUBROUTINE VECTOR

    ! Do inverse iteration to compute each of the eigenvectors
    ! and write these to the output file

    USE SourceReceiverPositions

    IMPLICIT NONE
    INTEGER          :: IPVT( MaxN ), ii, Info, j, Job, lda, ldb, Medium, ML, MU, NMat, NRHS, NTotal, NzTab
    REAL             :: z( MaxN )
    REAL    (KIND=8) :: A( 13, MaxN ), EVector( MaxN ), x
    COMPLEX          :: EVectorS( MaxN )
    COMPLEX (KIND=8) :: cx

    ! Tabulate Z-coordinates & compute NMat, NTotal
    NMat   = 0   ! size of the matrix
    NTotal = 0   ! number of z-points
    j      = 0

    DO Medium = 1, SSP%NMedia
       j      = j + 1
       z( j ) = SNGL( SSP%Depth( Medium ) )

       DO ii = 1, N( Medium )
          j      = j + 1
          z( j ) = z( j - 1 ) + SNGL( h( Medium ) )
       END DO

       IF ( SSP%Material( Medium ) == 'ACOUSTIC' ) THEN
          NMat = NMat +       N( Medium ) + 1
       ELSE
          NMat = NMat + 4 * ( N( Medium ) + 1 )
       END IF

       NTotal = NTotal + N( Medium ) + 1
    END DO

    ! Open MODFile and write header

    NzTab = NTotal   ! NzTab is the variable used in the other models giving the number of tabulated points
    
    IF ( ifreq == 1 ) THEN
       ! LRecordLength must not increase between profiles !!!
       LRecordLength = MAX( 2 * Nfreq, 2 * NzTab, 32, 3 * ( LastAcoustic - FirstAcoustic + 1 ) )   ! Logical record length in `longwords' (4 bytes)
       OPEN ( FILE = TRIM( FileRoot) //'.mod', UNIT = MODFile, ACCESS = 'DIRECT', STATUS = 'REPLACE', &
              RECL = 4 * LRecordLength, FORM = 'UNFORMATTED' )
    END IF

    IF ( ifreq == 1 ) THEN
       WRITE( MODFile, REC = IRecProfile     ) LRecordLength, Title, Nfreq, LastAcoustic - FirstAcoustic + 1, NzTab, NzTab
       WRITE( MODFile, REC = IRecProfile + 1 ) ( N( Medium ), SSP%Material( Medium ), Medium = FirstAcoustic, LastAcoustic )
       WRITE( MODFile, REC = IRecProfile + 2 ) ( REAL( SSP%Depth( Medium ) ), REAL( rho( Loc( Medium ) + 1 ) ), &
            Medium = FirstAcoustic, LastAcoustic )
       WRITE( MODFile, REC = IRecProfile + 3 ) freqVec( 1 : Nfreq )
       WRITE( MODFile, REC = IRecProfile + 4 ) z( 1 : NzTab )
       IRecProfile = IRecProfile + 5
    END IF

    ! top and bottom halfspace info (changes with frequency and possibly profile)
    WRITE( MODFile, REC = IRecProfile + 1 ) &
         HSTop%BC, CMPLX( HSTop%cP ), CMPLX( HSTop%cS ), REAL( HSTop%rho ), REAL( SSP%Depth( 1              ) ), &
         HSBot%BC, CMPLX( HSBot%cP ), CMPLX( HSBot%cS ), REAL( HSBot%rho ), REAL( SSP%Depth( SSP%NMedia + 1 ) )

    ! MAIN LOOP: For each eigenvalue call SINVIT to get eigenvector

    ModeLoop: DO Mode = 1, M
       x = 1.00000001 * EVMat( 1, Mode )

       CALL SETA( A, NMat, x ) ! Set up the matrix

       lda = 13
       ML  = 4
       MU  = 4
       CALL DGBTRF( NMat, NMat, ML, MU, A, lda, IPVT, INFO )  ! Factor the matrix

       IF ( INFO /= 0 ) THEN
          WRITE( PRTFile, * ) 'INFO = ', INFO
          STOP
       END IF

       JOB                 = 0
       EVector( 1 : NMat ) = 1.0
       NRHS                = 1      ! number of right hand sides
       ldb                 = NMat   ! first dimension of b
       CALL DGBTRS( 'N', NMat, ML, MU, NRHS, A, lda, IPVT, EVector, ldb, JOB )  ! Solve

       cx = x   ! convert to complex
       CALL Normalize( EVector, EVectorS, NMat, cx ) ! Normalize

       WRITE( MODFile, REC = 6 + Mode ) ( EVectorS( j ), j = 1, NMat )
    END DO ModeLoop

  END SUBROUTINE VECTOR

  !**********************************************************************

  SUBROUTINE SETA( A, NMat, x )

    ! ROUTINE FOR SETTING UP THE MATRIX A

    INTEGER,          INTENT( OUT ) :: NMat
    INTEGER                         :: ii, j, l, Medium, IPow
    REAL    (KIND=8), INTENT( OUT ) :: A( 13, MaxN )
    REAL    (KIND=8)                :: gammaP, gammaS, gammaP2, gammaS2, rmu, factor, h_rho, two_by_h
    REAL    (KIND=8), INTENT( IN  ) :: x
    COMPLEX (KIND=8)                :: F, G

    ! Compute size of matrix
    NMat = 0
    DO Medium = 1, SSP%NMedia
       IF ( SSP%Material( Medium ) == 'ACOUSTIC' ) THEN
          NMat = NMat +       N( Medium ) + 1
       ELSE
          NMat = NMat + 4 * ( N( Medium ) + 1 )
       END IF
    END DO

    A = 0.0 ! Zero out the matrix

    ! MAIN LOOP ( over each Medium )

    L = 0   ! check this (L was uninitialized in previous version)
    j = 0
    DO Medium = 1, SSP%NMedia
       IF ( Medium /= 1 ) THEN

          ! Acousto-elastic interface

          IF ( SSP%Material( Medium - 1 ) == 'ACOUSTIC' .AND. SSP%Material( Medium ) == 'ELASTIC' ) THEN
             ! R2 Condition
             L = L + 1
             j = j + 1

             h_rho = h( Medium - 1 ) * rho( Loc( Medium - 1 ) + 1 )

             A( 10, j - 1 ) = 1.0
             A(  9, j     ) = 0.5 * ( B1( L ) - h( Medium - 1 ) ** 2 * x )
             A(  7, j + 2 ) = omega2 * h_rho

             ! R3 Condition
             A( 10, j     ) = 0.0
             A(  7, j + 3 ) = 1.0
             ! R4 Condition
             A( 11, j     ) = 1.0
             A(  7, j + 4 ) = 1.0
          END IF

          ! Elasto-acoustic interface

          IF ( SSP%Material( Medium - 1 ) == 'ELASTIC' .AND. SSP%Material( Medium ) == 'ACOUSTIC' ) THEN
             L = Loc( Medium ) + 1
             j = j + 5
             ! R3 Condition
             A(  9, j - 2 ) = 1.0
             ! R4 Condition
             A(  9, j - 1 ) = 1.0
             A(  8, j     ) = 1.0
             ! R2 Condition
             h_rho = h( Medium ) * rho( Loc( Medium ) + 1 )
             A( 12, j - 3 ) = -h_rho * omega2
             A(  9, j     ) = 0.5 * ( B1( L ) - h( Medium ) ** 2 * x )
             A(  8, j + 1 ) = 1.0
          END IF

          ! Acoustic/acoustic interface

          IF ( SSP%Material( Medium - 1 ) == 'ACOUSTIC' .AND. SSP%Material( Medium ) == 'ACOUSTIC' ) THEN
             ! Continuity of normal displacement
             L = L + 1
             j = j + 1
             h_rho = h( Medium - 1 ) * rho( Loc( Medium - 1 ) + 1 )
             A( 10, j - 1 ) = 1.0 / h_rho
             A(  9, j     ) = 0.5 * ( B1( L ) - h( Medium - 1 ) ** 2 * x ) / h_rho

             L = L + 1
             h_rho = h( Medium ) * rho( Loc( Medium ) + 1 )
             A(  8, j + 1 ) = 0.5 * ( B1( L ) - h( Medium     ) ** 2 * x ) / h_rho
             A(  7, j + 2 ) = 1.0 / h_rho
             ! Continuity of pressure
             j = j + 1
             A( 10, j - 1 ) =  1.0
             A(  9, j     ) = -1.0
          END IF

          ! Elastic/elastic interface

          IF ( SSP%Material( Medium - 1 ) == 'ELASTIC' .AND. SSP%Material(Medium) == 'ELASTIC' ) THEN
             ! Continuity of everything
             j = j + 1
             A( 11, j     ) =  1.0
             A(  7, j + 4 ) = -1.0
             j = j + 1
             A( 11, j     ) =  1.0
             A(  7, j + 4 ) = -1.0
             j = j + 1
             A( 11, j     ) =  1.0
             A(  7, j + 4 ) = -1.0
             j = j + 1
             A( 11, j     ) =  1.0
             A(  7, j + 4 ) = -1.0
          END IF
       END IF

       IF ( SSP%Material( Medium ) == 'ACOUSTIC' ) THEN

          ! Acoustic section

          L = Loc( Medium ) + 1
          IF ( Medium == 1 ) j = 1
          DO ii = 1, N( Medium ) - 1
             j = j + 1
             L = L + 1
             A( 10, j - 1 ) = 1.0
             A(  9, j     ) = B1( L ) - h( Medium ) * h( Medium ) * x
             A(  8, j + 1 ) = 1.0
          END DO
       ELSE

          ! Elastic section

          L = Loc( Medium )
          two_by_h = 2.0 / h( Medium )
          DO ii = 1, N( Medium )
             L = L + 1
             j = j + 1
             A( 11, j ) = two_by_h
             A( 12, j ) = x * B4( L )
             A( 13, j ) = x * B3( L ) - rho( L )

             j = j + 1
             A( 10, j ) = -1.0
             A( 11, j ) = two_by_h
             A( 13, j ) = -rho( L )

             j = j + 1
             A(  9, j ) = B1( L )
             A( 11, j ) = two_by_h
             A( 12, j ) = x

             j = j + 1
             A(  9, j ) = B2( L )
             A( 10, j ) = -B4( L )
             A( 11, j ) = two_by_h

             L = L + 1
             j = j + 1
             A( 7, j ) = -two_by_h
             A( 8, j ) = x * B4( L )
             A( 9, j ) = x * B3( L ) - rho( L )

             j = j + 1
             A( 6, j ) = -1.0
             A( 7, j ) = -two_by_h
             A( 9, j ) = -rho( L )

             j = j + 1
             A( 5, j ) = B1( L )
             A( 7, j ) = -two_by_h
             A( 8, j ) = x

             j = j + 1
             A( 5, j ) =  B2( L )
             A( 6, j ) = -B4( L )
             A( 7, j ) = -two_by_h

             j = j - 4
             L = L - 1
          END DO
       END IF
    END DO   ! next medium

    ! Top BC

    SELECT CASE ( SSP%Material( 1 ) )
    CASE ( 'ELASTIC' )      ! Elastic medium
       SELECT CASE ( HSTop%BC )
       CASE ( 'R' )         ! Rigid top
          A( 9, 1 ) = 1.0
          A( 9, 2 ) = 1.0
       CASE ( 'V' )         ! Free (vacuum) top
          A( 7, 3 ) = 1.0
          A( 7, 4 ) = 1.0
       END SELECT
    CASE( 'ACOUSTIC' )      ! Acoustic medium
       CALL BCImpedance( x, 'TOP', HSTop, F, G, IPOW )   ! Top impedance
       IF ( G == 0.0 ) THEN ! Free (vacuum) top
          A( 9, 1 ) = 1.0
          A( 8, 2 ) = 0.0
       ELSE
          h_rho     = rho( 1 ) * h( 1 )
          A( 9, 1 ) = 0.5 * ( B1( 1 ) - h( 1 ) ** 2 * x ) / h_rho + F / G
          A( 8, 2 ) = 1.0 / h_rho
       END IF
    END SELECT

    ! Bottom BC

    SELECT CASE ( SSP%Material( SSP%NMedia ) )
    CASE ( 'ELASTIC' )     ! Elastic medium
       SELECT CASE ( HSBot%BC )
       CASE ( 'R' )        ! Rigid bottom
          A( 11, NMat - 3 ) = 1.0
          A( 11, NMat - 2 ) = 1.0
       CASE ( 'V' )        ! Free (vacuum) bottom
          A(  9, NMat - 1 ) = 1.0
          A(  9, NMat     ) = 1.0
       CASE ( 'A' )        ! Elastic bottom
          gammaP2 = x - omega2 / HSBot%cP ** 2
          gammaS2 = x - omega2 / HSBot%cS ** 2
          gammaP  = SQRT( gammaP2 )
          gammaS  = SQRT( gammaS2 )
          Rmu     = HSBot%rho * HSBot%cS ** 2
          factor  = -Rmu / ( gammaS * gammaP - x )

          A( 11, NMat - 3 ) = factor * gammaP * ( x - gammaS2 )
          A( 10, NMat - 2 ) = factor * ( 2.0 * gammaS * gammaP - ( gammaS2 + x ) )
          A( 12, NMat - 3 ) = x * A( 10, NMat - 2 )
          A( 11, NMat - 2 ) = factor * gammaS * ( x - gammaS2 )
          A(  9, NMat - 1 ) = 1.0
          A(  9, NMat     ) = 1.0
       END SELECT
    CASE( 'ACOUSTIC' )    ! Acoustic medium
       CALL BCImpedance( x, 'BOT', HSBot, F, G, IPOW )   ! Bottom impedance

       IF ( G == 0.0 ) THEN ! Free (vacuum) bottom
          A( 10, NMat - 1 ) = 0.0
          A(  9, NMat     ) = 1.0
       ELSE
          L = L + 1
          h_rho             = h( SSP%NMedia ) * rho( Loc( SSP%NMedia ) + 1 )
          A( 10, NMat - 1 ) = 1.0 / h_rho
          A(  9, NMat     ) = 0.5 * ( B1( L ) - h( SSP%NMedia ) ** 2 * x ) / h_rho - F / G
       END IF
    END SELECT

  END SUBROUTINE SETA

  !**********************************************************************

  SUBROUTINE BCImpedance( x, BotTop, HS, F, G, IPow )                                                   

    ! Compute Boundary Condition Impedance
    ! Same as KRAKENC except for ELASUP, ELASDN removed

    use RefCoef
    INTEGER, INTENT( OUT ) :: IPow
    INTEGER                :: IBot, ITop, iPower
    REAL    (KIND=8)       :: x, omega, rhoInside
    COMPLEX (KIND=8)       :: cx, kx, kz, gammaP, RCmplx, CInside
    COMPLEX (KIND=8), INTENT( OUT ) :: F, G
    CHARACTER      (LEN=3) :: BotTop
    TYPE( HSInfo )         :: HS
    TYPE(ReflectionCoef)   :: RInt

    cx   = x   ! complex version of x for routines that expect a complex argument
    IPow = 0 

    ! Get rho, c just Inside the boundary               
    ! There is at least one acoustic layer in the problem, except
    ! in the case where BOUNCE is used to get the refl. coef. for       
    ! a purely elastic stack of layers.                                 

    SELECT CASE ( BotTop )
    CASE ( 'TOP' ) 
       IF ( FirstAcoustic > 0 ) THEN
          Itop       = Loc( FirstAcoustic ) + N( LastAcoustic ) + 1
          rhoInside  = rho( Itop )
          cInside    = SQRT( omega2 * h( FirstAcoustic ) ** 2 / ( 2.0 + B1( Itop ) ) )
       END IF
    CASE ( 'BOT' )
       IF ( LastAcoustic > 0 ) THEN 
          Ibot       = Loc( LastAcoustic ) + N( LastAcoustic ) + 1 
          rhoInside  = rho( Ibot )
          cInside    = SQRT( omega2 * h( LastAcoustic ) ** 2 / ( 2.0 + B1( Ibot ) ) )
       END IF
    END SELECT

    SELECT CASE ( HS%BC )
    CASE ( 'V' )                   ! Vacuum
       F = 1.0
       ! G = -CI * PEKRT( omega2 / cInside ** 2 - x ) * sigma(1) ** 2
       G = 0.0
    CASE ( 'R' )                    ! Rigid
       F = 0.0
       G = 1.0
    CASE ( 'A' )                    ! Acousto-elastic half-space
       gammaP = SQRT( x - omega2 / HS%cP ** 2 )
       F = 1.0
       G = HS%rho / gammaP
    CASE ( 'F' )                    ! Tabulated reflection coefficient
       ! Compute the grazing angle THETA
       kx         = SQRT( x ) 
       kz         = SQRT( omega2 / CInside ** 2 - x ) 

       RInt%theta = RadDeg * DATAN2( DBLE( kz ), DBLE( kx ) )

       ! Evaluate R( ThetaInt )
       IF ( BotTop == 'TOP' ) THEN
          CALL InterpolateReflectionCoefficient( RInt, RTop, NTopPts, PRTFile )
       ELSE
          CALL InterpolateReflectionCoefficient( RInt, RBot, NBotPts, PRTFile )
       END IF

       ! Convert R(THETA) to ( f, g ) in Robin BC
       RCmplx = RInt%R * EXP( i * RInt%phi )
       F      = 1.0
       G      = ( 1.0 + RCmplx ) / ( i * kz * ( 1.0 - RCmplx ) )
    CASE ( 'P' )                    ! Precalculated reflection coef
       CALL InterpolateIRC( cx, F, G, IPower, xTab, FTab, GTab, ITab, NkTab )
    END SELECT

    IF ( BotTop == 'TOP' ) G = -G   ! A top BC has the sign flipped relative to a bottom BC

  END SUBROUTINE BCImpedance

  !**********************************************************************

  SUBROUTINE Normalize( EVector, EVectorS, NMat, x )

    ! Convert the EVector to ( u, w, tau_zx, tau_zz ) and normalize

    INTEGER          :: ii, Medium, NMat, j, kk
    REAL    (KIND=8) :: EVector( MaxN ), SqNorm2, rhoMedium
    COMPLEX          :: EVectorS( MaxN ), CIk
    COMPLEX (KIND=8) :: R1, R2, R3, R4, x, SqNorm, Contrib, RN, ScaleFac, gammaP

    CIk = i * SQRT( x )

    ! Scale down the mode

    SqNorm2 = 0.0

    DO ii = 1, NMat
       SqNorm2 = MAX( SqNorm2, ABS( EVector( ii ) ) )
    END DO

    EVector( 1 : NMat ) = EVector( 1 : NMat ) / SqNorm2

    ! Loop to compute norm of eigenvector

    SqNorm = 0.0
    j      = 1
    kk     = 1

    Media: DO Medium = 1, SSP%NMedia
       rhoMedium = rho( Loc( Medium ) + 1 )

       DO ii = 1, N( Medium ) + 1
          IF ( SSP%Material( Medium ) == 'ELASTIC' ) THEN
             R1 = EVector( kk     )
             R2 = EVector( kk + 1 )
             R3 = EVector( kk + 2 )
             R4 = EVector( kk + 3 )

             Contrib = h( Medium ) * ( -x * B3( j ) * R1 ** 2 + B4( j ) * R1 * R4 - R2 * R3 )

             EVectorS( kk     ) = CIk * EVector( kk     )
             EVectorS( kk + 1 ) =       EVector( kk + 1 )
             EVectorS( kk + 2 ) = CIk * EVector( kk + 2 )
             EVectorS( kk + 3 ) =       EVector( kk + 3 )

             kk = kk + 4
          ELSE
             Contrib        = -h( Medium ) * EVector( kk ) ** 2 / ( rhoMedium * omega2 )
             EVectorS( kk ) = -EVector( kk )   ! sign flip implies this is tau_zz not pressure
             kk             = kk + 1
          END IF

          IF ( ii == 1 .OR. ii == N( Medium ) + 1 ) Contrib = 0.5 * Contrib
          SqNorm = SqNorm + Contrib
          j      = j + 1
       END DO
    END DO Media

    ! Bottom half-space contribution

    IF ( HSBot%BC == 'A' ) THEN
       IF ( SSP%Material( SSP%NMedia + 1 ) == 'ELASTIC' ) THEN
          WRITE( PRTFile, * ) 'Elastic halfspace normalization not implemented'
       ELSE
          gammaP  = SQRT( x - omega2 / HSBot%cP ** 2 )
          Contrib = -EVector( NMat ) ** 2 / ( 2.0 * gammaP * HSBot%rho * omega2 )
       END IF
    END IF

    SqNorm = SqNorm + Contrib

    ! Scale the mode

    RN       = -omega2 * SqNorm
    ScaleFac = 1.0 / SQRT( RN )
    EVectorS( 1 : NMat ) = ScaleFac * EVectorS( 1 : NMat )

  END SUBROUTINE Normalize
END PROGRAM KRAKEL
