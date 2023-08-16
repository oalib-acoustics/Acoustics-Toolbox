PROGRAM FIELD

  ! Generate file with the complex pressure field (these are replica vectors in the terminology of matched-field processing)
  !
  ! Useage: Field FileRoot
  !
  ! where
  !    FileRoot.mod contains the mode  file (input)
  !    FileRoot.shd contains the shade file (output)

  USE MathConstants
  USE SourceReceiverPositions
  USE BeamPattern
  USE interpolation
  USE SubTabulate
  USE ReadModes
  USE EvaluateMod
  USE EvaluateCMMod
  USE EvaluateADMod
  USE RWSHDFile
  USE FatalError

  IMPLICIT NONE
  INTEGER,   PARAMETER :: FLPFile = 5, PRTFile = 6, SHDFile = 25, MaxM = 20000, MaxNfreq = 1000 ! MaxM also in Evaluate, EvaluateAD, EvaluateCM
  INTEGER              :: ifreq, iirz, iProf, NProf, NRro, IAllocStat, iRec, iS, M, MLimit, MSrc, IOStat
  REAL                 :: zMin, zMax
  REAL (KIND=8)        :: freq0, atten = 0, omega, c0
  COMPLEX              :: k( MaxM )
  CHARACTER   (LEN=50) :: Opt
  CHARACTER   (LEN=80) :: SHDTitle, Title, FileRoot
  CHARACTER   (LEN=1 ) :: Comp
  CHARACTER   (LEN=10) :: PlotType = '          '
  REAL (KIND=8), ALLOCATABLE :: rProf( : )   ! Profile ranges
  REAL (KIND=8), ALLOCATABLE :: Rro( : )     ! Receiver range-offsets (array tilt)
  COMPLEX,       ALLOCATABLE :: phiS( :, : ), phiR( :, : ), P( :, : ), C( : ), Ptemp( : )
  REAL (KIND=8), ALLOCATABLE :: kz2( : ), thetaT( : ), S( : )

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  ! Open the print file
  OPEN( UNIT = PRTFile, FILE = 'field.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )
  IF ( IOStat /= 0 ) THEN
     WRITE( PRTFile, * ) 'PRTFile = ', TRIM( FileRoot ) // '.prt'
     CALL ERROUT( 'FIELD', 'Unable to open PRTFile' )
  END IF

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) '__________________________________________________________________________'
  WRITE( PRTFile, * ) 'Running FIELD'
  WRITE( PRTFile, * ) 'Sums modes, producing pressure'
  WRITE( PRTFile, * )

  ! open the field paramaters file (FLPFile)
  OPEN( FILE = TRIM( FileRoot ) // '.flp', UNIT = FLPFile, STATUS = 'OLD', FORM = 'FORMATTED', IOSTAT = IOStat, ACTION = 'READ' )
  IF ( IOStat /= 0 ) THEN
     WRITE( PRTFile, * ) 'FLPFile = ', TRIM( FileRoot ) // '.flp'
     CALL ERROUT( 'FIELD', 'Unable to open FLPFile' )
  END IF

  SHDTitle( 1 : 1 ) = '$'

  READ( FLPFile, * ) SHDTitle
  READ( FLPFile, * ) Opt
  READ( FLPFile, * ) MLimit
  MLimit = MIN( MaxM, MLimit )

  SELECT CASE ( Opt( 1 : 1 ) )
  CASE ( 'X' )
     WRITE( PRTFile, * ) 'Line source'
  CASE ( 'R' )
     WRITE( PRTFile, * ) 'Point source'
  CASE ( 'S' )
     WRITE( PRTFile, * ) 'Scaled cylindrical (Point source with cylindrical spreading removed)'
  CASE DEFAULT
     CALL ERROUT( 'FIELD', 'Unknown value for Option( 1 : 1 ); should be X, R, or S' )
  END SELECT

  Comp = Opt( 3 : 3 )

  SELECT CASE ( Opt( 3 : 3 ) )
  CASE ( '*' )
     WRITE( PRTFile, * ) 'Using a source beam pattern'
  CASE ( 'O', ' ' )
     WRITE( PRTFile, * ) 'Omnidirectional source'
  CASE DEFAULT
     CALL ERROUT( 'FIELD', 'Unknown option for source beam pattern' )
  END SELECT

  SELECT CASE ( Opt( 4 : 4 ) )
  CASE ( 'C', ' ' )
     WRITE( PRTFile, * ) 'Coherent mode addition'
  CASE ( 'I' )
     WRITE( PRTFile, * ) 'Incoherent mode addition'
  CASE DEFAULT
     CALL ERROUT( 'FIELD', 'Unknown option for coherent vs. incoherent mode addition' )
  END SELECT

  SBPFlag = Opt( 3 : 3 )
  CALL ReadPAT( FileRoot, PRTFile )  ! Read Source Beam Pattern

  ! Read profile ranges (where a new mode set takes effect)

  CALL ReadVector( NProf, RProf, 'Profile ranges, RProf', 'km' )
  RProf = RProf / 1000.0   ! convert m back to km (undoing what ReadVector did)

  IF ( NProf == 1      ) THEN
     WRITE( PRTFile, * ) 'Range-independent calculation'
  ELSE
     WRITE( PRTFile, * ) 'Range-dependent calculation'
     SELECT CASE ( Opt( 2 : 2 ) )
     CASE ( 'C' )
        WRITE( PRTFile, * ) 'Coupled modes'
        IF ( Opt( 4 : 4 ) == 'I' ) &
           CALL ERROUT( 'FIELD', 'Coupled mode calculations do not support the incoherent addition of modes' )

     CASE ( 'A' )
        WRITE( PRTFile, * ) 'Adiabatic modes'
     CASE DEFAULT
        CALL ERROUT( 'FIELD', 'Unknown option for type of calculation; should be Adiabatic or Coupled modes' )
     END SELECT
  END IF

  ! EvaluateAD/EvaluateCM need a profile at zero range
  IF ( rProf( 1 ) /= 0.0 ) CALL ERROUT( 'FIELD', 'The first profile must be at a range of 0 km' )

  CALL ReadRcvrRanges           ! Read receiver ranges, Rr

  zMin = -HUGE( zMin )
  zMax = +HUGE( zMax )
  CALL ReadSzRz( zMin, zMax )   ! Read source/receiver depths, Sz, Rz

  CALL ReadVector( NRro, Rro, 'receiver range offsets (array tilt), Rro', 'm' )  ! Read receiver ranges (offsets from vertical)

  IF ( NRro /= Pos%Nrz ) THEN
     WRITE( PRTFile, * ) 'NRro, NRz = ', NRro, Pos%Nrz
     CALL ERROUT( 'FIELD', 'Number of range offsets, NRro, not equal to number of depths, NRz' )
  ENDIF

  WRITE( PRTFile, * )
  ALLOCATE( phiS( MaxM, Pos%NSz ), phiR( MaxM, Pos%NRz ), C( MaxM ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( 'FIELD', 'Dynamic memory allocation failed: Too many receiver depths' )

  ALLOCATE ( P( Pos%NRz, Pos%NRr ), Ptemp( Pos%NRr ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( 'FIELD', ' Insufficient memory to allocate P( NRz, NRr )' )

  !  *** Read in modes ***

  iProf = 1
  allocate( freqVec( MaxNfreq ) )   !!! must match dimension in GetModes.f90 contained in ReadModes.f90
  !!! the treatment of multiple frequencies is ugly here because we don't know the number of frequencies until we read the mode header
  Nfreq = MaxNfreq

  FreqLoop: DO ifreq = 1, MaxNfreq
     IF ( ifreq > Nfreq ) EXIT FreqLoop  !!! Nfreq is updated afer calling GetModes below
     CALL GetModes( FileRoot, iProf, ifreq, MaxM, Pos%Sz, Pos%NSz, 'N' , k, phiS, MSrc, freqVec, Nfreq, Title )
     CALL GetModes( FileRoot, iProf, ifreq, MaxM, Pos%Rz, Pos%NRz, Comp, k, phiR, MSrc, freqVec, Nfreq, Title )

     ! Generate header

     IF ( SHDTitle( 1 : 1 ) == '$' ) SHDTitle = Title   ! If a title was not given in the .flp file, take it from the mode file

     IF ( ifreq == 1 ) THEN
        CALL WriteHeader( TRIM( FileRoot ) // '.shd', SHDTitle, freq0, atten, PlotType )
        iRec = 10
     END IF

     ! *** For each source depth, evaluate and write the field ***

     SourceDepths: DO iS = 1, Pos%Nsz
        M = MIN( MLimit, MSrc )   ! Set number of propagating modes

        C( 1 : MSrc ) = phiS( 1 : MSrc, is )

        ! apply the source beam pattern
        IF ( SBPFlag == '*' .AND. iS == 1 ) THEN
           ALLOCATE( kz2( MSrc ), thetaT( MSrc ), S( MSrc ) )
           c0    = 1500   !!! reference sound speed, should be speed at the source depth
           omega = 2 * pi * freqVec( ifreq )
           kz2   = REAL( omega ** 2 / c0 ** 2 - k( 1 : MSrc ) ** 2 )      ! vertical wavenumber squared
           WHERE ( kz2 < 0 ) kz2 = 0                                      ! remove negative values

           thetaT = RadDeg * ATAN( SQRT( kz2 ) / REAL( k( 1 : MSrc ) ) )  ! calculate the angle in degrees
           CALL interp1( SrcBmPat( :, 1 ), SrcBmPat( :, 2 ), thetaT, S )
           C( 1 : Msrc ) = C( 1 : Msrc ) * REAL( S )         ! apply the shading
        END IF

        IF ( NProf == 1 ) THEN   ! Range-independent case
           CALL Evaluate(                                     C, phiR,         Pos%NRz, Pos%Rr, Pos%NRr, Rro, k, M, Opt, P )
        ELSE
           SELECT CASE ( Opt( 2 : 2 ) )
           CASE ( 'C' )   ! Coupled mode theory
              CALL EvaluateCM( FileRoot, rProf, NProf,        C, phiR, Pos%Rz, Pos%NRz, Pos%Rr, Pos%NRr,      k, M, Opt, P )
           CASE ( 'A' )   ! Adiabatic mode theory
              CALL EvaluateAD( FileRoot, rProf, NProf, ifreq, C, phiR, Pos%Rz, Pos%NRz, Pos%Rr, Pos%NRr,         M, Opt, P )
           END SELECT
        ENDIF

        ! write out the field
        RcvrDepths: DO iirz = 1, Pos%Nrz
           IRec  = IRec + 1
           Ptemp = P( iirz, 1 : Pos%Nrr )   ! temporary variable, avoids compiler warning
           WRITE( SHDFile, REC = IRec ) Ptemp
        END DO RcvrDepths

     END DO SourceDepths

  END DO FreqLoop

  ! do some clean up (that the compiler should really take care of automatically)
  CLOSE( SHDFile )
  IF ( ALLOCATED( kz2 ) ) DEALLOCATE( kz2, thetaT, S )

  WRITE( PRTfile, * ) 'Field completed successfully'

END PROGRAM FIELD
