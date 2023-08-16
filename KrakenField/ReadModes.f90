MODULE ReadModes

  USE PekRoot
  USE FatalError

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: PRTFile = 6
  INTEGER, PARAMETER, PRIVATE :: MODeFile = 30, MaxN = 100001, MaxMedium = 501, MaxNfreq = 1000
  REAL,    PARAMETER, PRIVATE :: pi = 3.141592
  INTEGER,            PRIVATE :: ir, iz
  COMPLEX,            PRIVATE :: kTop2, kBot2

CONTAINS
  
  SUBROUTINE GetModes( FileRoot, IProf, ifreq, MaxM, rd, Nrd, Comp, k, PhiR, M, freqVec, Nfreq, Title )                                   

  ! Read in modes and extract values at rcvr depths                   

  ! If IProf = 1 then the file is opened and IRecProfile is set to the first record
  ! For each new profile, the modefile starts with 5 descriptive header records
  ! Then, for each frequency:
  !   A record with M, the number of modes
  !   A record with the halfspace properties (attenuation may change with frequency)
  !   The actual modes with one per record
  !   The eigenvalues, k, folded to fit the record length

  USE calculateweights

  INTEGER, INTENT( IN  ) :: IProf, ifreq        ! profile and frequency indices for the mode set
  INTEGER, INTENT( IN  ) :: MaxM                ! row dimension of PhiR in calling program
  INTEGER, INTENT( IN  ) :: Nrd                 ! number of receiver depths (where modes are sampled)
  REAL,    INTENT( IN  ) :: rd( Nrd )           ! receiver depths
  CHARACTER (LEN=*), INTENT( IN  ) :: FileRoot  ! Name of the mode file
  INTEGER, INTENT( OUT ) :: M                   ! number of modes
  REAL  (KIND=8)         :: freqVec( MaxNfreq ) ! frequencies
  COMPLEX, INTENT( OUT ) :: k( MaxM )           ! modal wavenumbers
  COMPLEX, INTENT( OUT ) :: PhiR( MaxM, Nrd )   ! mode shapes
  CHARACTER (LEN=*), INTENT( OUT ) :: Title     ! title in mode file
  CHARACTER (LEN=1), INTENT( IN  ) :: Comp      ! component (vertical, horizontal,...) extracted (ignored if no shear)
  INTEGER                :: Nfreq, IRecProfile = 1, LRecl, Mode, NMat, NMedia, NTot, N( MaxMedium ), ird( Nrd ) 
  REAL                   :: Z( MaxN ), W( Nrd ), Depth( MaxMedium ), rho( MaxMedium ), Tolerance, wT
  REAL                   :: rhoT, rhoB, depthT, depthB   ! halfspace density, depth
  COMPLEX                :: cpT, csT, cpB, csB      ! halfspace pressure and shear speeds
  COMPLEX                :: PhiT( MaxN )
  CHARACTER (LEN=8)      :: Material( MaxMedium )
  CHARACTER (LEN=1)      :: BCTop, BCBot            ! boundary condition type for top and bottom
  CHARACTER (LEN=7)      :: Model                   ! KRAKEN or KRAKENC
  SAVE iRecProfile

  ! Read the mode header and do some checks
  ! Each new profile has its own mode header

  CALL ReadModeHeader( FileRoot, IProf, IRecProfile, LRecl, Title, freqVec, Nfreq, NMedia, NTot, NMat, N, Material, Depth, rho, Z )
  CALL Weight( Z, NTot, rd, Nrd, W, ird )   ! Locate indices of receiver points

  ! Read in eigenvalues, k( I ), skipping through frequencies until we get to the right one
  CALL ReadWavenumbers( IRecProfile, ifreq, k, M, MaxM, LRecL )

  IF ( M == 0 ) RETURN

  ! Abort if the storage is inadequate
  IF ( M > MaxM ) THEN 
     WRITE( PRTFile, * ) 'M = ', M, '   MaxM = ', MaxM
     CALL ERROUT( 'ReadModeHeader', 'Insufficient storage to read all the modes: increase MaxM' )
     M = MaxM
  ENDIF

  ! Read top and bottom halfspace info
  READ( ModeFile, REC = IRecProfile + 1 ) BCTop( 1 : 1 ), cPT, cST, rhoT, DepthT, BCBot( 1 : 1 ), cPB, cSB, rhoB, DepthB

  IF ( BCTop( 1 : 1 ) == 'A' ) kTop2 = ( 2.0 * pi * REAL( freqVec( ifreq ), KIND = 4 ) / cPT ) **2 
  IF ( BCBot( 1 : 1 ) == 'A' ) kBot2 = ( 2.0 * pi * REAL( freqVec( ifreq ), KIND = 4 ) / cPB ) **2
  
  ! Loop over receiver depths to check for safe interpolation
  ! Receivers must be within a fraction of a wavelength
  ! of tabulated pts. We accept one wavelength at 1500 m/s

  Tolerance = 1500.0 / REAL( freqVec( ifreq ), KIND = 4 )
  Model     = Title( 1 : 7 )   ! Should be KRAKEN or KRAKENC

  RcvrDepth: DO ir = 1, Nrd 

     iz = ird( ir ) 
     WT = ABS( MIN( W( ir ), 1.0 - W( ir ) ) ) 

     IF ( rd( ir ) < DepthT ) THEN        ! Receiver in upper halfspace
        IF ( cST /= 0.0 .OR. BCTop( 1 : 1 ) /= 'A' ) THEN 
           WRITE( PRTFile, * ) 'Receiver depth: ', rd( ir )
           WRITE( PRTFile, * ) 'Highest valid depth: ', DepthT
           WRITE( PRTFile, * ) 'Warning in GetMode : Receiver above depth of top'
        ENDIF

     ELSE IF ( rd( ir ) > DepthB ) THEN   ! Receiver in lower halfspace
        IF ( cSB /= 0.0 .OR. BCBot( 1 : 1 ) /= 'A' ) THEN 
           WRITE( PRTFile, * ) 'Receiver depth: ', rd( ir )
           WRITE( PRTFile, * ) 'Lowest valid depth: ', DepthB
           WRITE( PRTFile, * ) 'Warning in GetMode : Receiver below depth of bottom'
        ENDIF

     ELSE IF ( NTot > 1 ) THEN            ! Receiver between two grid points or large extrapolation
        IF ( WT * ( Z( iz + 1 ) - Z( iz ) ) > Tolerance ) THEN
           WRITE( PRTFile, * ) 'Receiver depth: ', rd( ir )
           WRITE( PRTFile, * ) 'Nearest depths: ', Z( iz ), Z( iz + 1 )
           WRITE( PRTFile, * ) 'Tolerance: ', Tolerance 
           WRITE( PRTFile, * ) 'Warning in GetMode : Modes not tabulated near requested pt.'
        ENDIF
     ELSE                                 ! Receiver near a single grid point
        IF ( ABS( rd( ir ) - Z( iz ) ) > Tolerance ) THEN 
           WRITE( PRTFile, * ) 'Rd, Tabulation depth ', rd( ir), Z( iz )
           WRITE( PRTFile, * ) 'Tolerance: ', Tolerance 
           WRITE( PRTFile, * ) 'Warning in GetMode : Modes not tabulated near requested pt.'
        ENDIF
     ENDIF

  END DO RcvrDepth

  ! Read in the modes
  ModeLoop: DO Mode = 1, M
     CALL ReadOneMode( Mode, IRecProfile, NTot, NMat, W, ird, N, Material, NMedia, Comp, &
          DepthT, BCTop, DepthB, BCBot, rd, Nrd, k, PhiT, Model )
     PhiR( Mode, 1 : Nrd ) = PhiT( 1 : Nrd )
  END DO ModeLoop

  ! advance the record number by M modes + the wavenumber records (pointing to the next set)
  IRecProfile = IRecProfile + 3 + M + ( 2 * M - 1 ) / LRecl

END SUBROUTINE GetModes

!**********************************************************************C

SUBROUTINE ReadModeHeader( FileRoot, iProf, IRecProfile, LRecl, Title, freqVec, Nfreq, &
     NMedia, NTot, NMat, N, Material, Depth, rho, Z )

  ! Reads the header information from ModeFile                          
  ! Note T suffix means top
  !      B suffix means bottom                                        

  ! iProf     is a profile number; it is not fully used
  ! FileRoot  is the user-provided file name                             
  ! FileNameT is the temporary name we build
  ! These have to be two separate variables to ensure there
  ! is space allocated to construct the file name even when           
  ! the user calls us with FileName = ' ' 

  ! IRecProfile must point to the first record of the profile                             

  INTEGER,            INTENT(IN ) :: iProf
  CHARACTER (LEN= *), INTENT(IN ) :: FileRoot
  INTEGER,            INTENT(OUT) :: N( * ) , NMedia, NTot, NMat
  INTEGER,            INTENT(OUT) :: LRecL, Nfreq
  INTEGER,            INTENT(INOUT) :: IRecProfile
  REAL,               INTENT(OUT) :: Depth( * ), Z( * ), rho( * )
  REAL      (KIND=8), INTENT(OUT) :: freqVec( * )
  CHARACTER (LEN= *), INTENT(OUT) :: Title
  CHARACTER (LEN= 8), INTENT(OUT) :: Material( * )
  LOGICAL              :: OpenFlag
  INTEGER              :: iostat, Medium
  CHARACTER (LEN=80)   :: FileNameT

  ! open ModeFile, if not already open
  FileNameT = TRIM( FileRoot ) // '.mod'
  
  INQUIRE( FILE = FileNameT, OPENED = OpenFlag )
  IF ( .NOT. OpenFlag ) THEN
     OPEN( UNIT = ModeFile, FILE = FileNameT, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', &
          RECL = 100, IOSTAT = iostat )
     IF ( IOSTAT /= 0 ) THEN
        WRITE( PRTFile, * ) 'Mode file = ', FileNameT
        CALL ERROUT( 'GetMode - ReadModeHeader', 'Unable to open the mode file' )
     END IF

     READ( ModeFile, REC = 1 ) LRecL
     CLOSE( UNIT = ModeFile )
     OPEN(  UNIT = ModeFile, FILE = FileNameT, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', &
          RECL = 4 * LRecL, IOSTAT = iostat )
  END IF

  ! If this is the first profile, reset the record counter to the beginning of the file
  IF ( iProf == 1 ) THEN
     IRecProfile = 1   ! set counter pointing to the first record to read
  END IF

  ! Read header info
  READ( ModeFile, REC = IRecProfile     ) LRecL, Title( 1 : 80 ), Nfreq, NMedia, NTot, NMat
  READ( ModeFile, REC = IRecProfile + 1 ) ( N( Medium ), Material( Medium ), Medium = 1, NMedia)
  READ( ModeFile, REC = IRecProfile + 2 ) ( Depth( Medium ), rho( Medium ),  Medium = 1, NMedia )
  READ( ModeFile, REC = IRecProfile + 3 ) freqVec( 1 : Nfreq )
  READ( ModeFile, REC = IRecProfile + 4 ) Z( 1 : NTot ) 
  IRecProfile = IRecProfile + 5

END SUBROUTINE ReadModeHeader

!**********************************************************************

SUBROUTINE ReadWavenumbers( IRecProfile, ifreq, k, M, MaxM, LRecL )

  ! Read in wavenumbers (eigenvalues), k( I )
  ! They come at the end of the ModeFile, because they're calculated after the eigenvectors

  INTEGER, INTENT( IN    ) :: MaxM, LRecl, ifreq
  INTEGER, INTENT( INOUT ) :: IRecProfile
  INTEGER, INTENT( OUT   ) :: M
  COMPLEX, INTENT( OUT   ) :: k( * )
  INTEGER                  :: IFirst, ILast, IRec, ifreqT

  ! Read the number of modes, M
  DO ifreqT = 1, ifreq
     IF ( ifreqT > 1 ) IRecProfile = IRecProfile + ( 3 + M + ( 2 * M - 1 ) / LRecl )   ! advance to appropriate record
     READ( ModeFile, REC = IRecProfile ) M
  END DO

  IF ( M > 0 ) THEN   ! Are there any modes?
     IFirst = 1 
     DO IRec = 1, 1 + ( 2 * MIN( M, MaxM ) - 1 ) / LRecL
        ILast = MIN( M, IFirst + LRecL / 2 - 1 )
        READ( ModeFile, REC = IRecProfile + 1 + M + IRec ) k( IFirst : ILast )
      IFirst = ILast + 1
     END DO
  END IF  
END SUBROUTINE ReadWavenumbers

!**********************************************************************

SUBROUTINE ReadOneMode( Mode, IRecProfile, NTot, NMat, W, ird, N, Material, NMedia, Comp, &
     &   DepthT, BCTop, DepthB, BCBot, rd, Nrd, k, PhiR , Model )

  ! Read in a single eigenfunction and extract receiver values
  ! Results are returned in PhiR

  INTEGER,           INTENT( IN ) :: Mode, NMedia, Nrd, N( NMedia ), ird( Nrd ), IRecProfile, NTot, NMat
  REAL,              INTENT( IN ) :: rd( Nrd ), W( Nrd ), DepthT, DepthB
  COMPLEX,           INTENT( IN ) :: k( * )
  CHARACTER (LEN=1), INTENT( IN ) :: Comp
  CHARACTER (LEN=8), INTENT( IN ) :: Material( NMedia )
  CHARACTER (LEN=1), INTENT( IN ) :: BCTop, BCBot
  CHARACTER (LEN=7), INTENT( IN ) :: Model                   ! KRAKEN or KRAKENC
  COMPLEX,           INTENT( OUT) :: PhiR( Nrd )
  LOGICAL                         :: TufLuk 
  INTEGER                         :: iMat
  COMPLEX                         :: Phi( NMat ), gammaT = 0.0, gammaB = 0.0
  COMPLEX  (KIND=8)               :: gamma2

  READ( ModeFile, REC = IRecProfile + 1 + Mode ) ( Phi( iMat ), iMat = 1, NMat ) 

  ! Is there an elastic medium in the problem?
  TufLuk = .FALSE. 
  IF ( ANY( Material( 1 : NMedia ) == 'ELASTIC' ) ) TufLuk = .TRUE.

  ! Extract the component specified by 'Comp'
  IF ( TufLuk ) CALL EXTRACT( Phi, N, Material, NMedia, Comp ) 

  ! Extract values at receiver depths
  ! n.b. should be using real( k(mode) ) for KRAKEN
  IF ( BCTop == 'A' ) THEN
     SELECT CASE( Model )
     CASE ( 'KRAKENC' )
        gamma2 = k( Mode ) ** 2 - kTop2
     CASE DEFAULT   ! KRAKEN
        gamma2 = REAL( k( Mode ) ) ** 2 - kTop2
     END SELECT
     gammaT = CMPLX( PekerisRoot( gamma2 ) )
  END IF

  IF ( BCBot == 'A' ) THEN 
     SELECT CASE( Model )
     CASE ( 'KRAKENC' )
        gamma2 =       k( Mode )   ** 2 - kBot2
     CASE DEFAULT   ! KRAKEN
        gamma2 = REAL( k( Mode ) ) ** 2 - kBot2
     END SELECT
     gammaB = CMPLX( PekerisRoot( gamma2 ) )
  END IF
  
  RcvrDepth: DO ir = 1, Nrd 
     IF ( rd( ir ) < DepthT ) THEN      ! Receiver in upper halfspace
        PhiR( ir ) = Phi( 1    ) * EXP( -gammaT * ( DepthT - rd( ir  ) ) )

     ELSE IF ( rd( ir ) > DepthB ) THEN ! Receiver in lower halfspace
        PhiR( ir ) = Phi( NTot ) * EXP( -gammaB * ( rd( ir ) - DepthB ) )

     ELSE IF ( NTot > 1 ) THEN 
        iz         = ird( ir )
        PhiR( ir ) = Phi( iz ) + W( ir ) * ( Phi( iz + 1 ) - Phi( iz ) )

     ELSE                               ! mode is tabulated at only one depth
        iz = ird( ir ) 
        PhiR( ir ) = Phi( iz ) 
     ENDIF

  END DO RcvrDepth

END SUBROUTINE ReadOneMode

!**********************************************************************

SUBROUTINE Extract( Phi, N, Material, NMedia, Comp ) 

  ! For elastic problems where a stress-displacement vector is output,
  ! extracts the desired component                                    

  INTEGER,           INTENT( IN    ) :: N( * ), NMedia
  CHARACTER (LEN=1), INTENT( IN    ) :: Comp
  CHARACTER (LEN=8), INTENT( IN    ) :: Material( * )
  COMPLEX,           INTENT( INOUT ) :: Phi( * ) 
  INTEGER                            :: i, j = 1, k = 1, Medium

  MediumLoop: DO Medium = 1, NMedia 

     DO i = 1, N( Medium ) + 1
        MaterialType: SELECT CASE ( Material( Medium ) )
        CASE ( 'ACOUSTIC' )
           Phi( j ) = Phi( k )
           k = k + 1
        CASE ( 'ELASTIC' )
           Component: SELECT CASE ( Comp )
           CASE ( 'H' )   ! horizontal displacement
              Phi( j ) = Phi( k     )
           CASE ( 'V' )   ! vertical   displacement
              Phi( j ) = Phi( k + 1 )
           CASE ( 'T' )   ! tangential stress
              Phi( j ) = Phi( k + 2 )
           CASE ( 'N' )   ! normal     stress
              Phi( j ) = Phi( k + 3 )
           END SELECT Component

           k = k + 4
        END SELECT MaterialType
        j = j + 1 
     END DO   ! next point in medium

  END DO MediumLoop

END SUBROUTINE Extract
END MODULE ReadModes
