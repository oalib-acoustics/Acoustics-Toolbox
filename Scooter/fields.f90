PROGRAM FIELDS

  ! Note: The Matlab version of this (fieldsco.m) is preferred unless memory efficiency is a problem
  ! Also, this has not been updated to process multiple frequencies

  ! Compute pressure from Green's function.

  ! The transform parameters are chosen to satisfy certain sampling requirements.

  ! It is assumed that the user chose the maximum possible 
  ! spacing in computing g(k) and so Delta_k is not allowed to be any larger.

  ! The other part of the kernel is the exp(ikr) term.
  ! This must be sampled on the k-axis finely enough to have about 6 pts per wavelength.

  ! Delta_r implies kmax. kmax is the upper limit of integration and beyond the last point computed.
  ! This assumes the user had chosen kmax as small as possible while still covering the support of G(k).

  USE BeamPattern
  USE SourceReceiverPositions
  IMPLICIT NONE

  SAVE
  INTEGER, PARAMETER   :: FLPFile = 5, PRTFile = 6, GRNFile = 20, SHDFile = 25
  REAL,    PARAMETER   :: pi = 3.14159265, RadDeg = 180.0 / pi
  COMPLEX (KIND=8), PARAMETER :: i  = ( 0.0, 1.0 )
  INTEGER              :: Nrr
  INTEGER              :: IAllocStat, Nk, NrrLast, Nt, Nt2
  REAL                 :: Atten, AttInt, Freq, kmax, Delta_k, Delta_kInterp, &
                          Rmin, Rmax, Delta_r
  CHARACTER   (LEN=80) :: PlotTitle, FileRoot
  CHARACTER   (LEN= 4) :: Option
  REAL,    ALLOCATABLE :: k( : ), kInterp( : ), cVec( : )
  COMPLEX, ALLOCATABLE :: G( : ), Ginterp( : ), P( : ), Temp( : )

  INTEGER :: IRatio, IRatioDelta_r, J, IOStat

  ! get the file root for naming all input and output files
  ! should add some checks here ...

  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )

  ! Open the print file
  OPEN( UNIT = PRTFile, FILE = 'fields.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )

  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) '__________________________________________________'
  WRITE( PRTFile, * ) 'Running FIELDS'
  WRITE( PRTFile, * ) 'Transforms Green''s function FILE, producing pressure'
  WRITE( PRTFile, * )

  ! open the field paramaters file (FLPFile)
  OPEN( FILE = TRIM( FileRoot ) // '.flp', UNIT = FLPFile, STATUS = 'OLD', FORM = 'FORMATTED', IOSTAT = IOStat, ACTION = 'READ' )

  IF ( IOStat /= 0 ) THEN
     WRITE( PRTFile, * ) 'FLPFile = ', TRIM( FileRoot ) // '.flp'
     CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Unable to open FLPFile' )
  END IF

  ! Begin by reading in the data
  READ( FLPFile, * ) Option
  SELECT CASE ( Option( 1 : 1 ) )
  CASE ( 'R' )
     WRITE( PRTFile, * ) 'Cylindrical coordinates/point source'
  CASE ( 'X' )
     WRITE( PRTFile, * ) 'Cartesian coordinates/line source'
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Option( 1 : 1 ) should select Line vs. point source' )
  END SELECT

  SELECT CASE ( Option( 2 : 2 ) )
  CASE ( 'P' )
     WRITE( PRTFile, * ) 'Using only positive part of wavenumber spectrum'
  CASE ( 'N' )
     WRITE( PRTFile, * ) 'Using only negative part of wavenumber spectrum'
  CASE ( 'B' )
     WRITE( PRTFile, * ) 'Using both positive and negative part of wavenumber spectrum'
  CASE DEFAULT
     CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Option(2:2) should select positive, negative, or both parts of the spectrum' )
  END SELECT

  SELECT CASE ( Option( 3 : 3 ) )
  CASE ( 'O' )
     WRITE( PRTFile, * ) 'Polynomial interpolation'
  CASE ( 'P' )
     WRITE( PRTFile, * ) 'Pade interpolation'
  CASE DEFAULT
     Option( 3 : 3 ) = 'O'   ! use Polynomial interpolation as default
     !CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Option(3:3) should select Polynomial (''O'') or Pade (''P'') interpolation' )
  END SELECT

  ! optionally read in a source beam pattern
  SBPFlag = Option( 4 : 4 )
  CALL ReadPAT( FileRoot, PRTFile )  ! Source Beam Pattern

  CALL ReadRcvrRanges   ! Read receiver ranges, Rr
  Nrr     = Pos%NRr
  Rmin    = Pos%Rr( 1   )
  Rmax    = Pos%Rr( Nrr )
  Delta_r = ( Rmax - Rmin ) / ( NRr - 1 )

  CALL ReadHeader   ! Read the header records from GRNFile

  IF ( Nfreq > 1 ) THEN
     CALL ERROUT( PRTFile, 'F', 'FIELDS', 'cannot handle broadband runs; use the Matlab version FIELDSCO.m' )
  END IF

  ! Set up for transform: need Delta_k, NT
  ! Delta_k is what scooter used; Delta_kInterp is for interpolation
  ! If Delta_r is too big, take submultiple

  kmax = k( 1 ) + 2.0 * pi / Delta_r
  IRatioDelta_r = 1
  IF ( kmax < k( Nk ) ) THEN
     IRatioDelta_r = INT( ( k( Nk ) - k( 1 ) ) / ( kmax - k( 1 ) ) ) + 1
     Delta_r       = Delta_r / IRatioDelta_r
     Nrr           = IRatioDelta_r * ( NRr - 1 ) + NRr
     kmax          = k( 1 ) + 2 * pi / Delta_r
     WRITE( PRTFile, * ) 'Number or ranges, Nrr, increased so that wavenumber limit exceeds kmax used by SCOOTER', Nrr
  END IF

  ! Compute NT based on Delta_k (= 1/Rmax) sampling requirement
  NT2 = NINT( Rmax * ( kmax - k( 1 ) ) / ( 2 * pi ) )       ! Delta_k limit
  IF ( NT2 > Nrr ) THEN
     WRITE( PRTFile, * ) 'NT bumped to', NT2
     WRITE( PRTFile, * ) 'Thus we are zero filling the wavenumber spectrum to effectively interpolate on a finer grid'
  END IF
  NT = MAX( Nrr, NT2 )

  ! bump Nt if necessary to make sure Delta_kInterp is not coarser than Delta_k grid
  Delta_k       = ( k( Nk ) - k( 1 ) ) / Nk;
  Delta_kInterp = 2 * pi / ( Nt * Delta_r );
  IF ( Delta_kInterp > Delta_k ) THEN
     IRatio = NINT( Delta_kInterp / Delta_k );
     Nt     = IRatio * Nt;
     WRITE( PRTFile, * ) 'Transform size, Nt, bumped to ensure Delta_k sampling is fine enough', Nt
  END IF

  ! Parameters when N must be a power of 2
  NT = 2 ** ( INT( LOG10( REAL( NT ) ) / 0.301 ) + 1 )
  Delta_kInterp = 2 * pi / ( NT * Delta_r )

  WRITE( PRTFile, * ) 'NT  used = ', NT
  WRITE( PRTFile, * ) 'Nrr used = ', Nrr
  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) 'Delta_k       = ', Delta_k
  WRITE( PRTFile, * ) 'Delta_kInterp = ', Delta_kInterp
  WRITE( PRTFile, * ) 'Delta_r       = ', Delta_r
  WRITE( PRTFile, * )

  ALLOCATE( kInterp( NT ), Ginterp( NT ), P( NT ), Temp( 2 * NT ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Insufficient memory to allocate kInterp, Ginterp, ...' )

  ! Compute optimal stabilizing attenuation.
  ! Atten is the value used by SCOOTER
  ! AttInt is the optimal value for the interpolated k-values

  AttInt = ( Delta_kInterp / Delta_k ) * Atten

  IF ( PlotTitle( 1 : 5 ) == 'SPARC' ) THEN
     WRITE( PRTFile, * ) 'SPARC RUN: AttInt SET TO ZERO'
     AttInt = 0.0
     IF ( Atten /= 0.0 ) CALL ERROUT( PRTFile, 'F', 'FIELDS', 'Stabilizing attenuation must vanish' )
  ENDIF

  kInterp = k( 1 ) + [ ( J, J = 0, Nt - 1 ) ] * Delta_kInterp   ! Set up vector of kInterp points
  kmax    = kInterp( Nt ) + Delta_kInterp   ! because the fft goes from 0 to kmax but G(kmax) is not computed
  Delta_r = 2 * pi / ( kmax - kInterp( 1 ) )
  NRr     = Nt
  NRrLast = min( NINT( ( Rmax - Rmin ) / Delta_r ) + 1, Nt )

  ! set up vector of range points (subsampled to satisfy user request)
  Pos%Nrr = ( NrrLast - 1 ) / IratioDelta_r + 1
  DEALLOCATE( Pos%Rr )
  ALLOCATE( Pos%Rr( Pos%Nrr ) )
  Pos%Rr = Rmin + [ ( J, J = 0, Pos%NRr - 1 ) ] * IratioDelta_r * Delta_r

  ! Construct shade file
  CALL SHADE( IratioDelta_r )

  WRITE( PRTfile, * ) 'Fields completed successfully'

CONTAINS

  !**********************************************************************!

  SUBROUTINE ReadHeader

    ! Routine to read header from disk file
    ! This routine is essentially the same as at/misc/ReadHeader except that the range vector is replaced
    ! by a wavenumber vector.

    USE SourceReceiverPositions

    INTEGER :: IOStat, LRecL

    OPEN( FILE = TRIM( FileRoot ) // '.grn', UNIT = GRNFile, STATUS = 'OLD', ACCESS = 'DIRECT', &
         FORM = 'UNFORMATTED', RECL = 10, IOSTAT = IOStat, ACTION = 'READ' )
    IF ( IOStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'FIELDS:RDHead', 'Unable to open GRNFile' )

    READ( GRNFile, REC = 1 ) LRECL
    CLOSE( GRNFile )
    OPEN( FILE = TRIM( FileRoot ) // '.grn', UNIT = GRNFile, STATUS = 'OLD', ACCESS = 'DIRECT', &
         FORM = 'UNFORMATTED', RECL = 4 * LRECL, ACTION = 'READ' )

    ! Read data
    READ( GRNFile, REC = 1  ) LRECL, PlotTitle
    !READ( GRNFile, REC = 2 ) PlotType, XS, YS
    READ( GRNFile, REC = 3  ) Nfreq, Pos%Ntheta, Pos%NSx, Pos%NSy, Pos%NSz, Pos%NRz, Pos%Nrr, atten
    Nk = Pos%Nrr

    ALLOCATE( k( Nk ), G( Nk ), cVec( Nk ), Stat = IAllocStat )
    ALLOCATE( FreqVec( Nfreq ), Pos%sx( Pos%NSx ), Pos%sy( Pos%NSy ), Pos%sz( Pos%NSz ), &
         Pos%rz( Pos%NRz ), Pos%theta( Pos%Ntheta ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadHeader', 'Too many source/receiver combinations' )

    READ( GRNFile, REC = 4  ) FreqVec
    READ( GRNFile, REC = 5  ) Pos%theta
    READ( GRNFile, REC = 6  ) Pos%Sx
    READ( GRNFile, REC = 7  ) Pos%Sy
    READ( GRNFile, REC = 8  ) Pos%Sz
    READ( GRNFile, REC = 9  ) Pos%Rz
    READ( GRNFile, REC = 10 ) cVec   ! Pos%r contains the phase speeds

    Freq    = FreqVec( 1 )
    k       = 2 * pi * Freq / cVec
    Delta_k = k( 2 ) - k( 1 )

  END SUBROUTINE ReadHeader

  !**********************************************************************!

  SUBROUTINE SHADE( IratioDelta_r )

    ! Performs the transforms to convert the Green's function file to a shade file
    ! Expects
    !    k, Nk: The wavenumber data
    !    sz, rz: Source depth, receiver depth data
    ! Returns
    !    Nothing

    USE SourceReceiverPositions
    USE BeamPattern
    USE interpolation
    INTEGER, INTENT( IN ) :: IratioDelta_r
    INTEGER               :: Isz, Irz, ISR, IRec
    CHARACTER ( LEN=10 )  :: PlotType = '          '
    REAL      ( KIND=8 )  :: kz2( Nk ), thetaT( Nk ), S( Nk )
    REAL      ( KIND=8 )  :: omega, c0

    Rmin = Pos%Rr( 1 )
    CALL WriteHeader( TRIM( FileRoot ) // '.shd', PlotTitle, Atten, PlotType )

    SourceDepth: DO Isz = 1, Pos%NSz
       WRITE( PRTFile, * ) 'Transform for source depth: ', Pos%sz( Isz )

       RcvrDepth: DO Irz  = 1, Pos%NRz
          ISR  = ( Isz - 1 ) * Pos%NRz + Irz   ! Index of source/receiver
          IRec = 10 + ISR
          READ( GRNFile, REC = IRec ) G( 1 : Nk )

          ! apply the source beam pattern
          IF ( SBPFlag == '*' ) THEN
             c0    = 1500;   ! reference sound speed, should be speed at the source depth
             omega = 2 * pi * freq
             kz2   = omega ** 2 / c0 ** 2 - k( 1 : Nk ) ** 2      ! vertical wavenumber squared
             WHERE ( kz2 < 0 ) kz2 = 0                            ! remove negative values

             thetaT = RadDeg * ATAN( SQRT( kz2 ) / k( 1 : Nk ) )  ! calculate the angle in degrees
             CALL interp1( SrcBmPat( :, 1 ), SrcBmPat( :, 2 ), thetaT, S )
             
             G( 1 : Nk ) = G( 1 : Nk ) * REAL( S )                ! apply the shading
          END IF

          CALL InterpolateG( k, Atten, AttInt, G, Nk, kInterp, Ginterp, NT, Option )

          ! Evaluate either Fourier or Hankel transform
          IF ( Option( 1 : 1 ) == 'X' ) THEN
             CALL FourierTransform( NT, kInterp( 1 ), Delta_kInterp, AttInt, Rmin, NRr, Delta_r, Ginterp, Temp, P )
          ELSE
             CALL HankelTransform(  NT, kInterp( 1 ), Delta_kInterp, AttInt, Rmin, NRr, Delta_r, Ginterp, Temp, P, Option )
          ENDIF

          WRITE( SHDFile, REC = IRec ) P( 1 : NRrLast : IratioDelta_r )   ! Write out the field

       END DO RcvrDepth
    END DO SourceDepth

  END SUBROUTINE SHADE

  !**********************************************************************!

  SUBROUTINE InterpolateG( k, Atten, AttInt, G, Nk, kInterp, Ginterp, NT, Option )

    ! Produces an evenly sampled kernel from input G

    USE PolyMod
    INTEGER, PARAMETER   :: iSize = 3
    COMPLEX, PARAMETER   :: i = ( 0.0, 1.0 )
    INTEGER, INTENT(IN)  :: Nk, Nt
    REAL,    INTENT(IN)  :: k( Nk ), kInterp( Nt ), Atten, AttInt
    COMPLEX, INTENT(IN)  :: G( Nk )
    COMPLEX, INTENT(OUT) :: Ginterp( Nt )
    INTEGER              :: ICenter, IRight, INew, It
    COMPLEX              :: xinterp, x( iSize ), f( iSize ), Pade
    CHARACTER   (LEN=80) :: ErrorMessage = '      '
    CHARACTER   (LEN=3 ) :: Option

    ! Initialize interpolation data
    ICenter = iSize / 2 + 1   ! counter to mark center of abscissas used for interpolation
    IRight  = iSize
    INew    = iSize

    x( 1 : iSize ) = ( k( 1 : iSize ) + i * Atten )**2   ! abscissa
    f( 1 : iSize ) = G( 1 : iSize )                      ! ordinate

    DO It = 1, NT ! Main loop: step through values in kInterp

       ! Desirable/possible to advance interpolation window?
       DO WHILE ( kInterp( It ) > k( ICenter ) .AND. IRight < Nk )
          ICenter = ICenter + 1
          IRight  = IRight  + 1
          INew    = MOD( INew, iSize ) + 1   ! this is the next open slot

          x( INew ) = ( k( IRight ) + i * Atten ) ** 2
          f( INew ) = G( IRight )
       END DO

       ! Interpolate or zero fill
       IF ( kInterp( It ) <= k( Nk ) ) THEN
          xinterp = ( kInterp( It ) + i * AttInt ) ** 2
          IF ( Option( 3 : 3 ) == 'O' ) THEN
             Ginterp( It ) = Poly( xinterp, x, f, iSize )   ! polynomial
          ELSE
             Ginterp( It ) = Pade( xinterp, x, f, iSize, ErrorMessage )   ! Pade
          ENDIF

          IF ( ErrorMessage( 1 : 6 ) /= '      ' ) WRITE( PRTFile, * ) ErrorMessage
       ELSE
          Ginterp( It ) = 0.0
       ENDIF
    END DO

  END SUBROUTINE InterpolateG

END PROGRAM FIELDS
