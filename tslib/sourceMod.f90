MODULE SourceMod

   USE FatalError

   IMPLICIT NONE
   INTEGER, PARAMETER :: PRTFIL = 6, STSFIL = 10
   INTEGER, PARAMETER :: MaxNt = 10000000, MaxNST = 10000000
   INTEGER, PRIVATE   :: it_Src = 1, it, is ! it_Src used to store the current step into the source time series

CONTAINS

  SUBROUTINE SOURCE( T, S, SD, Nsd, Ntout, omega, fmin, fmax, Pulse, PulseTitle, IniFlag )

    ! Evaluate the source time series (STS) at the time vector T
    ! This should be converted to a module

    ! stemp is a temporary workspace

    ! IniFlag, used to control filtering, has to be controlled on the outside
    ! since SPARC changes the filtering as a function of wavenumber.

    LOGICAL,   INTENT( INOUT ) :: IniFlag
    INTEGER,   INTENT( IN    ) :: Ntout
    INTEGER,   INTENT( INOUT ) :: Nsd
    REAL,      INTENT( IN    ) :: T( * ), omega, fmin, fmax
    REAL,      INTENT( OUT   ) :: SD( Nsd )
    CHARACTER, INTENT( IN    ) :: Pulse*4
    CHARACTER, INTENT( OUT   ) :: PulseTitle*( * )

    INTEGER            :: Nt
    REAL               :: TF( MaxNt ), stempR( MaxNt ), deltat, TStart
    COMPLEX            :: S( Nsd, * ), stemp( MaxNt )
    EQUIVALENCE ( stempR, stemp )
    COMPLEX, ALLOCATABLE :: SF( :, : )

    SAVE TF, SF, Nt

    ! If first call, tabulate the time series

    IF ( IniFlag ) THEN
       Nt = MaxNST / NSd
       IF ( .NOT. ALLOCATED( SF ) ) ALLOCATE( SF( Nsd, Nt ) )

       IF ( Pulse( 1 : 1 ) == 'F' .OR. Pulse( 1 : 1 ) == 'B' ) THEN
          ! From a file
          CALL ReadSTS( Pulse, PulseTitle, SD, Nsd, TF, stempR, SF, Nt )
          deltat = TF( 2 ) - TF( 1 )
          !!! should check that the time series is evenly sampled, i.e. deltat constant
          !!! affects the filtering
       ELSE
          ! Use one of the canned wavelets
          Nsd    = 1
          deltat = 0.2 / omega
          TStart = -200.0 * deltat
          Nt     = 1024

          DO it = 1, Nt
             TF( it ) = TStart + ( it - 1 ) * deltat
             CALL CANS( TF( it ), omega, Pulse, SF, Nsd, it, PulseTitle )
          END DO
       ENDIF

       ! Filter the time series

       DO is = 1, Nsd
          stemp( 1 : Nt ) = SF( is, 1 : Nt ) !  Copy time series into stemp         

          IF ( Pulse( 4 : 4 ) /= 'N' ) CALL BandPass( stemp, Nt, deltat, fmin, fmax )   ! Band-pass filter
          IF ( Pulse( 2 : 2 ) == 'H' ) CALL PreEnv(   stemp, Nt )   ! Form the pre-envelope
          IF ( Pulse( 2 : 2 ) == 'Q' ) CALL Hilbert(  stemp, Nt )   ! Hilbert transform

          SF( is, 1 : Nt ) = stemp( 1 : Nt )   ! Copy time series back
       END DO   ! next source depth

       IniFlag = .FALSE.
    ENDIF

    CALL EvaluateSTS( T, Pulse, TF, SF, Nsd, Nt, S, Ntout )

  END SUBROUTINE SOURCE

  !**********************************************************************C

  SUBROUTINE ReadSTS( Pulse, PulseTitle, SD, Nsd, TF, stemp, SF, Nt )

    ! Time series from file
    ! stemp is used to hold a temporary real vector

    INTEGER,   INTENT( INOUT ) :: NSd, Nt
    REAL,      INTENT( OUT   ) :: SD( NSd ), TF( * )
    COMPLEX,   INTENT( OUT   ) :: SF( Nsd, * )
    CHARACTER, INTENT( IN    ) :: Pulse*4
    CHARACTER, INTENT( OUT )   :: PulseTitle*( * )

    REAL :: stemp( Nsd ), tempV( Nsd ), temp, TMax

    OPEN( FILE = 'STSFIL', UNIT = STSFIL, STATUS = 'OLD', FORM = 'FORMATTED' )

    READ( STSFIL, * ) PulseTitle
    READ( STSFIL, * ) Nsd, SD   ! Number of source depths, vector of source depths

    ! Read in time series data

    Nt = 0
    DO it = 1, MaxNt

       READ( STSFIL, *, END = 2000 ) TF( it ), stemp

       IF ( it > 1 ) THEN
          IF ( TF( it ) < TF( it - 1 ) ) CALL ERROUT( 'SOURCE:SFILE', 'Time series not ordered in time' )
       ENDIF

       SF( :, it ) = stemp
       Nt = Nt + 1

       IF ( Nt * Nsd > MaxNST ) CALL ERROUT( 'SOURCE:SFILE', 'Too many time series points' )
    END DO

    ! normal exit is by EOF in file, not by falling through loop
    CALL ERROUT( 'SOURCE:SFILE', 'Too many time series points' )

    ! Time reversal
    ! note if Nt is odd, middle point doesn't move

2000 IF ( Pulse( 1 : 1 ) == 'B' ) THEN
       TMax = TF( Nt )
       DO it = 1, Nt / 2
          tempV                = REAL( SF( :, it ) )
          SF( :, it )          = SF( :, Nt - it + 1 )
          SF( :, Nt - it + 1 ) = tempV

          temp                  = TMax - TF( it )
          TF( it )              = TMax - TF( Nt - it + 1 )
          TF( Nt - it + 1 )     = temp
       END DO
    ENDIF

    CLOSE( STSFIL )

  END SUBROUTINE ReadSTS

  !**********************************************************************C

  SUBROUTINE EvaluateSTS( T, Pulse, TF, SF, Nsd, Nt, S, Ntout )

    ! Evaluate the source function by linear interpolation between time samples

    INTEGER,   INTENT( IN  ) :: Nsd, Nt, Ntout
    REAL,      INTENT( IN  ) :: T( * ), TF( * )
    COMPLEX,   INTENT( IN  ) :: SF( Nsd, * )
    COMPLEX,   INTENT( OUT ) :: S( Nsd, * )
    CHARACTER, INTENT( IN  ) :: Pulse*4   ! controls a pulse sign-flip

    INTEGER :: itout
    REAL    :: W

    ! it_Src is a saved variable containing the index into the source time series at the current step

    DO itout = 1, Ntout

       ! Necessary to advance time window?

       IF ( T( itout ) < TF( it_Src ) ) it_Src = 1 ! Case of T reset to origin for new march

       !!! It is possible that the following accesses TF( Nt + 1 )
       DO WHILE ( T( itout ) > TF( it_Src + 1 ) .AND. it_Src < Nt - 1 )
          it_Src = it_Src + 1
       END DO

       ! Linear interpolation in time

       S( :, itout ) = 0.0   ! in case T outside time window

       IF ( T( itout ) >= TF( 1  ) .AND. T( itout ) <= TF( Nt ) ) THEN
          W = ( T( itout ) - TF( it_Src ) ) / ( TF( it_Src + 1 ) - TF( it_Src ) )
          S( :, itout ) = SF( :, it_Src ) + W * ( SF( :, it_Src + 1 ) - SF( :, it_Src ) )

          IF ( Pulse( 3 : 3 ) == '-' ) S( :, itout ) = -S( :, itout )
       ENDIF

    END DO   ! next itout

  END SUBROUTINE EvaluateSTS

END MODULE SourceMod
