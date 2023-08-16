MODULE RefCoef

  ! reflection coefficient data

  USE FatalError
  SAVE
  INTEGER, PARAMETER            :: BRCFile = 31, TRCFile = 32, IRCFile = 12
  INTEGER                       :: NBotPts, NTopPts
  INTEGER                       :: NkTab
  INTEGER,          ALLOCATABLE :: iTab( : )
  REAL    (KIND=8), ALLOCATABLE :: xTab( : )
  COMPLEX (KIND=8), ALLOCATABLE :: fTab( : ), gTab( : )

  TYPE ReflectionCoef
      REAL(KIND=8) :: theta, R, phi
  END TYPE

  TYPE(ReflectionCoef), ALLOCATABLE :: RBot( : ), RTop( : )

CONTAINS

  SUBROUTINE ReadReflectionCoefficient( FileRoot, BotRC, TopRC, PRTFile )

    ! Optionally read in reflection coefficient for Top or Bottom boundary

    USE MathConstants
    IMPLICIT NONE
    INTEGER,            INTENT( IN ) :: PRTFile         ! unit number for print file
    CHARACTER (LEN=1 ), INTENT( IN ) :: BotRC, TopRC    ! flag set to 'F' if refl. coef. is to be read from a File
    CHARACTER (LEN=80), INTENT( IN ) :: FileRoot
    INTEGER            :: itheta, ik, IOStat, iAllocStat
    REAL  ( KIND = 8 ) :: freq
    CHARACTER (LEN=80) :: Title2

    IF ( BotRC == 'F' ) THEN
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Using tabulated bottom reflection coef.'
       OPEN( FilE = TRIM( FileRoot ) // '.brc', UNIT = BRCFile, STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
       IF ( IOStat /= 0 ) THEN
          WRITE( PRTFile, * ) 'BRCFile = ', TRIM( FileRoot ) // '.brc'
          CALL ERROUT( 'ReadReflectionCoefficient', 'Unable to open Bottom Reflection Coefficient file' )
       END IF

       READ(  BRCFile, * ) NBotPts
       WRITE( PRTFile, * ) 'Number of points in bottom reflection coefficient = ', NBotPts

       IF ( ALLOCATED( RBot ) ) DEALLOCATE( RBot )
       ALLOCATE( RBot( NBotPts ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 'ReadReflectionCoefficient', 'Insufficient memory for bot. refl. coef.: reduce # points'  )

       READ(  BRCFile, * ) ( RBot( itheta ), itheta = 1, NBotPts )
       CLOSE( BRCFile )
       RBot%phi = DegRad * RBot%phi   ! convert to radians

    ELSE   ! should allocate something anyway, since variable is passed
       IF ( ALLOCATED( RBot ) ) DEALLOCATE( RBot )
       ALLOCATE(  RBot( 1 ), Stat = IAllocStat )
    ENDIF

    ! Optionally read in top reflection coefficient

    IF ( TopRC == 'F' ) THEN
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Using tabulated top    reflection coef.'
       OPEN( FILE = TRIM( FileRoot ) // '.trc', UNIT = TRCFile, STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
       IF ( IOStat /= 0 ) THEN
          WRITE( PRTFile, * ) 'TRCFile = ', TRIM( FileRoot ) // '.trc'
          CALL ERROUT( 'ReadReflectionCoefficient', 'Unable to open Top Reflection Coefficient file' )
       END IF

       READ(  TRCFile, * ) NTopPts
       WRITE( PRTFile, * ) 'Number of points in top reflection coefficient = ', NTopPts

       IF ( ALLOCATED( RTop ) ) DEALLOCATE( RTop )
       ALLOCATE( RTop( NTopPts ), Stat = IAllocStat )
       IF ( iAllocStat /= 0 ) &
          CALL ERROUT( 'ReadReflectionCoefficient', 'Insufficient memory for top refl. coef.: reduce # points'  )

       READ(  TRCFile, * ) ( RTop( itheta ), itheta = 1, NTopPts )
       CLOSE( TRCFile )
       RTop%phi = DegRad *  RTop%phi   ! convert to radians
    ELSE   ! should allocate something anyway, since variable is passed
       IF ( ALLOCATED( RTop ) ) DEALLOCATE( RTop )
       ALLOCATE( RTop( 1 ), Stat = iAllocStat )
    ENDIF

    ! Optionally read in internal reflection coefficient data

    IF ( BotRC == 'P' ) THEN
       WRITE( PRTFile, * ) 'Reading precalculated refl. coeff. table'
       OPEN( FILE = TRIM( FileRoot ) // '.irc', UNIT = IRCFile, STATUS = 'OLD', IOSTAT = IOStat, ACTION = 'READ' )
       IF ( IOStat /= 0 ) CALL ERROUT( 'ReadReflectionCoefficient', &
                                       'Unable to open Internal Reflection Coefficient file' )

       READ(  IRCFile, * ) Title2, freq
       READ(  IRCFile, * ) NkTab
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of points in internal reflection coefficient = ', NkTab

       IF ( ALLOCATED( xTab ) )  DEALLOCATE( xTab, fTab, gTab, iTab )
       ALLOCATE( xTab( NkTab ), fTab( NkTab ), gTab( NkTab ), iTab( NkTab ), Stat = iAllocStat )
       IF ( iAllocStat /= 0 ) CALL ERROUT( 'ReadReflectionCoefficient', 'Too many points in reflection coefficient' )

       READ( IRCFile, FMT = "( 5G15.7, I5 )" ) ( xTab( ik ), fTab( ik ), gTab( ik ), iTab( ik ), ik = 1, NkTab )
       CLOSE( IRCFile )
    ENDIF

  END SUBROUTINE ReadReflectionCoefficient

  SUBROUTINE InterpolateReflectionCoefficient( RInt, R, NPts, PRTFile )

    ! Given an angle RInt%ThetaInt, returns the magnitude and
    ! phase of the reflection coefficient (RInt%R, RInt%phi).

    ! Uses linear interpolation using the two nearest abscissas
    ! Assumes phi has been unwrapped so that it varies smoothly.
    ! I tried modifying it to allow a complex angle of incidence but
    ! stopped when I realized I needed to fuss with a complex atan2 routine

    IMPLICIT NONE
    INTEGER,              INTENT( IN    ) :: PRTFile, NPts   ! unit number of print file, # pts in refl. coef.
    TYPE(ReflectionCoef), INTENT( IN    ) :: R( NPts )       ! Reflection coefficient table
    TYPE(ReflectionCoef), INTENT( INOUT ) :: RInt            ! interpolated value of refl. coef.
    INTEGER         :: iLeft, iRight, iMid
    REAL ( KIND=8 ) :: alpha, Thetaintr

    iLeft  = 1
    iRight = NPts

    thetaIntr = REAL( RInt%Theta )   ! This should be unnecessary? probably used when I was doing complex angles

    ! Three cases: ThetaInt left, in, or right of tabulated interval

    IF     ( thetaIntr < R( iLeft  )%theta ) THEN
       !iRight = 2
       RInt%R   = 0.0     ! R( iLeft  )%R
       RInt%phi = 0.0     ! R( iLeft  )%phi
       WRITE( PRTFile, * ) 'Warning in InterpolateReflectionCoefficient : Refl. Coef. being set to 0 outside tabulated domain'
       WRITE( PRTFile, * ) 'angle = ', thetaintr, 'lower limit = ', R( iLeft)%theta

    ELSE IF( thetaIntr > R( iRight )%theta ) THEN
       !iLeft = NPts - 1
       RInt%R   = 0.0     ! R( iRight )%R
       RInt%phi = 0.0     ! R( iRight )%phi
       ! WRITE( PRTFile, * ) 'Warning in InterpolateReflectionCoefficient : Refl. Coef. being set to 0 outside tabulated domain'
       ! WRITE( PRTFile, * ) 'angle = ', ThetaIntr, 'upper limit = ', R( iRight )%theta

    ELSE
       ! Search for bracketing abscissas: Log2( NPts ) stabs required for a bracket

       DO WHILE ( iLeft /= iRight - 1 )
          iMid = ( iLeft + iRight ) / 2
          IF ( R( iMid )%theta > thetaIntr ) THEN
             iRight = iMid
          ELSE
             iLeft  = iMid
          ENDIF
       END DO

       ! Linear interpolation for reflection coef

       alpha    = ( RInt%theta - R( iLeft )%theta ) / ( R( iRight )%theta - R( iLeft )%theta )
       RInt%R   = ( 1 - alpha ) * R( iLeft )%R   + alpha * R( iRight )%R
       RInt%phi = ( 1 - alpha ) * R( iLeft )%phi + alpha * R( iRight )%phi

    ENDIF

  END SUBROUTINE InterpolateReflectionCoefficient

  SUBROUTINE InterpolateIRC( x, f, g, iPower, xTab, fTab, gTab, iTab, NkTab )

    ! Internal reflection coefficient interpolator.
    ! Returns f, g, iPower for given x using tabulated values.
    ! Uses polynomial interpolation to approximate the function between the tabulated values

    USE PolyMod
    IMPLICIT NONE

    INTEGER, PARAMETER :: N = 3                    ! order of the polynomial for interpolation
    INTEGER,            INTENT( IN  ) :: NkTab, iTab( NkTab )
    REAL    ( KIND=8 ), INTENT( IN  ) :: xTab( NkTab )
    COMPLEX ( KIND=8 ), INTENT( IN  ) :: fTab( NkTab ), gTab( NkTab ), x
    INTEGER,            INTENT( OUT ) :: iPower
    COMPLEX ( KIND=8 ), INTENT( OUT ) :: f, g
    INTEGER            :: i, j, iDel, iLeft, iMid, iRight, NAct
    REAL    ( KIND=8 ) :: xReal
    COMPLEX ( KIND=8 ) :: xT( N ), fT( N ), gT( N )

    xReal  = DBLE( x )        ! taking the real part
    iLeft  = 1
    iRight = NkTab

    ! Three cases: x left, in, or right of tabulated interval

    IF (     xReal < xTab( iLeft  ) ) THEN
       f      = fTab( iLeft  )
       g      = gTab( iLeft  )
       iPower = iTab( iLeft  )
    ELSE IF( xReal > xTab( iRight ) ) THEN
       f      = fTab( iRight )
       g      = gTab( iRight )
       iPower = iTab( iRight )
    ELSE

       ! Search for bracketing abscissas:
       ! Log base 2 (NPts) stabs required for a bracket

       DO WHILE ( iLeft /= iRight-1 )
          iMid = ( iLeft + iRight )/2
          IF ( xTab( iMid ) > xReal ) THEN
             iRight = iMid
          ELSE
             iLeft  = iMid
          ENDIF
       END DO

       ! Extract the subset for interpolation and scale

       iLeft  = MAX( iLeft  - ( N - 2 ) / 2, 1     )
       iRight = MiN( iRight + ( N - 1 ) / 2, NkTab )

       NAct = iRight - iLeft + 1
       DO i = 1, NAct
          j       = i + iLeft - 1
          iDel    = iTab( j ) - iTab( iLeft )
          xT( i ) = xTab( j )
          fT( i ) = fTab( j ) * 10.0D0 ** iDel
          gT( i ) = gTab( j ) * 10.0D0 ** iDel
       END DO

       ! Construct the interpolate

       f      = Poly( x, xT, fT, NAct )
       g      = Poly( x, xT, gT, NAct )
       iPower = iTab( iLeft )
    ENDIF

  END SUBROUTINE InterpolateIRC
END MODULE RefCoef
