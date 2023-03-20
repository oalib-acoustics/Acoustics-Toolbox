MODULE sspmod

  ! holds SSP input by user and associated variables

  USE FatalError

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER, PRIVATE :: ENVFile = 5, PRTFile = 6
  INTEGER, PARAMETER          :: MaxSSP = 20001, MaxMedia = 501
  INTEGER, PRIVATE            :: N, iz, ILoc, Lay
  INTEGER                     :: ISSP
  REAL (KIND=8)               :: alphaR = 1500, betaR = 0, alphaI = 0, betaI = 0, rhoR = 1
  REAL (KIND=8), PRIVATE      :: h, z, R

  ! SSP
  TYPE SSPStructure
     INTEGER           :: Loc( MaxMedia ), NPts( MaxMedia ), NMedia
     REAL     (KIND=8) :: z( MaxSSP ), alphaR( MaxSSP ), alphaI( MaxSSP ), rho( MaxSSP ), betaR( MaxSSP ), betaI( MaxSSP )
     REAL     (KIND=8) :: Depth( MaxMedia ), sigma( MaxMedia ), beta( MaxMedia ), fT( MaxMedia )
     COMPLEX  (KIND=8) :: cp( MaxSSP ), cs( MaxSSP ), n2( MaxSSP ),  &
                          cpSpline( 4, MaxSSP ), csSpline( 4, MaxSSP ), rhoSpline( 4, MaxSSP )
     COMPLEX  (KIND=8) :: cpCoef( 4, MaxSSP ), csCoef( 4, MaxSSP ), rhoCoef( 4, MaxSSP ), csWork( 4, MaxSSP )   ! for the PCHIP coefficients
     COMPLEX  (KIND=8) :: rhoT( MaxSSP )   ! temporary values for calling PCHIP (ordinate argument must be complex valued!)
     CHARACTER (LEN=1) :: Type
     CHARACTER (LEN=2) :: AttenUnit
     CHARACTER (LEN=8) :: Material( MaxMedia )
  END TYPE SSPStructure

  TYPE( SSPStructure ) :: SSP

  TYPE HSInfo
     CHARACTER (LEN=1) :: BC                            ! Boundary condition type
     REAL     (KIND=8) :: alphaR, alphaI, betaR, betaI  ! P-wave, S-wave speeds (user units)
     REAL     (KIND=8) :: beta, fT                      ! power law and transition frequency
     COMPLEX  (KIND=8) :: cP, cS                        ! P-wave, S-wave speeds (neper/m loss)
     REAL     (KIND=8) :: rho, BumpDensity, eta, xi     ! density, boss parameters
  END TYPE

  TYPE( HSInfo )       :: HSTop, HSBot

CONTAINS

  SUBROUTINE EvaluateSSP( cP, cS, rho, Medium, N1, freq, Task )

    ! Call the particular SSP routine specified by SSPType
    ! Performs two Tasks:
    !    Task = 'TAB'  then tabulate cP, cS, rho
    !    Task = 'INIT' then initialize
    ! Note that freq is only needed if Task = 'INIT'

    INTEGER,           INTENT(IN)    :: Medium
    INTEGER,           INTENT(INOUT) :: N1
    REAL     (KIND=8), INTENT(OUT)   :: rho( * )
    REAL     (KIND=8), INTENT(IN)    :: freq
    COMPLEX  (KIND=8), INTENT(OUT)   :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)    :: Task
    COMPLEX  (KIND=8)                :: cPT, cST

    SELECT CASE ( SSP%Type )
    CASE ( 'A' )  !  Analytic profile option 
       IF ( Task( 1 : 4 ) == 'INIT' ) THEN
          N1 = 21

          CALL ANALYT( cP, cS, rho, Medium, N1 )
          h = ( SSP%Depth( Medium + 1 ) - SSP%Depth( Medium ) ) / ( N1 - 1 )

          DO iz = 1, N1
             z   = SSP%Depth( Medium ) + ( iz - 1 ) * h
             cPT =  cP( iz )
             cST =  cS( iz )
             WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
                  z,  REAL( cPT ),  REAL( cST ), rho( iz ), AIMAG( cPT ), AIMAG( cST )
          END DO
       ELSE
          CALL ANALYT( cP, cS, rho, Medium, N1 )
       ENDIF
    CASE ( 'N' )  !  N2-linear profile option
       CALL n2Linear( cP, cS, rho, Medium, N1, Task )
    CASE ( 'C' )  !  C-linear profile option 
       CALL cLinear(  cP, cS, rho, Medium, N1, Task )
    CASE ( 'P' )  !  monotone PCHIP ACS profile option 
       CALL cPCHIP(   cP, cS, rho, Medium, N1, Task )
    CASE ( 'S' )  !  Cubic spline profile option 
       CALL cCubic(   cP, cS, rho, Medium, N1, Task )
    CASE DEFAULT  !  Non-existent profile option 
       WRITE( PRTFile, * ) 'Profile option: ', SSP%Type
       CALL ERROUT( 'EvaluateSSP', 'Unknown profile option' )
    END SELECT

    RETURN
  END SUBROUTINE EvaluateSSP

  !**********************************************************************!

  SUBROUTINE n2Linear( cP, cS, rho, Medium, N1, Task )

    ! Tabulate cP, cS, rho for specified Medium
    ! Uses N2-linear segments for P and S-wave speeds
    ! Uses rho-linear segments for density

    INTEGER,           INTENT( INOUT ) :: N1
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task
    COMPLEX  (KIND=8)              :: N2Bot, N2Top

    ! If Task = 'INIT' then this is the first call and SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization
       CALL ReadSSP( Medium, N1 )
    ELSE   ! Task = 'TABULATE'
       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          z = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) z = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot

          DO WHILE ( z > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          iSSP = ILoc + Lay
          R = ( z - SSP%z( iSSP ) ) / ( SSP%z( iSSP + 1 ) - SSP%z( iSSP ) )

          ! P-wave

          N2Top    = 1.0 / SSP%cp( iSSP     )**2
          N2Bot    = 1.0 / SSP%cp( iSSP + 1 )**2
          cP( iz ) = 1.0 / SQRT( ( 1.0 - R ) * N2Top + R * N2Bot )

          ! S-wave

          IF ( SSP%cs( iSSP ) /= 0.0 ) THEN
             N2Top    = 1.0 / SSP%cs( iSSP     )**2
             N2Bot    = 1.0 / SSP%cs( iSSP + 1 )**2
             cS( iz ) = 1.0 / SQRT( ( 1.0 - R ) * N2Top + R * N2Bot )
          ELSE
             cS( iz ) = 0.0
          ENDIF

          rho( iz ) = ( 1.0 - R ) * SSP%rho( iSSP ) + R * SSP%rho( iSSP + 1 )
       END DO

    ENDIF

    RETURN
  END SUBROUTINE n2Linear
  
  !**********************************************************************!
  
  SUBROUTINE cLinear( cP, cS, rho, Medium, N1, Task )

    ! Tabulate cP, cS, rho for specified Medium
    ! Uses c-linear segments for P and S-wave speeds
    ! Uses rho-linear segments for density

    INTEGER,           INTENT( INOUT ) :: N1
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task

    ! If Task = 'INIT' then this is the first call and SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization
       CALL ReadSSP( Medium, N1 )
    ELSE   ! Task = 'TABULATE'
       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          z = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) z = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot

          DO WHILE ( z > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          iSSP = ILoc + Lay
          R = ( z - SSP%z( iSSP ) ) / ( SSP%z( iSSP + 1 ) - SSP%z( iSSP ) )
          cP(  iz ) = ( 1.0 - R ) * SSP%cp(  iSSP ) + R * SSP%cp(  iSSP + 1 )
          cS(  iz ) = ( 1.0 - R ) * SSP%cs(  iSSP ) + R * SSP%cs(  iSSP + 1 )
          rho( iz ) = ( 1.0 - R ) * SSP%rho( iSSP ) + R * SSP%rho( iSSP + 1 )
       END DO
    ENDIF

    RETURN
  END SUBROUTINE cLinear

  !**********************************************************************!
  
  SUBROUTINE cPCHIP( cP, cS, rho, Medium, N1, Task )

    ! Tabulate cP, cS, rho values for the specified Medium
    ! Uses PCHIP segments for P, S-wave speeds and density rho

    USE pchipMod
    INTEGER,           INTENT( INOUT ) :: N1
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task
    REAL     (KIND=8)              :: xt

    ! If Task = 'INIT' then this is the first call and the SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization
       CALL ReadSSP( Medium, N1 )
    ELSE   ! Task = 'TABULATE'
       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          z = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) z = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot
          DO WHILE ( z > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          iSSP = ILoc + Lay
          xt   = z - SSP%z( iSSP )

          cP(  iz ) =       SSP%cpCoef( 1,  iSSP ) &
                        + ( SSP%cpCoef( 2,  iSSP ) &
                        + ( SSP%cpCoef( 3,  iSSP ) &
                        +   SSP%cpCoef( 4,  iSSP ) * xt ) * xt ) * xt
          cS(  iz ) =       SSP%csCoef( 1,  iSSP ) &
                        + ( SSP%csCoef( 2,  iSSP ) &
                        + ( SSP%csCoef( 3,  iSSP ) &
                        +   SSP%csCoef( 4,  iSSP ) * xt ) * xt ) * xt
          rho( iz ) = DBLE( SSP%rhoCoef( 1, iSSP ) &
                        + ( SSP%rhoCoef( 2, iSSP ) &
                        + ( SSP%rhoCoef( 3, iSSP ) &
                        +   SSP%rhoCoef( 4, iSSP ) * xt ) * xt ) * xt )
       END DO
    ENDIF

    RETURN
  END SUBROUTINE cPCHIP

  !**********************************************************************!
  
  SUBROUTINE cCubic( cP, cS, rho, Medium, N1, Task )

    ! Tabulate cP, cS, rho for specified Medium
    ! Uses cubic spline for P, S-wave speeds and density rho

    INTEGER,           INTENT( INOUT ) :: N1
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task
    REAL     (KIND=8)              :: HSPLNE
    COMPLEX  (KIND=8)              :: SPLINE

    ! If Task = 'INIT' then this is the first call and SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization
       CALL ReadSSP( Medium, N1 )
    ELSE   ! Task = 'TABULATE'
       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          z = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) z = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot
          DO WHILE ( z > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          iSSP = ILoc + Lay
          HSPLNE = z - SSP%z( iSSP )

          cP(  iz ) =       SPLINE( SSP%cpSpline(  1, iSSP ), HSPLNE )
          cS(  iz ) =       SPLINE( SSP%csSpline(  1, iSSP ), HSPLNE )
          rho( iz ) = DBLE( SPLINE( SSP%rhoSpline( 1, iSSP ), HSPLNE ) )

       END DO
    ENDIF

    RETURN
  END SUBROUTINE cCubic

!**********************************************************************!

  SUBROUTINE ReadSSP( Medium, N1 )

    ! reads the SSP data from the environmental file for a given medium
    
    INTEGER, INTENT( IN    ) :: Medium
    INTEGER, INTENT( INOUT ) :: N1
    INTEGER                  :: iSSP

    WRITE( PRTFile, * )

    IF ( Medium == 1 ) THEN
       WRITE( PRTFile, * ) 'Sound speed profile:'
       WRITE( PRTFile, "( '      z         alphaR      betaR     rho        alphaI     betaI'    )" )
       WRITE( PRTFile, "( '     (m)         (m/s)      (m/s)   (g/cm^3)      (m/s)     (m/s)', / )" )
    END IF

    SSP%NPts( Medium ) = N1

    ! The variable SSP%Loc( Medium ) points to the starting point for the
    ! data in the arrays z, alpha, beta and rho
    IF ( Medium == 1 ) THEN
       SSP%Loc( Medium ) = 0
    ELSE
       SSP%Loc( Medium ) = SSP%Loc( Medium - 1 ) + SSP%NPts( Medium - 1 )
    ENDIF
    ILoc = SSP%Loc( Medium )

    !  Read in data and convert attenuation to Nepers/m 
    N1 = 1
    DO iSSP = 1, MaxSSP
       iz = SSP%Loc( Medium ) + iSSP

       READ(  ENVFile, *    ) SSP%z( iz ), alphaR, betaR, rhoR, alphaI, betaI
       WRITE( PRTFile, FMT="( F10.2,      3X, 2F10.2,       3X, F6.2, 3X, 2F10.4 )" ) &
                              SSP%z( iz ), alphaR, betaR, rhoR, alphaI, betaI

       ! Verify that the depths are monotonically increasing
       IF ( iSSP > 1 ) THEN
          IF ( SSP%z( iz ) .LE. SSP%z( iz - 1 ) ) THEN
              WRITE( PRTFile, * ) 'Bad depth in SSP: ', SSP%z( iz )
              CALL ERROUT( 'ReadSSP', 'The depths in the SSP must be monotonically increasing' )
          END IF
       END IF

       SSP%alphaR( iz ) = alphaR
       SSP%alphaI( iz ) = alphaI
       SSP%rho(    iz ) = rhoR
       SSP%betaR(  iz ) = betaR
       SSP%betaI(  iz ) = betaI

       ! Did we read the last point?
       IF ( ABS( SSP%z( iz ) - SSP%Depth( Medium + 1 ) ) < 100. * EPSILON( 1.0e0 ) ) THEN
          SSP%NPts( Medium ) = N1
          IF ( Medium == 1 ) SSP%Depth( 1 ) = SSP%z( 1 )
          IF ( SSP%NPts( Medium ) == 1 ) THEN
              WRITE( PRTFile, * ) '#SSP points: ', SSP%NPts( Medium )
              CALL ERROUT( 'ReadSSP', 'The SSP must have at least 2 points in each layer' )
          END IF

          RETURN
       ENDIF

       N1 = N1 + 1
    END DO

    ! Fall through means too many points in the profile
    WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
    CALL ERROUT( 'ReadSSP', 'Number of SSP points exceeds limit' )

  END SUBROUTINE ReadSSP

  !**********************************************************************!

  SUBROUTINE UpdateSSPLoss( freq, freq0 )
    ! Updates the imaginary part of the sound speed based on the frequency
    ! The depth of the SSP point is also used if there is a bio layer
    USE AttenMod
    USE pchipMod

    REAL     (KIND=8), INTENT(IN) :: freq, freq0   ! freq0 is the reference frequency where dB/m was specified
    INTEGER                       :: Medium
    INTEGER                       :: IBCBeg, IBCEnd

    DO Medium = 1, SSP%NMedia
       ILoc = SSP%Loc( Medium )

       DO iSSP = 1, SSP%NPts( Medium )
          iz = SSP%Loc( Medium ) + iSSP
          SSP%cp( iz ) = CRCI( SSP%z( iz ), SSP%alphaR( iz ),  SSP%alphaI( iz ), freq, freq0, &
             SSP%AttenUnit, SSP%beta( Medium), SSP%ft( Medium ) )
          SSP%cs( iz ) = CRCI( SSP%z( iz ), SSP%betaR(  iz ),  SSP%betaI(  iz ), freq, freq0, &
             SSP%AttenUnit, SSP%beta( Medium), SSP%ft( Medium ) )

          SSP%rhoT( iz ) = SSP%rho( iz )

          SSP%cpSpline(  1, iz ) = SSP%cp(  iz )
          SSP%csSpline(  1, iz ) = SSP%cs(  iz )
          SSP%rhoSpline( 1, iz ) = SSP%rho( iz )
       END DO

       ! Compute std polynomial coefs if PCHIP interpolation has been selected
       IF ( SSP%Type == 'P' ) THEN
          CALL PCHIP( SSP%z( ILoc + 1 ), SSP%cp(   ILoc + 1 ), SSP%NPts( Medium ), SSP%cpCoef(  1, ILoc + 1 ), SSP%csWork  )
          CALL PCHIP( SSP%z( ILoc + 1 ), SSP%cs(   ILoc + 1 ), SSP%NPts( Medium ), SSP%csCoef(  1, ILoc + 1 ), SSP%csWork )
          CALL PCHIP( SSP%z( ILoc + 1 ), SSP%rhoT( ILoc + 1 ), SSP%NPts( Medium ), SSP%rhoCoef( 1, ILoc + 1 ), SSP%csWork )
       END IF

       ! Compute spline coefs if spline interpolation has been selected
       IF ( SSP%Type == 'S' ) THEN
          IBCBeg = 0
          IBCEnd = 0
          CALL CSPLINE( SSP%z( ILoc + 1 ), SSP%cpSpline(  1, ILoc + 1 ), SSP%NPts( Medium ), IBCBeg, IBCEnd, SSP%NPts( Medium ) )
          CALL CSPLINE( SSP%z( ILoc + 1 ), SSP%csSpline(  1, ILoc + 1 ), SSP%NPts( Medium ), IBCBeg, IBCEnd, SSP%NPts( Medium ) )
          CALL CSPLINE( SSP%z( ILoc + 1 ), SSP%rhoSpline( 1, ILoc + 1 ), SSP%NPts( Medium ), IBCBeg, IBCEnd, SSP%NPts( Medium ) )
       END IF
    END DO
  END SUBROUTINE UpdateSSPLoss

  
  !**********************************************************************!

  SUBROUTINE UpdateHSLoss( freq, freq0 )

     ! update loss in the halfspaces based on the frequency
     ! depth of HUGE ensures bio loss in the halfspace is not included

     USE AttenMod
     REAL (KIND=8), INTENT(IN) :: freq, freq0   ! freq0 is the reference frequency where dB/m was specified

     IF ( HSTop%BC == 'A' ) THEN
        HSTop%cP  = CRCI( HUGE( 1D0 ), HSTop%alphaR, HSTop%alphaI, freq, freq0, SSP%AttenUnit, HSTop%beta, HSTop%fT )
        HSTop%cS  = CRCI( HUGE( 1D0 ), HSTop%betaR,  HSTop%betaI,  freq, freq0, SSP%AttenUnit, HSTop%beta, HSTop%fT )
     END IF

     IF ( HSBot%BC == 'A' ) THEN
        HSBot%cP  = CRCI( HUGE( 1D0 ), HSBot%alphaR, HSBot%alphaI, freq, freq0, SSP%AttenUnit, HSBot%beta, HSBot%fT )
        HSBot%cS  = CRCI( HUGE( 1D0 ), HSBot%betaR,  HSBot%betaI,  freq, freq0, SSP%AttenUnit, HSBot%beta, HSBot%fT )
     END IF

  END SUBROUTINE UpdateHSLoss

END MODULE sspmod
