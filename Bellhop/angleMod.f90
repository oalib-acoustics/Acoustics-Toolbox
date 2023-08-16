MODULE anglemod

  USE MathConstants
  USE SubTabulate
  USE SourceReceiverPositions
  USE SortMod
  USE FatalError

  IMPLICIT NONE
  SAVE

  INTEGER          :: ialpha, ibeta
  INTEGER, PRIVATE :: AllocateStatus
  INTEGER,       PRIVATE, PARAMETER :: ENVFile = 5, PRTFile = 6
  REAL (KIND=8), PRIVATE, PARAMETER :: c0 = 1500.0

  TYPE AnglesStructure
     INTEGER       :: Nalpha = 0, Nbeta = 1, iSingle_alpha = 0, iSingle_beta = 0
     REAL (KIND=8) :: Dalpha, Dbeta
     REAL (KIND=8), ALLOCATABLE:: alpha( : ), beta( : )
  END TYPE AnglesStructure

  Type( AnglesStructure ) :: Angles

CONTAINS
  SUBROUTINE ReadRayElevationAngles( freq, Depth, TopOpt, RunType )

    REAL      (KIND=8), INTENT( IN  ) :: freq, Depth
    CHARACTER (LEN= 6), INTENT( IN  ) :: TopOpt, RunType
    REAL      (KIND=8)                :: d_theta_recommended

    IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
       READ( ENVFile, * ) Angles%Nalpha, Angles%iSingle_alpha ! option to trace a single beam
    ELSE
       READ( ENVFile, * ) Angles%Nalpha
    END IF

    IF ( Angles%Nalpha == 0 ) THEN   ! automatically estimate Nalpha to use
       IF ( RunType( 1 : 1 ) == 'R' ) THEN
          Angles%Nalpha = 50         ! For a ray trace plot, we don't want too many rays ...
       ELSE
          ! you're letting ME choose? OK: ideas based on an isospeed ocean
          ! limit based on phase of adjacent beams at maximum range
          Angles%Nalpha = MAX( INT( 0.3 * Pos%Rr( Pos%NRr ) * freq / c0 ), 300 )

          ! limit based on having beams that are thin with respect to the water depth
          ! assumes also a full 360 degree angular spread of rays
          ! Should check which Depth is used here, in case where there is a variable bathymetry
          d_theta_recommended = ATAN( Depth / ( 10.0 * Pos%Rr( Pos%NRr ) ) )
          Angles%Nalpha = MAX( INT( pi / d_theta_recommended ), Angles%Nalpha )
       END IF
    END IF

    ALLOCATE( Angles%alpha( MAX( 3, Angles%Nalpha ) ), STAT = AllocateStatus )
    IF ( AllocateStatus /= 0 ) CALL ERROUT( 'ReadRayElevationAngles', 'Insufficient memory to store beam angles'  )

    IF ( Angles%Nalpha > 2 ) Angles%alpha( 3 ) = -999.9
    READ( ENVFile, * ) Angles%alpha

    CALL SubTab( Angles%alpha, Angles%Nalpha )
    CALL Sort(   Angles%alpha, Angles%Nalpha )

    ! full 360-degree sweep? remove duplicate beam
    IF ( Angles%Nalpha > 1 .AND. ABS( MOD( Angles%alpha( Angles%Nalpha ) - Angles%alpha( 1 ), 360.0D0 ) ) < &
         10.0 * SPACING( 360.0D0 ) ) &
       Angles%Nalpha = Angles%Nalpha - 1

    WRITE( PRTFile, * ) '__________________________________________________________________________'
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) '   Number of beams in elevation   = ', Angles%Nalpha
    IF ( Angles%iSingle_alpha > 0 ) WRITE( PRTFile, * ) 'Trace only beam number ', Angles%iSingle_alpha
    WRITE( PRTFile, * ) '   Beam take-off angles (degrees)'

    IF ( Angles%Nalpha >= 1 ) WRITE( PRTFile, "( 5G14.6 )" ) Angles%alpha( 1 : MIN( Angles%Nalpha, Number_to_Echo ) )
    IF ( Angles%Nalpha > Number_to_Echo ) WRITE( PRTFile,  "( G14.6 )" ) ' ... ', Angles%alpha( Angles%Nalpha )

    IF ( Angles%Nalpha > 1 .AND. Angles%alpha( Angles%Nalpha ) == Angles%alpha( 1 ) ) &
         CALL ERROUT( 'ReadRayElevationAngles', 'First and last beam take-off angle are identical' )

    IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
       IF ( Angles%iSingle_alpha < 1 .OR. Angles%iSingle_alpha > Angles%Nalpha ) &
            CALL ERROUT( 'ReadRayElevationAngles', 'Selected beam, iSingle_alpha not in [ 1, Angles%Nalpha ]' )
    END IF

  END SUBROUTINE ReadRayElevationAngles

  !**********************************************************************!

  SUBROUTINE ReadRayBearingAngles( freq, TopOpt, RunType )

    REAL      (KIND=8), INTENT( IN ) :: freq
    CHARACTER (LEN= 6), INTENT( IN ) :: TopOpt, RunType

    IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
       READ( ENVFile, * ) Angles%Nbeta, Angles%iSingle_beta ! option to trace a single beam
    ELSE
       READ( ENVFile, * ) Angles%Nbeta
    END IF

    IF ( Angles%Nbeta == 0 ) THEN   ! automatically estimate Nbeta to use
       IF ( RunType( 1 : 1 ) == 'R' ) THEN
          Angles%Nbeta = 50         ! For a ray trace plot, we don't want too many rays ...
       ELSE
          Angles%Nbeta = MAX( INT( 0.1 * Pos%rr( Pos%NRr ) * freq / c0 ), 300 )
       END IF
    END IF

    ALLOCATE( Angles%beta( MAX( 3, Angles%Nbeta ) ), STAT = AllocateStatus )
    IF ( AllocateStatus /= 0 ) CALL ERROUT( 'ReadRayBearingAngles', 'Insufficient memory to store beam angles'  )

    IF ( Angles%Nbeta > 2 ) Angles%beta( 3 ) = -999.9
    READ( ENVFile, * ) Angles%beta

    CALL SubTab( Angles%beta, Angles%Nbeta )
    CALL Sort(   Angles%beta, Angles%Nbeta )

    ! full 360-degree sweep? remove duplicate beam
    IF ( Angles%Nbeta > 1 .AND. ABS( MOD( Angles%beta( Angles%Nbeta ) - Angles%beta( 1 ), 360.0D0 ) ) < 10.0 * SPACING( 1.0D0 ) ) &
         Angles%Nbeta = Angles%Nbeta - 1

    ! Nx2D CASE: beams must lie on rcvr radials--- replace beta with theta
    IF ( RunType( 6 : 6 ) == '2' .AND. RunType( 1 : 1 ) /= 'R' ) THEN
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Replacing beam take-off angles, beta, with receiver bearing lines, theta'
       DEALLOCATE( Angles%beta )

       Angles%Nbeta = Pos%Ntheta
       ALLOCATE( Angles%beta( MAX( 3, Angles%Nbeta ) ), STAT = AllocateStatus )
       IF ( AllocateStatus /= 0 ) CALL ERROUT( 'ReadRayBearingAngles', 'Insufficient memory to store beam angles'  )
       Angles%beta( 1 : Angles%Nbeta ) = Pos%theta( 1 : Pos%Ntheta )   ! Nbeta should = Ntheta
    END IF

    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) '   Number of beams in bearing   = ', Angles%Nbeta
    IF ( Angles%iSingle_beta > 0 ) WRITE( PRTFile, * ) 'Trace only beam number ', Angles%iSingle_beta
    WRITE( PRTFile, * ) '   Beam take-off angles (degrees)'

    IF ( Angles%Nbeta >= 1 ) WRITE( PRTFile, "( 5G14.6 )" ) Angles%beta( 1 : MIN( Angles%Nbeta, Number_to_Echo ) )
    IF ( Angles%Nbeta > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )" ) ' ... ', Angles%beta( Angles%Nbeta )

    IF ( Angles%Nbeta > 1 .AND. Angles%beta( Angles%Nbeta ) == Angles%beta( 1 ) ) &
         CALL ERROUT( 'ReadRayBearingAngles', 'First and last beam take-off angle are identical' )

    IF ( TopOpt( 6 : 6 ) == 'I' ) THEN
       IF ( Angles%iSingle_beta < 1 .OR. Angles%iSingle_beta > Angles%Nbeta ) &
            CALL ERROUT( 'ReadRayBearingAngles', 'Selected beam, iSingle_beta not in [ 1, Angles%Nbeta ]' )
    END IF
    Angles%beta  = DegRad * Angles%beta   ! convert to radians

    Angles%Dbeta = 0.0
    IF ( Angles%Nbeta /= 1 ) Angles%Dbeta = ( Angles%beta( Angles%NBeta ) - Angles%beta( 1 ) ) / ( Angles%Nbeta - 1 )

  END SUBROUTINE ReadRayBearingAngles

END MODULE anglemod
