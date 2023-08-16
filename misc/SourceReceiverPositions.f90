MODULE SourceReceiverPositions
  ! Reads in source depths, receiver depths, receiver ranges, and receiver bearings

  USE monotonicMod
  USE SortMod
  USE SubTabulate
  USE FatalError

  IMPLICIT NONE
  INTEGER, PARAMETER          :: Number_to_Echo = 10
  INTEGER, PRIVATE            :: IAllocStat     ! used to capture status after allocation
  INTEGER, PRIVATE, PARAMETER :: ENVFile = 5, PRTFile = 6   ! unit 5 is usually (not always) the ENVFile
  INTEGER                     :: Nfreq = 1      ! number of frequencies
  REAL (KIND=8), ALLOCATABLE  :: freqVec( : )   ! frequency vector for braodband runs

  ! Only position depths and their weights are left in single precision below
  ! These are used for mode interpolation
  TYPE Position
     INTEGER              :: NSx = 1, NSy = 1, NSz = 1, NRz = 1, NRr = 1, Ntheta = 1    ! number of x, y, z, r, theta coordinates
     REAL (KIND=8)              :: Delta_r, Delta_theta
     INTEGER,       ALLOCATABLE :: iSz( : ), iRz( : )
     REAL (KIND=8), ALLOCATABLE :: Sx( : ), Sy( : )            ! Source   x, y coordinates
     REAL (KIND=4), ALLOCATABLE :: Sz( : )                     ! Source   z coordinates
     REAL (KIND=8), ALLOCATABLE :: Rr( : )                     ! Receiver r coordinates
     REAL (KIND=4), ALLOCATABLE :: Rz( : )                     ! Receiver z coordinates
     REAL (KIND=4), ALLOCATABLE :: ws( : ), wr( : )            ! weights for source/receiver interpolation in depth
     REAL (KIND=8), ALLOCATABLE :: theta( : )                  ! Receiver bearings
  END TYPE Position

  TYPE ( Position ) :: Pos   ! structure containing source and receiver positions

  INTERFACE ReadVector
     MODULE PROCEDURE ReadVector_sngl, ReadVector_dble
  END INTERFACE ReadVector

CONTAINS

  SUBROUTINE ReadfreqVec( freq0, BroadbandOption )

    ! Optionally reads a vector of source frequencies for a broadband run
    ! If the broadband option is not selected, then the input freq (a scalar) is stored in the frequency vector

    REAL (KIND=8), INTENT( IN ) :: freq0             ! Nominal or carrier frequency
    CHARACTER,     INTENT( IN ) :: BroadbandOption*( 1 )
    INTEGER                     :: ifreq

    ! Broadband run?
    IF ( BroadbandOption == 'B' ) THEN
       READ( ENVFile, * ) Nfreq
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) '   Number of frequencies =', Nfreq
       IF ( Nfreq <= 0 ) CALL ERROUT( 'ReadEnvironment', 'Number of frequencies must be positive'  )
    END IF

    IF ( ALLOCATED( freqVec ) ) DEALLOCATE( freqVec )
    ALLOCATE( freqVec( MAX( 3, Nfreq ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( 'ReadEnvironment', 'Too many frequencies'  )

    IF ( BroadbandOption == 'B' ) THEN
       WRITE( PRTFile, * ) '   Frequencies (Hz)'
       freqVec( 2 ) = -999.9
       freqVec( 3 ) = -999.9
       READ(  ENVFile, * ) freqVec( 1 : Nfreq )
       CALL SubTab( freqVec, Nfreq )

       WRITE( PRTFile, "( 5G14.6 )" ) ( freqVec( ifreq ), ifreq = 1, MIN( Nfreq, Number_to_Echo ) )
       IF ( Nfreq > Number_to_Echo ) WRITE( PRTFile,  "( G14.6 )" ) ' ... ', freqVec( Nfreq )
    ELSE
       freqVec( 1 ) = freq0
    END IF

    RETURN

  END SUBROUTINE ReadfreqVec

  !********************************************************************!

  SUBROUTINE ReadSxSy( ThreeD )

    ! Reads source x-y coordinates

    LOGICAL, INTENT( IN ) :: ThreeD   ! flag indicating whether this is a 3D run

    IF ( ThreeD ) THEN
       CALL ReadVector( Pos%NSx, Pos%Sx, 'Source   x-coordinates, Sx', 'km' )
       CALL ReadVector( Pos%NSy, Pos%Sy, 'Source   y-coordinates, Sy', 'km' )
    ELSE
       ALLOCATE( Pos%Sx( 1 ), Pos%Sy( 1 ) )
       Pos%Sx( 1 ) = 0.
       Pos%Sy( 1 ) = 0.
    END IF

    RETURN
  END SUBROUTINE ReadSxSy

  !********************************************************************!

  SUBROUTINE ReadSzRz( zMin, zMax )

    ! Reads source and receiver z-coordinates (depths)
    ! zMin and zMax are limits for those depths; sources and receivers are shifted to be within those limits

    REAL,    INTENT( IN ) :: zMin, zMax
    !LOGICAL               :: monotonic

    CALL ReadVector( Pos%NSz, Pos%Sz, 'Source   z-coordinates, Sz', 'm' )
    CALL ReadVector( Pos%NRz, Pos%Rz, 'Receiver z-coordinates, Rz', 'm' )

    IF ( ALLOCATED( Pos%ws ) ) DEALLOCATE( Pos%ws, Pos%iSz )
    ALLOCATE( Pos%ws( Pos%NSz ), Pos%iSz( Pos%NSz ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( 'ReadSzRz', 'Too many sources'  )

    IF ( ALLOCATED( Pos%wr ) ) DEALLOCATE( Pos%wr, Pos%iRz )
    ALLOCATE( Pos%wr( Pos%NRz ), Pos%iRz( Pos%NRz ), Stat = IAllocStat  )
    IF ( IAllocStat /= 0 ) CALL ERROUT( 'ReadSzRz', 'Too many receivers'  )

    ! *** Check for Sz/Rz in upper or lower halfspace ***

    IF ( ANY( Pos%Sz( 1 : Pos%NSz ) < zMin ) ) THEN
       WHERE ( Pos%Sz < zMin ) Pos%Sz = zMin
       WRITE( PRTFile, * ) 'Warning in ReadSzRz : Source above or too near the top bdry has been moved down'
    END IF

    IF ( ANY( Pos%Sz( 1 : Pos%NSz ) > zMax ) ) THEN
       WHERE( Pos%Sz > zMax ) Pos%Sz = zMax
       WRITE( PRTFile, * ) 'Warning in ReadSzRz : Source below or too near the bottom bdry has been moved up'
    END IF

    IF ( ANY( Pos%Rz( 1 : Pos%NRz ) < zMin ) ) THEN
       WHERE( Pos%Rz < zMin ) Pos%Rz = zMin
       WRITE( PRTFile, * ) 'Warning in ReadSzRz : Receiver above or too near the top bdry has been moved down'
    END IF

    IF ( ANY( Pos%Rz( 1 : Pos%NRz ) > zMax ) ) THEN
       WHERE( Pos%Rz > zMax ) Pos%Rz = zMax
       WRITE( PRTFile, * ) 'Warning in ReadSzRz : Receiver below or too near the bottom bdry has been moved up'
    END IF

!!$    IF ( .NOT. monotonic( Pos%sz, Pos%NSz ) ) THEN
!!$       CALL ERROUT( 'SzRzRMod', 'Source depths are not monotonically increasing' )
!!$    END IF 
!!$ 
!!$    IF ( .NOT. monotonic( Pos%rz, Pos%NRz ) ) THEN
!!$       CALL ERROUT( 'SzRzRMod', 'Receiver depths are not monotonically increasing' )
!!$    END IF 

    RETURN
  END SUBROUTINE ReadSzRz

  !********************************************************************!

  SUBROUTINE ReadRcvrRanges

    CALL ReadVector( Pos%NRr, Pos%Rr, 'Receiver r-coordinates, Rr', 'km' )

    ! calculate range spacing
    Pos%delta_r = 0.0
    IF ( Pos%NRr /= 1 ) Pos%delta_r = Pos%Rr( Pos%NRr ) - Pos%Rr( Pos%NRr - 1 )

    IF ( .NOT. monotonic( Pos%rr, Pos%NRr ) ) THEN
       CALL ERROUT( 'ReadRcvrRanges', 'Receiver ranges are not monotonically increasing' )
    END IF 
 
    RETURN
  END SUBROUTINE ReadRcvrRanges

  !********************************************************************!

  SUBROUTINE ReadRcvrBearings

    CALL ReadVector( Pos%Ntheta, Pos%theta, 'Receiver bearings, theta', 'degrees' )

    ! full 360-degree sweep? remove duplicate angle
    IF ( Pos%Ntheta > 1 ) THEN
       IF ( ABS( MOD( Pos%theta( Pos%Ntheta ) - Pos%theta( 1 ), 360.0D0 ) ) < 10.0 * TINY( 1.0D0 ) ) &
          Pos%Ntheta = Pos%Ntheta - 1
    END IF

    ! calculate angular spacing
    Pos%Delta_theta = 0.0
    IF ( Pos%Ntheta /= 1 ) Pos%Delta_theta = Pos%theta( Pos%Ntheta ) - Pos%theta( Pos%Ntheta - 1 )

    IF ( .NOT. monotonic( Pos%theta, Pos%Ntheta ) ) THEN
       CALL ERROUT( 'ReadRcvrBearings', 'Receiver bearings are not monotonically increasing' )
    END IF 
 
    RETURN
  END SUBROUTINE ReadRcvrBearings

  !********************************************************************!

  SUBROUTINE ReadVector_sngl( Nx, x, Description, Units )

    ! Read a vector x
    ! Description is something like 'receiver ranges'
    ! Units       is something like 'km'
 
    INTEGER,   INTENT( OUT ) :: Nx
    REAL,      ALLOCATABLE, INTENT( OUT ) :: x( : )
    CHARACTER, INTENT( IN  ) :: Description*( * ), Units*( * )
    INTEGER                  :: ix
   
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) '__________________________________________________________________________'
    WRITE( PRTFile, * )

    READ(  ENVFile, * ) Nx
    WRITE( PRTFile, * ) '   Number of ' // Description // ' = ', Nx

    IF ( Nx <= 0 ) CALL ERROUT( 'ReadVector', 'Number of ' // Description // 'must be positive'  )

    IF ( ALLOCATED( x ) ) DEALLOCATE( x )
    ALLOCATE( x( MAX( 3, Nx ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( 'ReadVector', 'Too many ' // Description )

    WRITE( PRTFile, * ) '   ', Description // ' (' // Units // ')'
    x( 2 ) = -999.9
    x( 3 ) = -999.9
    READ( ENVFile, * ) x( 1 : Nx )

    CALL SubTab( x, Nx )
    CALL Sort(   x, Nx )

    WRITE( PRTFile, "( 5G14.6 )" ) '   ', ( x( ix ), ix = 1, MIN( Nx, Number_to_Echo ) )
    IF ( Nx > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )" ) ' ... ', x( Nx )

    WRITE( PRTFile, * )

    ! Vectors in km should be converted to m for internal use
    IF ( LEN_TRIM( Units ) >= 2 ) THEN
       IF ( Units( 1 : 2 ) == 'km' ) x = 1000.0 * x
    END IF

  END SUBROUTINE ReadVector_sngl

  SUBROUTINE ReadVector_dble( Nx, x, Description, Units )

    ! Read a vector x
    ! Description is something like 'receiver ranges'
    ! Units       is something like 'km'
 
    INTEGER,   INTENT( OUT ) :: Nx
    REAL ( KIND = 8 ), ALLOCATABLE, INTENT( OUT ) :: x( : )
    CHARACTER, INTENT( IN  ) :: Description*( * ), Units*( * )
    INTEGER                  :: ix
   
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) '__________________________________________________________________________'
    WRITE( PRTFile, * )

    READ(  ENVFile, * ) Nx
    WRITE( PRTFile, * ) '   Number of ' // Description // ' = ', Nx

    IF ( Nx <= 0 ) CALL ERROUT( 'ReadVector', 'Number of ' // Description // 'must be positive'  )

    IF ( ALLOCATED( x ) ) DEALLOCATE( x )
    ALLOCATE( x( MAX( 3, Nx ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( 'ReadVector', 'Too many ' // Description )

    WRITE( PRTFile, * ) '   ', Description // ' (' // Units // ')'
    x( 2 ) = -999.9
    x( 3 ) = -999.9
    READ( ENVFile, * ) x( 1 : Nx )

    CALL SubTab( x, Nx )
    CALL Sort(   x, Nx )

    WRITE( PRTFile, "( 5G14.6 )" ) '   ', ( x( ix ), ix = 1, MIN( Nx, Number_to_Echo ) )
    IF ( Nx > Number_to_Echo ) WRITE( PRTFile, "( G14.6 )" ) ' ... ', x( Nx )

    WRITE( PRTFile, * )

    ! Vectors in km should be converted to m for internal use
    IF ( LEN_TRIM( Units ) >= 2 ) THEN
       IF ( Units( 1 : 2 ) == 'km' ) x = 1000.0 * x
    END IF

  END SUBROUTINE ReadVector_dble
  
END MODULE SourceReceiverPositions
