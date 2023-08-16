MODULE InverseIterationMod

  !     Michael B. Porter 7/1/85                                          

  !     This subroutine is a modified version of TINVIT in EISPACK
  !     specialized to the case where we only want a single vector.

  !     WARNING: Very probably this routine will fail for
  !     certain difficult cases for which the original tinvit would       
  !     not.  For instance, problems with nearly degenerate eigenvectors  
  !     or separable problems.

  !     Refer to the ALGOL procedure TRISTURM by Peters and Wilkinson.    
  !     Handbook for Auto. Comp., VOL.II-Linear Algebra, 418-439(1971).   

  !     This subroutine finds the eigenvector of a singular, tridiagonal  
  !     symmetric matrix corresponding to the zero eigenvalue, using inverse iteration.

  !     On input
  !        N is the order of the matrix.                                  
  !        D contains the diagonal elements of the input matrix.          
  !        E contains the subdiagonal elements of the input matrix        
  !          in its last N-1 positions.  E( 1 ), E( N + 1 ) are arbitrary.      

  !     On output                                                         
  !        All input arrays are unaltered.                                
  !        EigenVector contains the eigenvector.                                  
  !          Any vector which fails to converge is set to zero.           
  !        iError is set to                                                 
  !           0       For normal return,                                  
  !          -1       If the eigenvector fails to converge in MaxIteration iterations

  !        RV1, RV2, RV3, RV4, and EigenVector are temporary storage arrays.      

  IMPLICIT NONE
  INTEGER, PARAMETER :: MaxIteration = 3  ! Max number of inverse iterations
  INTEGER            :: I, Iteration
  REAL    (KIND=8)   :: EPS3, EPS4, NORM

  INTERFACE inverseiteration
     MODULE PROCEDURE inverseiterationD, inverseiterationZ
  END INTERFACE inverseiteration

CONTAINS

  SUBROUTINE InverseIterationD( N, D, E, iError, EigenVector ) 

    ! double precision, real version
    INTEGER,       INTENT(  IN ) :: N
    INTEGER,       INTENT( OUT ) :: iError
    REAL (KIND=8), INTENT(  IN ) :: D( N ), E( N + 1 )
    REAL (KIND=8), INTENT( OUT ) :: EigenVector( N )
    REAL (KIND=8)                :: U, UK, V, XU
    REAL (KIND=8)                :: RV1( N ), RV2( N ), RV3( N ), RV4( N )

    iError = 0 

    ! Compute (infinity) norm of matrix
    ! (this could be pre-calculated for addl speed ...)
    NORM = SUM( ABS( D ) ) + SUM( ABS( E( 2 : N ) ) )

    ! added a factor of 100 below because some Scholte modes weren't getting amplified enough
    ! mbp: 8/2010

    EPS3 = 100.0D0 * EPSILON( NORM ) * NORM  ! small number that replaces zero pivots
    UK   = N 
    EPS4 = UK * EPS3               ! small number used in scaling down the iterates
    UK   = EPS4 / DSQRT( UK ) 

    ! *** Elimination with interchanges ***                             

    XU = 1.0D0
    U  = D( 1 )
    V  = E( 2 )

    DO I = 2, N       
       IF ( ABS( E( I ) ) >= ABS( U ) ) THEN
          XU = U / E( I )
          RV4( I     ) = XU 
          RV1( I - 1 ) = E( I ) 
          RV2( I - 1 ) = D( I ) 
          RV3( I - 1 ) = E( I + 1 ) 
          U = V - XU * RV2( I - 1 ) 
          V =   - XU * RV3( I - 1 ) 
       ELSE 
          XU = E( I ) / U 
          RV4( I     ) = XU
          RV1( I - 1 ) = U
          RV2( I - 1 ) = V
          RV3( I - 1 ) = 0.0D0
          U = D( I ) - XU * V
          V = E( I + 1 )
       END IF
    END DO
    IF ( U == 0.0D0 ) U = EPS3 

    RV3( N - 1 ) = 0.0D0
    RV1( N     ) = U
    RV2( N     ) = 0.0D0
    RV3( N     ) = 0.0D0

    ! *** Main loop of inverse iteration ***

    EigenVector = UK

    DO Iteration = 1, MaxIteration 

       ! ***  Back substitution                                         
       DO I = N, 1, -1
          EigenVector( I ) = ( EigenVector( I ) - U * RV2( I ) - V * RV3( I ) ) / RV1( I )
          V = U 
          U = EigenVector( I )
       END DO

       ! *** Compute norm of vector and test for doneness
       NORM = SUM( ABS( EigenVector ) )
       IF ( NORM >= 1.0D0 ) RETURN   ! convergence so return

       XU  = EPS4 / NORM   ! Scale the vector down
       EigenVector = EigenVector * XU

       ! *** Forward elimination
       DO I = 2, N 
          U = EigenVector( I )

          IF ( RV1( I - 1 ) == E( I ) ) THEN   ! rows switched during triangularization
             U                    = EigenVector( I - 1 ) 
             EigenVector( I - 1 ) = EigenVector( I     ) 
          END IF
          EigenVector( I ) = U - RV4( I ) * EigenVector( I - 1 )
       END DO

    END DO   ! next iteration

    ! *** Fall through the loop means failure to converge               
    iError = -1
  END SUBROUTINE InverseIterationD

  ! _____________________________________

  SUBROUTINE InverseIterationZ( N, D, E, iError, EigenVector ) 

    ! double precision, complex version
    
    INTEGER,          INTENT(  IN ) :: N
    INTEGER,          INTENT( OUT ) :: iError
    COMPLEX (KIND=8), INTENT(  IN ) :: D( N ), E( N + 1 )
    COMPLEX (KIND=8), INTENT( OUT ) :: EigenVector( N )
    REAL    (KIND=8)                :: UK
    COMPLEX (KIND=8)                :: U, V, XU
    COMPLEX (KIND=8)                :: RV1( N ), RV2( N ), RV3( N ), RV4( N )

    iError = 0 

    ! Compute (infinity) norm of matrix
    ! (this could be pre-calculated for addl speed ...)
    NORM = SUM( ABS( REAL( D          ) ) ) + SUM( ABS( AIMAG( D          ) ) )  + &
      &    SUM( ABS( REAL( E( 2 : N ) ) ) ) + SUM( ABS( AIMAG( E( 2 : N ) ) ) )

    ! added a factor of 100 below because some Scholte modes weren't getting amplified enough
    ! mbp: 8/2010

    EPS3 = 100.0D0 * EPSILON( NORM ) * NORM  ! small number that replaces zero pivots
    UK   = N 
    EPS4 = UK * EPS3               ! small number used in scaling down the iterates
    UK   = EPS4 / DSQRT( UK ) 

    ! *** Elimination with interchanges ***                             

    XU = 1.0D0
    U  = D( 1 )
    V  = E( 2 )

    DO I = 2, N
       ! can replace abs with test on real, imag parts for speed ...
       IF ( ABS( E( I ) ) >= ABS( U ) ) THEN
          XU = U / E( I ) 
          RV4( I     ) = XU 
          RV1( I - 1 ) = E( I ) 
          RV2( I - 1 ) = D( I ) 
          RV3( I - 1 ) = E( I + 1 ) 
          U = V - XU * RV2( I - 1 ) 
          V =   - XU * RV3( I - 1 ) 
       ELSE 
          XU = E( I ) / U 
          RV4( I     ) = XU
          RV1( I - 1 ) = U
          RV2( I - 1 ) = V
          RV3( I - 1 ) = 0.0D0
          U = D( I ) - XU * V
          V = E( I + 1 )
       END IF
    END DO

    IF ( U == 0.0D0 ) U = EPS3 

    RV3( N - 1 ) = 0.0D0
    RV1( N     ) = U
    RV2( N     ) = 0.0D0
    RV3( N     ) = 0.0D0

    ! *** Main loop of inverse iteration ***

    EigenVector = UK

    DO Iteration = 1, MaxIteration 

       ! ***  Back substitution                                         
       DO I = N, 1, -1
          EigenVector( I ) = ( EigenVector( I ) - U * RV2( I ) - V * RV3( I ) ) / RV1( I )
          V = U 
          U = EigenVector( I )
       END DO

       ! *** Compute norm of vector and test for doneness
       NORM = SUM( ABS( REAL( EigenVector ) ) ) + SUM( ABS( AIMAG( EigenVector ) ) )
       IF ( NORM >= 1.0D0 ) RETURN   ! convergence so return

       XU  = EPS4 / NORM   ! Scale the vector down
       EigenVector = EigenVector * XU

       ! *** Forward elimination
       DO I = 2, N 
          U = EigenVector( I )

          IF ( RV1( I - 1 ) == E( I ) ) THEN   ! rows switched during triangularization
             U            = EigenVector( I - 1 ) 
             EigenVector( I - 1 ) = EigenVector( I     ) 
          END IF
          EigenVector( I ) = U - RV4( I ) * EigenVector( I - 1 )
       END DO

    END DO   ! next iteration

    ! *** Fall through the loop means failure to converge               

    iError = -1
  END SUBROUTINE InverseIterationZ

END MODULE InverseIterationMod

