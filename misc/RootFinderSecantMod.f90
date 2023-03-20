MODULE RootFinderSecantMod

  ! Secant method

  ! Input:
  !   x2 is an initial guess
  !   Tolerance is an error bound for the returned root
  !   MaxIteration is the maximum allowable number of iterations

  ! Output:
  !   x2 is the estimated root
  !   Iteration is the number of iterations used
  !   ErrorMessage is empty unless there was a problem

  IMPLICIT NONE
  INTEGER, PRIVATE :: iPower0, iPower1

  INTERFACE RootFinderSecant
     MODULE PROCEDURE ZsecantX, ZsecantCX
  END INTERFACE RootFinderSecant

CONTAINS

  SUBROUTINE ZSecantX( x2, Tolerance, Iteration, MaxIteration, ErrorMessage, Funct ) 

    INTEGER,            INTENT( IN    ) :: MaxIteration
    INTEGER,            INTENT( OUT   ) :: Iteration
    REAL ( KIND=8 ),    INTENT( IN    ) :: Tolerance
    REAL ( KIND=8 ),    INTENT( INOUT ) :: x2
    CHARACTER (LEN=80), INTENT( OUT   ) :: ErrorMessage
    EXTERNAL :: Funct   ! Name of the external subroutine that supplies the function whose roots are sought
    REAL ( KIND=8 )                     :: x0, x1, shift, f0, f1, cNum, cDen

    ErrorMessage = ' '
    IF ( Tolerance <= 0.0D0 ) THEN 
       ErrorMessage = 'Non-positive tolerance specified' 
       STOP 
    END IF

    x1 = x2 + 10.0D0 * Tolerance
    CALL FUNCT( x1, f1, iPower1 )
    !WRITE( *, * )
    !WRITE( *, FMT="( 2G24.16, I5 )" ) SQRT( x1 ), f1, iPower1

    DO Iteration = 1, MaxIteration 
       x0      = x1
       f0      = f1 
       iPower0 = iPower1 
       x1      = x2

       CALL FUNCT( x1, f1, iPower1 )

!!$     IF ( F1 == 0.0D0 ) THEN 
!!$        shift = 0.0D0
!!$     ELSE 
!!$        shift = ( x1 - x0 ) / ( 1.0D0 - F0 / F1 * 10.0D0 ** ( iPower0 - iPower1 ) )
!!$     ENDIF

       ! ugly stuff to block overflows by forcing shift to be bounded
       cNum = f1 * ( x1 - x0 )
       cDen = f1 - f0 * 10.0D0 ** ( iPower0 - iPower1 )

       IF ( ABS( cNum ) >= ABS( cDen * x1 ) ) THEN 
          shift = 0.1D0 * Tolerance
       ELSE
          shift = cNum / cDen 
       ENDIF

       x2 = x1 - shift
       !WRITE( *, FMT="( 3G24.16, I5, G24.16 )" ) SQRT( x1 ), x1, F1, iPower1, x2

       ! The following convergence test is very restrictive
       ! Found it was necessary to check all 3 points (x0, x1, x2) for certain problems ...
       IF ( ABS( x2 - x1 ) + ABS( x2 - x0 ) < Tolerance ) THEN
          ! WRITE( *, * ) 'converged', x0, x1, x2, Tolerance
          RETURN
       END IF
    END DO

    ErrorMessage = 'Failure to converge in RootFinderSecant'

  END SUBROUTINE ZSecantX

  ! _____________________________________

  SUBROUTINE ZSecantCX( x2, Tolerance, Iteration, MaxIteration, ErrorMessage, Funct )

    INTEGER,            INTENT( IN    ) :: MaxIteration
    INTEGER,            INTENT( OUT   ) :: Iteration
    REAL    ( KIND=8 ), INTENT( IN    ) :: Tolerance
    COMPLEX ( KIND=8 ), INTENT( INOUT ) :: x2
    CHARACTER (LEN=80), INTENT( OUT   ) :: ErrorMessage 
    EXTERNAL :: Funct   ! Name of the external subroutine that supplies the function whose roots are sought
    COMPLEX ( KIND=8 )                  :: x0, x1, shift, f0, f1, cNum, cDen

    ErrorMessage = ' '
    IF ( Tolerance <= 0.0D0 ) THEN 
       ErrorMessage = 'Non-positive tolerance specified' 
       STOP 
    END IF

    x1 = x2 + 100.0D0 * Tolerance
    CALL FUNCT( x1, f1, iPower1 )
    !WRITE( *, * )
    !WRITE( *, "( 4G17.9, I5 )" ) SQRT( x1 ), f1, iPower1

    DO Iteration = 1, MaxIteration 
       x0      = x1
       f0      = f1 
       iPower0 = iPower1 
       x1      = x2

       CALL FUNCT( x1, f1, iPower1 )

       ! ugly stuff to block overflows by forcing shift to be bounded
       cNum = f1 * ( x1 - x0 )
       cDen = f1 - f0 * 10.0D0 ** ( iPower0 - iPower1 )

       IF ( ABS( cNum ) >= ABS( cDen * x1 ) ) THEN 
          shift = 0.1D0 * Tolerance
       ELSE
          shift = cNum / cDen 
       ENDIF

       x2 = x1 - shift
       !WRITE( *, "( 6G24.16, I5 )" ) SQRT( x1 ), x1, f1, iPower1

       ! The following convergence test is very restrictive
       ! Found it was necessary to check all 3 points (x0, x1, x2) for certain problems ...
       IF ( ABS( x2 - x1 ) + ABS( x2 - x0 ) < Tolerance ) THEN
          !WRITE( *, * ) 'converged', sqrt( x2 ) ! x0, x1, x2, Tolerance
          RETURN
       END IF
    END DO

    ErrorMessage = 'Failure to converge in RootFinderSecant'

  END SUBROUTINE ZSecantCX

END MODULE RootFinderSecantMod

