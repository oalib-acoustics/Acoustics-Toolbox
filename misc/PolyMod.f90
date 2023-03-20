MODULE PolyMod

  ! Polynomial approximant of order N at X0
  ! mbp 7/2015 incorporating subroutines from decades past

  IMPLICIT NONE
  INTEGER, PRIVATE :: i, j

  INTERFACE Poly
     MODULE PROCEDURE PolyR, PolyC, PolyZ
  END INTERFACE Poly

CONTAINS

  FUNCTION PolyR( X0, X, F, N )

    INTEGER, INTENT( IN ) :: N   ! order of the polynomial
    REAL,    INTENT( IN ) :: x0, x( N ), f( N )  ! x, y values of the polynomial
    REAL                  :: ft( N ), h( N ), PolyR

    ! Initialize arrays

    h  = x - x0
    ft = f

    ! Recursion for solution
    IF ( N >= 2) THEN
       DO i = 1, N - 1
          DO j = 1, N - i
             ft( j ) = ( h( j + i ) * ft( i ) - h( i ) * ft( i + 1 ) ) / &
                  &            ( h( j + i ) - h( j ) )
          END DO
       END DO
    ENDIF

    PolyR = ft( 1 )

  END FUNCTION PolyR

  !__________________________________________________________________________

  COMPLEX FUNCTION PolyC( x0, x, f, N )

    INTEGER, INTENT( IN ) :: N   ! order of the polynomial
    COMPLEX, INTENT( IN ) :: x0, x( N ), f( N )  ! x, y values of the polynomial
    COMPLEX               :: ft( N ), h( N )

    ! Initialize arrays
    h  = x - x0
    ft = f

    ! Recursion for solution
    IF ( N >= 2) THEN
       DO i = 1, N-1
          DO j = 1, N-I
             !         ft( J ) = ( h( J+I ) * ft( J ) - h( J ) * ft( J+1 ) ) / &
             !                                   ( h( J+I ) - h( J ) )
             ft( J ) = ft( J ) + h( J ) * ( ft( J )  - ft( J+1 ) ) / &
                  &                       ( h( J+I ) - h(  J   ) )
          END DO
       END DO
    ENDIF

    PolyC = ft( 1 )

  END FUNCTION PolyC

  !__________________________________________________________________________

  FUNCTION PolyZ( x0, x, F, N )

    INTEGER,            INTENT( IN ) :: N  ! order of the polynomial
    COMPLEX ( KIND=8 ), INTENT( IN ) :: x0, x( N ), f( N ) ! x, y values of the polynomial
    COMPLEX ( KIND=8 )    :: PolyZ
    COMPLEX ( KIND=8 )    :: fT( N ), h( N )

    ! Initialize arrays
    h  = x - x0
    fT = f

    ! Recursion for solution
    IF ( N >= 2 ) THEN
       DO i = 1, N - 1
          DO j = 1, N - i
             fT( j ) = ( h( j + i ) * fT( j ) - h( j ) * fT( j + 1 ) ) / ( h( j + i ) - h( j ) )
          END DO
       END DO
    ENDIF
    PolyZ = fT( 1 )

  END FUNCTION PolyZ

END MODULE PolyMod

