COMPLEX FUNCTION Pade( x0, x, f, N, ErrorMessage )

  !     Computes the [MU, NU] Pade approximant at x0, where
  !        MU = [M/2], NU = M - MU and M = N - 1
  !
  !     FORMULAS TAKEN FROM BURLISCH AND STOER
  !        NUMERISCHE MATHEMATIK 8:1-13(1966)
  !
  !     N IS THE NUMBER OF POINTS USED.
  !     x CONTAINS THE INDEPENDENT VARIABLE
  !     f CONTAINS THE FUNCTIONAL VALUES AT EACH x
  !     x0 IS THE POINT AT WHICH THE APPROXIMATION IS TO BE CALCULATED
  !
  !     Note much wasted space in this implementation ...
  !
  !     Michael Porter, 11 September 1985

  IMPLICIT NONE
  INTEGER, INTENT( IN ) :: N   ! order of the Pade approximant
  COMPLEX, INTENT( IN ) :: x0, x( N ), f( N )  ! x, y values of the polynomila
  CHARACTER (LEN=80), INTENT( OUT ) :: ErrorMessage
  !INTEGER              :: II( 1 )
  INTEGER               :: K, M, Nearest
  COMPLEX               :: RDT, WDen, CNum, CDen, D( N, N ), C( N, N ), H( N )

  ErrorMessage = '      '
  H = x - x0   ! distance between x0 and every point in x

  ! Quick return if x lies "on" a data point
  Nearest = 1
  !II      = MINLOC( ABS( H ) )
  !Nearest = II( 1 )
  !IF ( ABS( H( Nearest ) ) < 100 * SPACING( REAL( x0 ) ) ) THEN
  !   Pade = f( Nearest )
  !   RETURN
  !ENDIF

  ! Recursion for solution
  DO M = 1, N
     D( 1, M ) = f( M ) + EPSILON( REAL( x0 ) )  ! perturbation to avoid 0/0 condition
     C( 1, M ) = f( M )

     IF ( M >= 2 ) THEN

        DO K = 1, M - 1
           RDT  = H( M - K ) / H( M ) * D( K, M - K )
           CNum = C( K, M - K + 1 ) - D( K, M - K )
           CDen = RDT - C( K, M - K + 1 )

           IF ( ABS( CNum ) > 1.0E10 * ABS( CDen ) ) THEN
              ErrorMessage = 'ERROR: Nearly-singular matrix in Pade'
              Pade = f( Nearest )
              RETURN
           ENDIF

           WDen = CNum / CDen
           D( K + 1, M - K ) = C( K, M - K + 1 ) * WDen
           C( K + 1, M - K ) = RDT * WDen
        END DO
     ENDIF
  END DO

  ! Sum diagonal elements to get interpolate
  Pade = 0.0
  DO K = 1, N
     Pade = Pade + D( K, N - K + 1 )
  END DO

END FUNCTION Pade
