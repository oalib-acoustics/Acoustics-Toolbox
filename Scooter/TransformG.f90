SUBROUTINE HankelTransform( Nt, WKMin, WKDel, Stabil, RMin, NR, RDel, X, Y, Z, OPTION )

  !     Hankel Transform by series approximation - Currently first term only.
  !     Formulas from Abramowitz and Stegun, Eqs. 9.2.5-9.2.10
  !     X      = input k space spectra, Nt points.
  !     Z      = output transmission loss, Nr points.
  !     Nt     = number of points in transform.
  !     WKMin  = minimum wavenumber of input series.
  !     WKDel  = spacing of input series.
  !     Stabil = amount integration has been moved off of the real axis.
  !     Nr     = number of output transmission loss points.(WITH Nr < 2*Nt)
  !     RDel   = spacing of transmission loss points.
  !     RMin   = Starting point for TL

  !     Adapted from Steve Wales' original program

  IMPLICIT NONE
  REAL      pi, sqrt2pi
  COMPLEX   i, ipiD4
  PARAMETER ( i = (0.0, 1.0), pi = 3.141592653589793 )
  INTEGER   ir, ik, Nt, Nr
  REAL      R( Nr ), WKVEC( Nt ), RMin, WKMin, Rdel, WKdel, Stabil
  COMPLEX   X( Nt ), Y( Nt ), Z( Nr ), CK, iStab, iRMin
  CHARACTER (LEN=3) :: OPTION

  ! Initialization

  ipiD4   = i * pi / 4.0
  iStab   = i * Stabil
  iRMin   = i * RMin
  Sqrt2pi = 1.0 / SQRT ( 2. * pi )
  CK      = CMPLX ( +Stabil, -WKMin )

  ! Following should be pulled out of the subroutine for pre-computation
  ! It's here as a safety against passing an irregularly spaced vector

  R(     1 : Nr ) = RMin + [ ( IR, IR = 0, Nr - 1 ) ] * RDel
  WKVEC( 1 : Nt ) =        [ ( IK, IK = 0, Nt - 1 ) ] * WKDel

  IF ( WKMin > 0.5 * WKDel )  THEN
     X( 1 ) = WkDel * Sqrt2pi *                    SQRT( WKMin + iStab ) * X( 1 )
  ELSE
     X( 1 ) = 0.5 * ( iStab + WKDel ) * Sqrt2pi *  SQRT(         iStab ) * X( 1 )
  END IF
  X( 2 : Nt ) = WkDel * Sqrt2pi * SQRT(   WKMin + WKVEC( 2 : Nt )  + iStab ) * X( 2 : Nt )

  ! First term, lower-half of cosine transform
  IF ( OPTION(2 : 2) == 'P' .OR. OPTION(2 : 2) == 'B' ) THEN   ! (P)ositive or (B)oth positive and negative spectrum
     Y = X * EXP( -iRMin * WKVEC + ipiD4 )
     CALL CFFT( Y, Nt, 1 )                 ! EXP(-ikX) transform
     Z = Y( 1 : Nr ) * EXP ( CK * R )
  ENDIF

  !  Second term, second-half of cosine transform
  IF ( OPTION(2 : 2) == 'N' .OR. OPTION(2 : 2) == 'B' ) THEN   ! (N)egative or (B)oth positive and negative spectrum
     Y = X * EXP( +iRMin * WKVEC - ipiD4 )
     CALL CFFT( Y, Nt, -1 )                ! EXP(+ikX) transform
  ENDIF

  IF      ( OPTION(2 : 2) == 'B' ) THEN   ! both sides of spectrum
     Z = Z + Y( 1 : Nr ) * EXP ( -CK * R )
  ELSE IF ( OPTION(2 : 2) == 'N' ) THEN   ! negative spectrum only
     Z =     Y( 1 : Nr ) * EXP ( -CK * R )
  ENDIF

  IF ( OPTION(1 : 1) == 'R' ) Z = Z / SQRT( R )   ! cylindrical spreading

END SUBROUTINE HankelTransform

!**********************************************************************

SUBROUTINE FourierTransform( Nt, WKMin, WKDel, Stabil, RMin, Nr, RDel, X, Y, Z )

  !     Fourier Transform
  !     X      = input k space spectra, Nt points.
  !     Z      = output transmission loss, Nr points.
  !     Nt     = number of points in transform.
  !     WKMin  = minimum wavenumber of input series.
  !     WKDel  = spacing of input series.
  !     Stabil = amount integration has been moved off of the real axis.
  !     Nr     = number of output transmission loss points.(WITH Nr <= 2*Nt)
  !     RDel   = spacing of transmission loss points.
  !     RMin   = Starting point for TL

  IMPLICIT NONE
  REAL      pi, sqrt2pi
  COMPLEX   i
  PARAMETER ( i = (0.0, 1.0), pi = 3.141592653589793 )
  INTEGER   ir, ik, Nt, Nr
  REAL      R( Nr ), WKVEC( Nt ), RMin, WKMin, Rdel, WKdel, Stabil
  COMPLEX   X( Nt ), Y( Nt ), Z( Nr ), CK, iStab, iRMin

  ! Initialization
  iStab   = i * Stabil
  iRMin   = i * RMin
  Sqrt2pi = 1. / SQRT ( 2. * pi )

  ! Following should be pulled out of the subroutine for pre-computation
  ! It's here as a safety against passing an irregularly spaced vector

  R(     1 : Nr ) = RMin + [ ( IR, IR = 0, Nr - 1 ) ] * RDel
  WKVEC( 1 : Nt ) =        [ ( IK, IK = 0, Nt - 1 ) ] * WKDel

  IF ( WKMin > 0.5 * WKDel )  THEN
     X( 1 ) = WKDel * Sqrt2pi * X( 1 )
  ELSE
     X( 1 ) = ( iStab + WKDel ) / 2.0 * Sqrt2pi * X( 1 )
  END IF
  X( 2 : Nt ) = WKDel * Sqrt2pi * X( 2 : Nt )

  ! First term, lower-half of cosine transform
  Y = X * EXP( -iRMin * WKVEC )
  CALL CFFT ( Y, Nt, +1 )                   ! EXP(-ikX) transform
  CK = CMPLX( +Stabil, -WKMin )
  Z =             Y( 1 : Nr ) * EXP ( CK * R )

  ! Second term, second-half of cosine transform
  Y = X * EXP( +iRMin * WKVEC )
  CALL CFFT( Y, Nt, -1 )                    ! EXP(+ikX) transform
  Z = Z( 1 : Nr ) + Y( 1 : Nr ) * EXP ( -CK * R )

END SUBROUTINE FourierTransform
