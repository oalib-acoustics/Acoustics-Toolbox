      SUBROUTINE CSPLINE (TAU, C, N, IBCBEG, IBCEND, NDIM)

!  TAKEN FROM "A PRACTICAL GUIDE TO SPLINES", BY CARL DE BOOR. 1978.
!  SPRINGER-VERLAG.  THE INPUT PARAMETER "NDIM" HAS BEEN ADDED TO
!  ALLOW FOR MULTIPLE CALLS WITH DIFFERENT VALUES OF N. - DENNIS DUNDORE

!  SUBSTANTIAL MODIFICATIONS MADE BY STEVE WALES, APRIL 1983, 
!  PRINCIPALLY TO HANDLE COMPLEX NUMBERS (C) & UPDATE THE FORTRAN.

! *****************************  I N P U T  ****************************

!  N = NUMBER OF DATA POINTS.  ASSUMED TO BE .GE. 2.

!  (TAU(I), C(1,I),I=1,...,N) = ABSCISSAE AND ORDINATES OF THE DATA
!      POINTS.  TAU IS ASSUMED TO BE STRICTLY INCREASING.

!  IBCBEG, IBCEND = BOUNDARY CONDITION INDICATORS, AND
!  C(2,1), C(2,N) = BOUNDARY CONDITION INFORMATION.  SPECIFICALLY,
!      IBCBEG = 0 IMPLIES NO BOUNDARY CONDITION AT TAU(1) IS GIVEN.
!            IN THIS CASE, THE "NOT-A-KNOT" CONDITION IS USED, IE THE
!            JUMP IN THE 3-RD DERIVATIVE ACROSS TAU(2) IS FORCED TO 0.,
!            THUS THE 1-ST AND 2-ND POLYNOMIAL PIECES ARE MADE TO
!            COINCIDE.
!      IBCBEG = 1 IMPLIES THAT THE SLOPE AT TAU(1) IS SET EQUAL TO C(2,1)
!            INPUT BY USER.
!      IBCBEG = 2 IMPLIES THAT THE 2-ND DERIVATIVE AT TAU(1) IS SET EQUAL
!            TO C(2,1), SUPPLIED BY INPUT.
!      IBCEND = 0, 1, OR 2 HAS ANALOGOUS MEANING CONCERNING THE BOUNDARY
!            CONDITION AT TAU(N), WITH INFORMATION SUPPLIED BY C(2,N).

!  NDIM = ROW DIMENSION OF C MATRIX:  C(4,NDIM)

! **************************  O U T P U T  ****************************

!  C(J,I), J=1,...,4;  I=1,...,L=N-1  =  THE POLY COEFFS OF THE CUBIC
!      SPLINE WITH INTERIOR KNOTS TAU(2),...,TAU(N-1).  PRECISELY, IN THE
!      INTERVAL (TAU(I), TAU(I+1)), THE SPLINE F IS GIVEN BY

!       F(X) = C(1,I) + H*(C(2,I) + H*(C(3,I)/2. + H*C(4,I)/6.))

!      WHERE H = X - TAU(I).

!     THE COEFFICIENTS CALCULATED ARE, 1) THE VALUE, 2) THE SLOPE, AND
!     3) THE CURVATURE AT EACH OF THE KNOTS 1 TO N-1, AND 4) THE RATE OF
!     CHANGE OF THE CURVATURE OVER THE FOLLOWING INTERVAL.  IN ADDITION,
!     WE HAVE THE VALUE AND THE SLOPE AT THE LAST POINT. THE LAST TWO
!     POSTIONS AT THE LAST POINT ARE THEN SET TO THE CURVATURE AT THAT
!     POINT (IN C(3,N)) AND THE MEAN VALUE OVER THE ENTIRE INTERVAL,
!     CALCULATED AS THE INTEGRAL OVER THE INTERVAL DIVIDED BY THE LENGTH
!     OF THE INTERVAL (IN C(4,N)).

! **********************************************************************

      IMPLICIT REAL (A-H,O-Z)
      REAL ::    TAU(N)
      REAL :: C(4,NDIM), G, DTAU, DIVDF1, DIVDF3

      L = N - 1

      DO 10 M = 2,N
          C(3,M) = TAU(M) - TAU(M-1)
          C(4,M) = (C(1,M) - C(1,M-1)) / C(3,M)
   10 CONTINUE

!   * BEGINNING BOUNDARY CONDITION SECTION *

      IF (IBCBEG.EQ.0)  THEN   ! IBCBEG = 0
          IF (N.GT.2)  THEN    !     N > 2
              C(4,1) = C(3,3)
              C(3,1) = C(3,2) + C(3,3)
              C(2,1) = ((C(3,2) + 2.0*C(3,1))*C(4,2)*C(3,3) + &
                             C(3,2)**2 * C(4,3)) / C(3,1)
            ELSE                                        !     N = 2
              C(4,1) = (1.0,0.0)
              C(3,1) = (1.0,0.0)
              C(2,1) = 2.0 * C(4,2)
          END IF

      ELSE IF (IBCBEG.EQ.1)  THEN   ! IBCBEG = 1
          C(4,1) = (1.0,0.0)
          C(3,1) = (0.0,0.0)

      ELSE IF (IBCBEG.EQ.2)  THEN   ! IBCBEG = 2
          C(4,1) = (2.0,0.0)
          C(3,1) = (1.0,0.0)
          C(2,1) = 3.0*C(4,2) - C(2,1)*C(3,2)/2.0
      END IF

!   * RUNNING CALCULATIONS TO N-1 - LOOP IS NOT EXECUTED IF N = 2 *

      DO 20 M = 2,L
          G = -C(3,M+1) / C(4,M-1)
          C(2,M) = G*C(2,M-1) + 3.0*(C(3,M)*C(4,M+1) + C(3,M+1)*C(4,M))
          C(4,M) = G*C(3,M-1) + 2.0*(C(3,M) + C(3,M+1))
   20 CONTINUE

!   * ENDING BOUNDARY CONDITION SECTION *

      IF (IBCEND /= 1)  THEN
          IF (IBCEND.EQ.0)  THEN
                IF (N.EQ.2 .AND. IBCBEG.EQ.0)  THEN
                  C(2,N) = C(4,N)
              ELSE IF ((N.EQ.3 .AND. IBCBEG.EQ.0) .OR. N.EQ.2)  THEN
                  C(2,N) = 2.0 * C(4,N)
                  C(4,N) = (1.,0.)
                  G = -1.0 / C(4,N-1)
              ELSE
                  G = C(3,N-1) + C(3,N)
                  C(2,N) = ((C(3,N) + 2.0*G) * C(4,N)*C(3,N-1) + &
                        C(3,N)**2 * (C(1,N-1)-C(1,N-2)) / C(3,N-1)) / G
                  G = -G / C(4,N-1)
                  C(4,N) = C(3,N-1)
              END IF

          ELSE IF (IBCEND.EQ.2)  THEN
              C(2,N) = 3.0 * C(4,N) + C(2,N)*C(3,N)/2.0
              C(4,N) = (2.0,0.0)
              G = -1.0 / C(4,N-1)
          END IF

          IF ( IBCBEG.GT.0 .OR. N.GT.2)  THEN
              C(4,N) = G*C(3,N-1) + C(4,N)
              C(2,N) = (G*C(2,N-1) + C(2,N)) / C(4,N)
          END IF
      END IF

!   * RUN THE ENDING BOUNDARY EFFECT BACK THROUGH THE EQUATIONS *

      DO 40 J = L,1,-1
          C(2,J) = (C(2,J) - C(3,J)*C(2,J+1)) / C(4,J)
   40 CONTINUE

!   * FINAL CALCULATIONS *

      DO 50 I = 2,N
          DTAU = C(3,I)
          DIVDF1 = (C(1,I)-C(1,I-1)) / DTAU
          DIVDF3 = C(2,I-1) + C(2,I) - 2.0*DIVDF1
          C(3,I-1) = 2.0 * (DIVDF1-C(2,I-1)-DIVDF3) / DTAU
          C(4,I-1) = (DIVDF3/DTAU) * (6.0/DTAU)
   50 CONTINUE

!     * ADD THE CURVATURE AT THE LAST POINT IN THE THIRD POSITION OF THE 
!       LAST NODE *

      C(3,N) = C(3,L) + (TAU(N)-TAU(L)) * C(4,L)


!     * ADD THE MEAN VALUE OF THE ENTIRE INTERVAL IN THE FOURTH POSITION OF 
!       THE LAST NODE.  MEAN VALUE IS CALCULATED AS THE INTEGRAL OVER THE
!       INTERVAL DIVIDED BY THE LENGTH OF THE INTERVAL. *

      C(4,N) = (0.0,0.0)
      DO 60 I = 1,L   ! INTEGRATE OVER THE INTERVAL
          DTAU = TAU(I+1) - TAU(I)
          C(4,N) = C(4,N) + DTAU*(C(1,I) + DTAU*(C(2,I)/2.0 +     &
                            DTAU*(C(3,I)/6.0 + DTAU*C(4,I)/24.0)))
   60 CONTINUE
      C(4,N) = C(4,N) / (TAU(N)-TAU(1))   ! DIVIDE BY LENGTH OF INTERVAL

      RETURN
      END



      SUBROUTINE VSPLINE (TAU, C, M, MDIM, F, N)

!     VSPLINE CALCULATES THE CUBIC SPLINE VALUES FOR A SET OF N POINTS  
!     IN F FROM THE M-POINT CUBI! SPLINE FIT IN ! AND THE NODES IN TAU.
!     THE POINTS ARE RETURNED IN F.  ALL OF THE POINTS IN F MUST LIE
!     BETWEEN TAU(1) AND TAU(M).

!     * * * * * * * * * * * * *   WARNINGS   * * * * * * * * * * * * *

!     POINTS OUTSIDE OF THE SPLINE FIT REGION ARE EXTRAPOLATED FROM THE END
!     INTERVALS.  THIS CAN RESULT IN WILD VALUES IF EXTRAPOLATED TOO FAR.
!     ALSO THE POINTS MUST BE IN STRICTLY ASCENDING ORDER, IF NOT THE
!     POINTS WHICH ARE OUT OF ORDER WILL BE EXTRAPOLATED FROM THE CURRENT
!     INTERVAL AGAIN RESULTING IN WILD VALUES.

      IMPLICIT REAL (A-H,O-Z)
      REAL ::     TAU(M)
      REAL :: C(4,MDIM), F(N), SPLINE

      J = 1
      DO 20 I = 1,N
   10 J1 = J + 1
      IF (TAU(J1).LT.REAL(F(I)) .AND. J1.LT.M) THEN ! CHECK TO MAKE SURE
           J = J + 1   ! THIS POINT IS NOT
           GO TO 10    ! IN THE NEXT INTERVAL.
      END IF
      H = DBLE (F(I)) - TAU(J)   ! DISTANCE FROM START OF INTERVAL
      F(I) = SPLINE (C(1,J), H)
   20 CONTINUE
      RETURN
      END
!**********************************************************************C
      FUNCTION SPLINE ( C, H )

!     THIS FUNCTION EVALUATES THE SPLINE AT THE POINT H

      IMPLICIT REAL ( A-H, O-Z )
      REAL C(4), SPLINE

      SPLINE = C(1) + H * ( C(2) + H * ( C(3) / 2.0 + H * C(4) / 6.0 ) )
      RETURN
      END

      FUNCTION SPLINEX ( C, H )

!     THIS FUNCTION EVALUATES THE SPLINE DERIVATIVE AT THE POINT H

      IMPLICIT REAL ( A-H, O-Z )
      REAL :: C(4), SPLINEX

      SPLINEX = C(2) + H * ( C(3) + H * C(4) / 2.0 )
      RETURN
      END

      FUNCTION SPLINEXX ( C, H )

!     THIS FUNCTION EVALUATES THE SPLINE 2ND DERIVATIVE AT THE POINT H

      IMPLICIT REAL ( A-H, O-Z )
      REAL :: C(4), SPLINEXX

      SPLINEXX = C(3) + H * C(4)
      RETURN
      END

!**********************************************************************C

      SUBROUTINE SPLINEALL ( C, H, F, FX, FXX )

!     THIS ROUTINE EVALUATES THE
!        SPLINE,
!        SPLINE DERIVATIVE, AND
!        SPLINE 2ND DERIVATIVE AT THE POINT H

      IMPLICIT REAL ( A-H, O-Z )
      PARAMETER ( HALF = 0.5, SIXTH = 1.0 / 6.0 )
      REAL :: C(4), F, FX, FXX

      F   = C(1) + H * ( C(2) + H * ( HALF * C(3) + SIXTH * H * C(4) ) )
      FX  = C(2) + H * ( C(3) + H * HALF * C(4) )
      FXX = C(3) + H * C(4)

      RETURN
      END
