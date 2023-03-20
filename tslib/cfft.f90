SUBROUTINE CFFT( DATA, N, IFORW )

  ! complex FFT

  IMPLICIT NONE
  INTEGER I1, I2A, I2B, I3, I3Rev, IP1, IP2, N, ISign, IForw
  REAL    theta, sinth
  COMPLEX DATA( N ), TEMP, W, WSTP

  ISIGN = -IFORW
  I3REV = 1

  DO I3 = 1, N
     IF ( I3 < I3REV ) THEN   ! switch values
        TEMP          = DATA( I3 )
        DATA( I3 )    = DATA( I3REV )
        DATA( I3REV ) = TEMP
     ENDIF
     ! following loop is just to compute I3REV
     IP1 = N / 2
     DO WHILE (  I3REV > IP1 )
        IF ( IP1 <= 1 ) EXIT
        I3REV = I3REV - IP1
        IP1   = IP1 / 2
     END DO
     I3REV = I3REV + IP1
  END DO

  IP1 = 1

  DO WHILE  ( IP1 < N )
     IP2   = IP1 * 2
     THETA = 6.283185307 / FLOAT( ISIGN * IP2 )
     SINTH = SIN( THETA / 2.)
     WSTP  = CMPLX( -2. * SINTH * SINTH, SIN( THETA ) )
     W     = 1.

     DO I1 = 1, IP1
        DO I3 = I1, N, IP2
           I2A = I3
           I2B = I2A + IP1
           TEMP      = W * DATA( I2B )
           DATA( I2B ) = DATA( I2A ) - TEMP
           DATA( I2A ) = DATA( I2A ) + TEMP
        END DO

        W = W * WSTP + W
     END DO

     IP1 = IP2
  END DO

  RETURN
END SUBROUTINE CFFT
