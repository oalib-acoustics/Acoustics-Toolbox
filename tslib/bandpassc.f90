SUBROUTINE BandPass( X, N, deltat, fmin, fmax )

  ! Band-pass filter of a complex time series, X( 1 : N )
  ! deltat is the time step
  ! fmin, fmax are the frequencies of the band
  
  ! This version assumes time series was real and therefore filters symmetically

  REAL, PARAMETER :: PI = 3.141592, alpha = 0.54
  COMPLEX, INTENT( INOUT ) :: X( N )

  ! Compute indices for frequency band

  deltaf = 1.0 / ( N * deltat )
  iMin   = INT( fmin / deltaf )
  iMax   = INT( fmax / deltaf )

  ! Quick exit if full band is passed:
  IF ( iMin <= 1 .AND. iMax >= N / 2 + 1 ) RETURN

  ! Check N is positive and a power of 2
  IF ( N <= 0 ) STOP 'FATAL ERROR in BandPass: N must be positive'
  
  NT = 2**( INT( LOG10( REAL( N ) ) / 0.30104 ) + 1 )
  IF ( NT /= N ) STOP 'FATAL ERROR in BandPass: N must be a power of 2'

  CALL CFFT( X, N, 1 )   ! forward transform
  X = X / N              ! scaling ...

  ! Zero out-of-band components (Hanning window)

  DO I = 1, N / 2 + 1

     IF ( I >= iMin .AND. I <= iMax ) THEN
        weight = alpha + ( 1.0 - alpha ) * COS( PI * ( I - 1 ) / ( iMax - 1 ) )
        weight = 1.0   ! disable windowing
     ELSE
        weight = 0.0
     ENDIF

     X( I ) = weight * X( I )
     IF ( I > 1 .AND. I < N/2 + 1 ) X( N + 2 - I ) = weight * X( N + 2 - I )
  END DO

  CALL CFFT( X, N, -1 )   ! inverse transform

END SUBROUTINE BandPass
