SUBROUTINE PREENV( X, N )

  ! Forms the pre-envelope of the function
  ! Real( X( N ) ) is the input time series
  ! The output vector is complex
  ! N must be a power of 2 and <16384

  COMPLEX X( N )

  !  Check N is postive and a power of 2

  IF ( N <= 0 ) STOP 'FATAL ERROR in PREENV: N must be positive'

  NT = 2**( INT( LOG10( REAL( N ) ) / 0.30104 ) + 1 )
  IF ( NT /= N ) STOP 'FATAL ERROR in PREENV: N must be a power of 2'

  CALL CFFT( X, N, 1 )     ! forward Fourier transform
  X = X / N                ! scaled appropriately

  IMID            = N / 2
  X( IMID+1 : N ) = 0.0    ! Zero-out the negative spectrum (the upper N/2 positions)

  CALL CFFT( X, N, -1 )    ! inverse Fourier transform

END SUBROUTINE PREENV
