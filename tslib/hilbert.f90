SUBROUTINE Hilbert( X, N )

  ! Forms the Hilbert transform of the function
  ! Real( X( N ) ) is the input time series
  ! N must be a power of 2 and <16384

  COMPLEX, PARAMETER :: CI = ( 0.0, 1.0 )
  COMPLEX X( N )

  ! Check N is positve and a power of 2

  IF ( N <= 0 ) STOP 'FATAL ERROR in HILBERT: N must be positive'

  NT = 2**( INT( LOG10( REAL( N ) ) / 0.30104 ) + 1 )
  IF ( NT /= N ) STOP 'FATAL ERROR in HILBERT: N must be a power of 2'

  ! Fourier transform

  CALL CFFT( X, N, 1 )
  X = X / N   ! scaling ...

  !  Multiply by i * sgn( f )

  IMID = N / 2
  X( 1 : IMID - 1 ) =  CI * X( 1 : IMID - 1 )   ! pos. spectrum
  X( IMID     )     = 0.0                       ! d.c. component
  X( IMID + 1 : N ) = -CI * X( IMID + 1 : N )   ! neg. spectrum

  ! Inverse Fourier transform

  CALL CFFT( X, N, -1 )

END SUBROUTINE Hilbert
