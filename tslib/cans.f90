SUBROUTINE CANS( T, omega, Pulse, SF, NSD, IT, PulseTitle )

  !     Computes the source time series

  !     T      is the time
  !     omega  is some frequency characterizing the pulse
  !     Pulse  is the letter indicating which pulse
  !     S      is the time series

  IMPLICIT NONE
  REAL, PARAMETER :: pi = 3.141592

  INTEGER,   INTENT( IN  ) :: NSD, IT
  REAL,      INTENT( IN  ) :: T, omega
  CHARACTER, INTENT( IN  ) :: Pulse*1
  COMPLEX,   INTENT( OUT ) :: SF( NSD, * )
  CHARACTER, INTENT( OUT ) :: PulseTitle*( * )
  INTEGER  :: NSIG
  REAL     :: S, F, A, T0, TS, TC, U

  S = 0.0
  F = omega / ( 2.0 * pi )

  IF ( T > 0.0 ) THEN    ! Evaluate selected Pulse
     SELECT CASE ( Pulse )
     CASE ( 'P' )   ! Pseudo gaussian
        ! (peak at 0, support [0, 3F]
        IF ( T <= 1.0 / F ) S = 0.75 - COS( omega * T ) + 0.25*COS( 2.0 * omega * T )
        PulseTitle = 'Pseudo gaussian'

     CASE ( 'R' )   ! Ricker wavelet
        ! (peak at F, support [0, 2F]
        U = omega * T - 5.0
        S = 0.5 * ( 0.25*U*U - 0.5 ) * SQRT( pi ) * EXP( -0.25*U*U )
        PulseTitle = 'Ricker wavelet'

     CASE ( 'A' )    ! Approximate Ricker wavelet
        ! (peak at F, support [0, 2.5F]
        TC = 1.55 / F
        IF ( T <= TC ) &
             &   S = 0.48829 *    COS( 2.0 * pi * T / TC ) &
             &      -0.14128 * 4* COS( 4.0 * pi * T / TC ) &
             &      +0.01168 * 9* COS( 6.0 * pi * T / TC )
        PulseTitle = 'Approximate Ricker wavelet'

     CASE ( 'S' )   ! Single sine
        ! (peak at F, support [0, infinity], nulls at nF
        IF ( T <= 1.0 / F ) S = SIN( omega * T )
        PulseTitle = 'Single sine'

     CASE ( 'H' )   ! Hanning weighted four sine
        ! (peak at F, support [0, infinity], first null at about 1.5F
        IF ( T <= 4.0 / F ) S = 0.5 * SIN( omega * T ) * ( 1 - COS( omega * T / 4.0 ) )
        PulseTitle = 'Hanning weighted four sine'

     CASE ( 'N' )   ! N-wave
        ! (peak at 1 F, support [0, 4F], [0,3F] OK
        IF ( T <= 1.0 / F ) S = SIN( omega * T ) - 0.5 * SIN( 2.0 * omega * T )
        PulseTitle = 'N-wave'

     CASE ( 'M' )   ! Miracle wave (has a Hilbert transform?)
        ! (peak at 0, support [0, infinity]
        A  = 1.0 / ( 6.0 * F )
        T0 = 6.0 * A
        TS = ( T - T0 ) / A
        S  = 1.0 / ( 1.0 + TS * TS )
        ! HS is the Hilbert transform of the time series
        ! HS = TS / ( 1.0 + TS * TS )
        PulseTitle = 'Miracle wave'

     CASE ( 'G' )   ! Gaussian
        ! (peak at 0, support [0, infinity]
        ! T0 is the peak time
        ! A is the 3dB down time
        ! NOTE S(0) = EXP( -NSIG ** 2 )

        NSIG = 3
        A    = 1.0 / F / ( 2.0 * NSIG )
        T0   = NSIG * A
        S    = EXP( -( ( T - T0 ) / A ) ** 2 )
        PulseTitle = 'Gaussian'
     CASE ( 'T' )   ! Tone
        ! (peak at F, support [0, infinity]

        IF ( T <= 0.4 ) S = SIN( omega * T )
        PulseTitle = 'Tone'
     CASE ( 'C' )   ! Sinc
        ! ( uniform spectrum from [0, F]

        S = SIN( omega * T ) / ( omega * T )
        PulseTitle = 'Sinc'
     END SELECT
  END IF

  SF( :, IT ) = S ! Copy into source vector (waveform duplicated for every source depth)

END SUBROUTINE CANS
