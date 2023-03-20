MODULE AttenMod

  ! Attenuation module
  ! Routines to convert a sound speed and attenuation in user units to a complex sound speed
  ! Includes a formula for volume attenuation

  USE FatalError
  IMPLICIT NONE
  INTEGER, PRIVATE, PARAMETER      :: PRTFile = 6
  INTEGER, PARAMETER               :: MaxBioLayers = 200
  INTEGER                          :: iBio, NBioLayers
  REAL    (KIND=8) :: T = 20, Salinity = 35, pH = 8, z_bar = 0, FG   ! Francois-Garrison volume attenuation; temperature, salinity, ...

  TYPE bioStructure
     REAL (KIND=8) :: Z1, Z2, f0, Q, a0
  END TYPE bioStructure

  TYPE( bioStructure ) :: bio( MaxBioLayers )

CONTAINS
  FUNCTION CRCI( z, c, alpha, freq, freq0, AttenUnit, beta, fT )

    ! Converts real wave speed and attenuation to a single
    !  complex wave speed (with positive imaginary part)

    ! AttenUnit
    ! 6 Cases:    N Nepers/meter
    !             M dB/meter      (M for Meters)
    !             m dB/meter with a power law
    !             F dB/m-kHz      (F for frequency dependent)
    !             W dB/wavelength (W for Wavelength)
    !             Q Q
    !             L Loss parameter
    !
    ! second letter adds volume attenuation according to standard laws:
    !             T for Thorp
    !             F for Francois Garrison
    !             B for biological
    !
    ! freq is the current frequency
    ! freq0 is the reference frequency for which the dB/meter was specified (used only for 'm')

    ! Returns
    ! c     real      part of sound speed
    ! alpha imaginary part of sound speed

    USE MathConstants

    REAL     (KIND=8), INTENT( IN )  :: freq, freq0, alpha, c, z, beta, fT
    CHARACTER (LEN=2), INTENT( IN )  :: AttenUnit
    REAL     (KIND=8)                :: f2, omega, alphaT, Thorp, a, FG
    COMPLEX  (KIND=8)                :: CRCI

    omega = 2.0 * pi * freq

    !  Convert to Nepers/m 
    alphaT = 0.0
    SELECT CASE ( AttenUnit( 1 : 1 ) )
    CASE ( 'N' )
       alphaT = alpha
    CASE ( 'M' )   ! dB/m
       alphaT = alpha / 8.6858896D0
    CASE ( 'm' )   ! dB/m with power law
       alphaT = alpha / 8.6858896D0
       IF ( freq < fT ) THEN   ! frequency raised to the power beta
          alphaT = alphaT * ( freq / freq0 ) ** beta
       ELSE                    ! linear in frequency
          alphaT = alphaT * ( freq / freq0 ) * ( fT / freq0 ) ** ( beta - 1 )
       END IF
    CASE ( 'F' )   ! dB/(m kHz)
       alphaT = alpha * freq / 8685.8896D0
    CASE ( 'W' )   ! dB/wavelength
       IF ( c /= 0.0         ) alphaT = alpha * freq / ( 8.6858896D0 * c )
       !        The following lines give f^1.25 frequency dependence
       !        FAC = SQRT( SQRT( freq / 50.0 ) )
       !        IF ( c /= 0.0 ) alphaT = FAC * alpha * freq / ( 8.6858896D0 * c )
    CASE ( 'Q' )   ! Quality factor
       IF ( c * alpha /= 0.0 ) alphaT = omega / ( 2.0 * c * alpha )
    CASE ( 'L' )   ! loss parameter
       IF ( c /= 0.0         ) alphaT = alpha * omega / c
    END SELECT

    ! added volume attenuation
    SELECT CASE ( AttenUnit( 2 : 2 ) )
    CASE ( 'T' )
       f2 = ( freq / 1000.0 ) ** 2

       ! Original formula from Thorp 1967
       ! Thorp = 40.0 * f2 / ( 4100.0 + f2 ) + 0.1 * f2 / ( 1.0 + f2 )   ! dB/kyard
       ! Thorp = Thorp / 914.4D0                 ! dB / m
       ! Thorp = Thorp / 8.6858896D0             ! Nepers / m

       ! Updated formula from JKPS Eq. 1.34
       Thorp = 3.3d-3 + 0.11 * f2 / ( 1.0 + f2 ) + 44.0 * f2 / ( 4100.0 + f2 ) + 3d-4 * f2   ! dB/km
       Thorp = Thorp / 8685.8896 ! Nepers / m

       alphaT = alphaT + Thorp
    CASE ( 'F' )   ! Francois-Garrison
       FG     = Franc_Garr( freq / 1000 );   ! dB/km
       FG     = FG / 8685.8896;                           ! Nepers / m
       alphaT = alphaT + FG
    CASE ( 'B' )   ! biological attenuation per Orest Diachok
       DO iBio = 1, NBioLayers
          IF ( z >= bio( iBio )%Z1 .AND. z <= bio( iBio )%Z2 ) THEN
             a = bio( iBio )%a0 / ( ( 1.0 - bio( iBio )%f0 ** 2 / freq ** 2  ) ** 2 + 1.0 / bio( iBio )%Q ** 2 )   ! dB/km
             a = a / 8685.8896   ! Nepers / m
             alphaT = alphaT + a
          END IF
       END DO
    END SELECT

    ! Convert Nepers/m to equivalent imaginary sound speed 
    alphaT = alphaT * c * c / omega
    CRCI   = CMPLX( c, alphaT, KIND=8 )

    IF ( alphaT > c ) THEN
       WRITE( PRTFile, * ) 'Complex sound speed: ', CRCI
       WRITE( PRTFile, * ) 'Usually this means you have an attenuation that is way too high'
       CALL ERROUT( 'AttenMod : CRCI ', 'The complex sound speed has an imaginary part > real part' )
    END IF

    RETURN
  END FUNCTION CRCI

  !**********************************************************************!

  FUNCTION Franc_Garr( f )

    ! Francois Garrison formulas for attenuation
    ! Based on a Matlab version by D. Jackson APL-UW

    ! mbp Feb. 2019
    ! Verified using F-G Table IV

    ! alpha = attenuation   (dB/km)
    ! f     = frequency     (kHz)
    ! T     = temperature   (deg C)
    ! S     = salinity      (psu)
    ! pH    = 7 for neutral water
    ! z_bar = depth         (m)

    !     Returns
    !        alpha = volume attenuation in dB/km

    REAL (KIND=8) :: f, Franc_Garr
    REAL (KIND=8) :: c, A1, A2, A3, P1, P2, P3, f1, f2

    c = 1412 + 3.21 * T + 1.19 * Salinity + 0.0167 * z_bar

    ! Boric acid contribution
    A1 = 8.86 / c * 10 ** ( 0.78 * pH - 5 )
    P1 = 1
    f1 = 2.8 * sqrt( Salinity / 35 ) * 10 ** ( 4 - 1245 / ( T + 273 ) )

    ! Magnesium sulfate contribution
    A2 = 21.44 * Salinity / c * ( 1 + 0.025 * T )
    P2 = 1 - 1.37D-4 * z_bar + 6.2D-9 * z_bar ** 2
    f2 = 8.17 * 10 ** ( 8 - 1990 / ( T + 273 ) ) / ( 1 + 0.0018 * ( Salinity - 35 ) )

    ! Viscosity
    P3 = 1 - 3.83D-5 * z_bar + 4.9D-10 * z_bar ** 2
    if ( T < 20 ) THEN
       A3 = 4.937D-4 - 2.59D-5 * T + 9.11D-7 * T ** 2 - 1.5D-8  * T ** 3
    else
       A3 = 3.964D-4 -1.146D-5 * T + 1.45D-7 * T ** 2 - 6.5D-10 * T ** 3
    end if

    Franc_Garr = A1 * P1 * ( f1 * f ** 2 ) / ( f1 ** 2 + f ** 2 ) + A2 * P2 * ( f2 * f ** 2 ) / ( f2 ** 2 + f ** 2 ) + &
         A3 * P3 * f ** 2

  END FUNCTION Franc_Garr

END MODULE AttenMod
