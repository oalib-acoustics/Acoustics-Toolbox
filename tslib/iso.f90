PROGRAM ISO

  ! Generates the time series for an isovelocity problem

  INTEGER STKFIL
  PARAMETER ( STKFIL = 6, PI = 3.141592 )

  REAL      RD( 51 ), STACK( 51, 1001 ), SD( 1 )
  COMPLEX   ST
  CHARACTER PULSE*4, PULSTITL*80

  C0    = 1500.0   ! Reference soundspeed
  PULSE = 'P   '   ! Pulse type
  FREQ  = 50.0
  FMIN  =   0.0
  FMAX  =  10.0 * FREQ
  OMEGA =   2.0 * PI * FREQ
  R     =  200.0   ! Range in meters

  TMIN = 0.0
  TMAX = 1.0
  NT   = 512
  DELT = TMAX / ( NT - 1 )

  ! Generate the vector of source/receiver depths

  NSD = 1
  SD( 1 ) = 0

  NRD = 11
  RD = [ ( IRD, IRD = 1, NRD ) ] * 10.0

  ! Construct the time series

  DO IR = 1, NRD
     SLANTR = SQRT( RD( IR ) ** 2 + R ** 2 )
     DO IT = 1, NT
        TIME = TMIN + ( IT - 1 ) * DELT
        T = TIME - SLANTR / C0    ! Reduced time

        !  calculate the source time series

        CALL SOURCE( T, ST, SD, NSD, NT, OMEGA, FMIN, FMAX, PULSE, PULSTITL )

        STACK( IR, IT ) = -ST / SLANTR
     END DO   ! Next TIME
  END DO  ! Next depth

  ! Write out the stack

  WRITE( STKFIL, * ) ''''//'Analytic- '//PULSTITL(1:69)//''''
  WRITE( STKFIL, * ) NRD, ( RD( IR ), IR = 1, NRD )

  DO IT = 1, NT
     TIME = TMIN + ( IT - 1 ) * DELT
     WRITE( STKFIL, * ) TIME, ( STACK( IR, IT ), IR = 1, NRD )
  END DO

  STOP
END PROGRAM ISO
