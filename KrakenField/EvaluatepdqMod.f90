MODULE EvaluatePDQMod
  IMPLICIT NONE

CONTAINS
  SUBROUTINE EvaluatePDQ( freq, k, phiR, phiS, M, maxM, IelementSource, xs, ys, theta, Ntheta, RminM, RmaxM, Nr, MLimit, P )

    ! Computes 3-D pressure field using adiabatic mode theory.           
    ! Normalized to pressure of point source at 1 meter.                 
    ! Note RminM must be zero.
    ! Compared to the standard version (Evaluate3d), this one uses the wavenumber for the first mode to scale all the modes.
    ! That approximation is poor if the field is dominated by steep angle paths.
    ! This often happens, for instance, for high frequency wind noise calculations near the equator
    ! mike porter, 1988

    USE ElementMod
    USE BeamPattern
    USE interpolation
    USE Evaluate3DMod

    INTEGER, INTENT( IN ) :: M( * ), MLimit, maxM             ! number of modes, limit on number of modes to propagate
    INTEGER, INTENT( IN ) :: Nr, Ntheta, IElementSource       ! number of receiver ranges and bearing lines
    REAL (KIND=8), INTENT( IN ) :: xs, ys                     ! source coordinate
    REAL (KIND=8), INTENT( IN ) :: freq                       ! source frequency
    REAL (KIND=8), INTENT( IN ) :: RminM, RmaxM               ! minimum and maximum receiver ranges in m
    REAL (KIND=8), INTENT( IN ) :: theta( Ntheta )            ! bearing angles for receivers
    COMPLEX, INTENT( IN ) :: k( maxM, * )                     ! wavenumbers
    COMPLEX, INTENT( IN ) :: PhiR( maxM, * ), PhiS( maxM, * ) ! source/receiver mode shapes
    COMPLEX, INTENT( OUT) :: P( Ntheta, Nr )                  ! pressure field
    INTEGER               :: Ielement, ir, Itheta, MProp, NewElement, Outside
    REAL (KIND=8)         :: delta_r, Rin, Rout, RM
    COMPLEX               :: PhiIn( maxM ), PhiOut( maxM ), const( maxM )
    COMPLEX               :: kIn( maxM ),  kOut( maxM )
    COMPLEX               :: kAverage( maxM )
    COMPLEX (KIND=8)      :: PhaseInc( maxM )
    REAL (KIND=8), ALLOCATABLE :: kz2( : ), thetaT( : ), S( : )
    REAL (KIND=8)         :: omega

    P       = 0    ! zero out the pressure field
    delta_r = ( RmaxM - RminM ) / ( Nr - 1 )

    !  *** Loop over angle ***                                           

    Bearing: DO Itheta = 1, Ntheta
       tsx = COS( DegRad * theta( Itheta ) )
       tsy = SIN( DegRad * theta( Itheta ) )

       ! Get modal values                                          
       Ielement = IelementSource
       ! WRITE( *, * )
       ! WRITE( *, * ) 'Itheta, tsx, tsy', Itheta, tsx, tsy

       CALL SourceElement( Ielement, outside, Rin, Rout, xs, ys, Mprop, M, maxM, &
            k, phiR, phiS, const, kin, phiin, kout, phiout )
       Mprop = MIN( MLimit, Mprop )

       IF ( MProp > 0 ) THEN   ! any propagating modes at the source?

          ! apply the source beam pattern
          ! using kIn for wavenumber; would be more precise to use interpolated value at source
          IF ( SBPFlag == '*' .AND. Itheta == 1 ) THEN
             ALLOCATE( kz2( MProp ), thetaT( MProp ), S( MProp ) )
             omega = 2 * pi * freq
             kz2   = REAL( omega ** 2 / c0 ** 2 - kIn( 1 : MProp ) ** 2 )      ! vertical wavenumber squared
             WHERE ( kz2 < 0 ) kz2 = 0                                         ! remove negative values

             thetaT = RadDeg * ATAN( sqrt( kz2 ) / REAL( kIn( 1 : MProp ) ) )  ! calculate the angle in degrees
             CALL interp1( SrcBmPat( :, 1 ), SrcBmPat( :, 2 ), thetaT, S )
             const( 1 : MProp ) = const( 1 : MProp ) * REAL( S )               ! apply the shading
          END IF

          const(    1 : Mprop ) = i * SQRT( 2.0 * pi ) * EXP( i * pi / 4.0 ) * const( 1 : Mprop )                         
          kAverage( 1 : Mprop ) = 0.5 * ( kin( 1 : Mprop ) + kout( 1 : Mprop ) ) 
          PhaseInc( 1 : Mprop ) = EXP( -i * kAverage( 1 : Mprop ) * delta_r )

          RangeStepping: DO ir = 1, Nr ! March forward in range 
             RM = RminM + ( ir - 1 ) * delta_r 
             IF ( RM == 0.0 ) RM = MIN( 1.0D0, delta_r )

             ! Crossing into new element?                                  
             EltLoop: DO WHILE ( RM > Rout )

                ! Copy outside info to inside                          
                NewElement = AdjacentElement( outside, Ielement )
                ! IF ( NewElement == 0 ) EXIT RangeStepping   ! jump out if there is no adjacent element

                Rin = Rout
                kin(      1 : Mprop ) = kout(   1 : Mprop ) 
                phiin(    1 : Mprop ) = phiout( 1 : Mprop ) 
                kAverage( 1 : Mprop ) = 0.5 * ( kin( 1 : Mprop ) + kout( 1 : Mprop ) ) 
                PhaseInc( 1 : Mprop ) = EXP( -i * kAverage( 1 : Mprop ) * delta_r )

                ! Get new outside info                                 
                CALL OUT( Ielement, NewElement, outside, Rout, xs, ys, Mprop, M, maxM, k, phiR, kout, phiout )

                IF ( MProp <= 0 ) EXIT RangeStepping   ! jump out if there are no propagating modes
                Ielement = NewElement 

             END DO EltLoop

             ! Compute modal contribution at this range

             const( 1 : Mprop ) = CMPLX( const( 1 : Mprop ) * PhaseInc( 1 : Mprop ) )   ! advance the phase
             P( Itheta, ir )    = CMPLX( SUM( const( 1 : Mprop ) * phiin( 1 : Mprop ) ) / SQRT( RM * kin( 1 ) ) )

          END DO RangeStepping
       END IF
    END DO Bearing

    IF ( ALLOCATED( kz2 ) ) DEALLOCATE( kz2, thetaT, S )

  END SUBROUTINE EvaluatePDQ
END MODULE EvaluatePDQMod
