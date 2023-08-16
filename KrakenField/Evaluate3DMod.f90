MODULE Evaluate3DMod

  USE ElementMod 
  IMPLICIT NONE
  REAL,    PARAMETER    :: pi = 3.141592, DegRad = pi / 180.0, RadDeg = 180.0 / pi, c0 = 1500   ! reference sound speed, should be speed at the source depth
  REAL (KIND=8)         :: tsx, tsy                        ! tangent vector from source
  COMPLEX, PARAMETER    :: i  = ( 0.0, 1.0 )

CONTAINS
  SUBROUTINE Evaluate3D( freq, k, PhiR, PhiS, M, maxM, IElementSource, xs, ys, theta, Ntheta, RminM, RmaxM, Nr, MLimit, P )

    ! Computes 3-D pressure field using adiabatic mode theory           
    ! Normalized to pressure of point source at 1 meter                 
    ! Note RminM must be zero
    ! mike porter, 1987

    USE BeamPattern
    USE interpolation

    INTEGER, PARAMETER    :: KBarFile = 55, ZBarFile = 56
    INTEGER, INTENT( IN ) :: M( * ), MLimit, maxM             ! number of modes, limit on number of modes to propagate
    INTEGER, INTENT( IN ) :: Nr, Ntheta, IElementSource       ! number of receiver ranges and bearing lines
    REAL (KIND=8), INTENT( IN ) :: xs, ys                     ! source coordinate
    REAL (KIND=8), INTENT( IN ) :: freq                             ! source frequency
    REAL (KIND=8), INTENT( IN ) :: RminM, RmaxM               ! minimum and maximum receiver ranges in m
    REAL (KIND=8), INTENT( IN ) :: theta( Ntheta )            ! bearing angles for receivers
    COMPLEX, INTENT( IN ) :: k( maxM, * )                     ! wavenumbers
    COMPLEX, INTENT( IN ) :: PhiR( maxM, * ), PhiS( maxM, * ) ! source/receiver mode shapes
    COMPLEX, INTENT( OUT) :: P( Ntheta, Nr )                  ! pressure field
    INTEGER               :: iElement, ir, itheta, L, MProp, NewElement, Outside
    REAL (KIND=8)         :: alpha, delta_r, Rin, Rout, RM
    COMPLEX               :: PhiIn( maxM ), PhiOut( maxM ), const( maxM ), PhiInt
    COMPLEX               :: kIn(   maxM ),   kOut( maxM ), T, kInt
    COMPLEX               :: sumk(  MaxM )
    REAL (KIND=8), ALLOCATABLE :: kz2( : ), thetaT( : ), S( : )
    REAL (KIND=8)         :: omega

    ! Open file for eigenfunctions                                  
    ! OPEN ( FILE = 'ZBarFile', UNIT = ZBarFile, STATUS = 'UNKNOWN' ) 

    P       = 0    ! zero out the pressure field
    delta_r = ( RmaxM - RminM ) / ( Nr - 1 )

    !  *** Loop over angle ***                                           

    Bearing: DO itheta = 1, Ntheta 
       tsx = COS( DegRad * theta( itheta ) )
       tsy = SIN( DegRad * theta( itheta ) ) 

       ! Get modal values
       iElement = IElementSource 
       ! WRITE( *, * ) 'Tracing bearing ', itheta, ATAN2( tsy, tsx ) / DegRad

       CALL SourceElement( iElement, Outside, RIn, ROut, xs, ys, MProp, M, maxM, &
            k, PhiR, PhiS, const, kIn, PhiIn, kOut, PhiOut )
       MProp = MIN( MLimit, MProp ) 

       IF ( MProp > 0 ) THEN   ! any propagating modes at the source?
          ! write modes at first range                                 
          ! IF ( itheta == 1 ) THEN 
          !    WRITE( ZBarFile, * ) MProp 
          !    WRITE( ZBarFile, * ) const( 1 : MProp ) 
          ! END IF

          ! apply the source beam pattern
          ! using kIn for wavenumber; would be more precise to use interpolated value at source
          IF ( SBPFlag == '*' .AND. itheta == 1 ) THEN
             ALLOCATE( kz2( MProp ), thetaT( MProp ), S( MProp ) )
             omega = 2 * pi * freq
             kz2   = REAL( omega ** 2 / c0 ** 2 - kIn( 1 : MProp ) ** 2 )      ! vertical wavenumber squared
             WHERE ( kz2 < 0 ) kz2 = 0                                         ! remove negative values

             thetaT = RadDeg * ATAN( sqrt( kz2 ) / REAL( kIn( 1 : MProp ) ) )  ! calculate the angle in degrees
             CALL interp1( SrcBmPat( :, 1 ), SrcBmPat( :, 2 ), thetaT, S )
             const( 1 : MProp ) = const( 1 : MProp ) * REAL( S )               ! apply the shading
          END IF

          const( 1 : MProp ) = i * SQRT( 2.0 * pi ) * EXP( i * pi / 4.0 ) * const( 1 : MProp )
          sumk = 0.0

          ! *** March forward in range ***                                 

          RangeStepping: DO ir = 1, Nr 
             RM = RminM + ( ir - 1 ) * delta_r 
             IF ( RM == 0.0 ) RM = MIN( 1.0D0, delta_r )

             ! Crossing into new element?                                  
             EltLoop: DO WHILE ( RM > ROut ) 

                ! Copy outside info to inside                          
                NewElement         = AdjacentElement( Outside, iElement ) 
                RIn                =   ROut 
                kIn(   1 : MProp ) =   kOut( 1 : MProp ) 
                PhiIn( 1 : MProp ) = PhiOut( 1 : MProp ) 

                ! Get new outside info                                 
                CALL OUT( iElement, NewElement, Outside, ROut, xs, ys, MProp, M, maxM, k, PhiR, kOut, PhiOut )
                IF ( MProp <= 0 ) EXIT RangeStepping   ! jump out if there are no propagating modes
                iElement = NewElement 

             END DO EltLoop

             ! *** Compute modal contribution at this range ***            

             T = 0.0 

             IF ( RIn /= ROut ) THEN 
                alpha = ( RM - RIn ) / ( ROut - RIn ) 
                alpha = MIN( MAX( alpha,  0.0D0 ), 1.0D0 )
             ELSE 
                alpha = 0.0 
             END IF

             ! IF ( ir == Nr ) WRITE( ZBarFile, * ) MProp
             Mode: DO L = 1, MProp 
                kInt      = CMPLX(   kIn( L ) + alpha * (   kOut( L ) -   kIn( L ) ) )
                PhiInt    = CMPLX( PhiIn( L ) + alpha * ( PhiOut( L ) - PhiIn( L ) ) )
                sumk( L ) = CMPLX(  sumk( L ) + delta_r * kInt )
                T         = T + PhiInt * const( L ) * EXP( -i * sumk( L ) ) / SQRT( kInt )    
                ! IF ( ir == Nr ) WRITE( ZBarFile, * ) PhiInt  ! write mode at last range     
             END DO Mode

             P( itheta, ir ) = CMPLX( T / SQRT( RM ) )

          END DO RangeStepping

          ! Write average wavenumber to file                           
          ! OPEN ( FILE = 'KBarFile', UNIT = KBarFile, STATUS = 'UNKNOWN' ) 
          ! WRITE( KBarFile, * ) MProp, RM 
          ! WRITE( KBarFile, * ) sumk( 1 : MProp ) / RM

       END IF
    END DO Bearing

    IF ( ALLOCATED( kz2 ) ) DEALLOCATE( kz2, thetaT, S )

  END SUBROUTINE Evaluate3D

  !**********************************************************************!

  SUBROUTINE SourceElement( iElement, Outside, RIn, ROut, xs, ys, MProp, M, maxM, k, PhiR, PhiS, const, kIn, PhiIn, kOut, PhiOut )

    ! Given an element number for the source                            
    ! Computes modal information

    INTEGER, INTENT( IN  ) :: M( * ), maxM                                ! number of modes
    REAL (KIND=8), INTENT( IN  ) :: xs, ys                                ! source coordinates
    COMPLEX, INTENT( IN  ) :: PhiS( maxM, * ), PhiR( maxM, * )            ! mode shapes
    COMPLEX, INTENT( IN  ) :: k(   maxM, * )
    COMPLEX, INTENT( OUT ) :: kIn( maxM ), kOut( maxM ) ! wavenumbers
    INTEGER, INTENT( OUT ) :: MProp                     ! number of modes propagating
    INTEGER, INTENT( OUT ) :: Outside                   ! side through which radial exits
    REAL (KIND=8), INTENT( OUT ) :: RIn, ROut                 ! in and out ranges for the radial
    COMPLEX, INTENT( OUT ) :: PhiIn( * ), PhiOut( * )   ! mode shapes
    COMPLEX, INTENT( OUT ) :: Const( * )                ! interpolated mode excitation coefficients
    INTEGER :: Inside, IBad, IGood1, IGood2, &
         iCorner1, iCorner2, iCorner3, iCorner4, iElement, &
         iSide, iSet1, iSet2, L, node1, node2
    REAL (KIND=8) :: R,  RV( 3 ), SV( 3 ), RVC( 3 ), SIn, SOut, x1S, y1S, xCenter, yCenter, &
         Delta, tx, ty, txC, tyC, alpha

    MProp = HUGE( MProp )

    ! Coordinates of the centroid of the source element                 

    xCenter = SUM( x( node( 1 : 3, iElement ) ) ) / 3.0
    yCenter = SUM( y( node( 1 : 3, iElement ) ) ) / 3.0 

    Side: DO iSide = 1, 3 

       node1 = node( iCorner( iSide, 1 ), iElement )
       node2 = node( iCorner( iSide, 2 ), iElement )
       iSet1 = iSet( node1 )
       iSet2 = iSet( node2 ) 
       MProp = MIN( MProp, M( iSet1 ), M( iSet2 ) ) 

       x1S = x( node1 ) - xs
       y1S = y( node1 ) - ys
       txC = x( node1 ) - xCenter
       tyC = y( node1 ) - yCenter                                 
       tx  = x( node2 ) - x( node1 )
       ty  = y( node2 ) - y( node1 ) 

       Delta = tsx * ty - tsy * tx 

       ! *** Radial parallel to side? ***                               

       IF ( Delta == 0.0 ) THEN 
          SV(  iSide ) = HUGE( SV( iSide ) )
       ELSE 
          RVC( iSide ) = ( txC * ty  - tyC * tx  ) / Delta 
          RV(  iSide ) = ( x1S * ty  - y1S * tx  ) / Delta 
          SV(  iSide ) = ( x1S * tsy - y1S * tsx ) / Delta 
       ENDIF
    END DO Side

    ! Identify two good sides and one bad side based on the intercept point

    IBad = 1 
    IF ( ABS( SV( 2 ) - 0.5 ) > ABS( SV( IBad ) - 0.5 ) ) IBad = 2                                                       
    IF ( ABS( SV( 3 ) - 0.5 ) > ABS( SV( IBad ) - 0.5 ) ) IBad = 3                                                     

    IGood1 = 1
    IGood2 = 2
    IF ( IBad == 1 ) THEN 
       IGood1 = 3 
    ELSE IF ( IBad == 2 ) THEN 
       IGood2 = 3 
    END IF

    ! The side with the lesser RVC is the inside                    

    IF ( RVC( IGood1 ) < RVC( IGood2 ) ) THEN 
       Inside  = IGood1
       Outside = IGood2 
    ELSE 
       Inside  = IGood2
       Outside = IGood1 
    ENDIF

    sIn      = SV( Inside ) 
    sIn      = MIN( MAX( sIn,  0.0D0 ), 1.0D0 ) 
    RIn      = RV( Inside ) 

    SOut     = SV( Outside ) 
    SOut     = MIN( MAX( SOut, 0.0D0 ), 1.0D0 ) 
    ROut     = RV( Outside ) 

    ! Get values of modes at Z = sz and (x, y) = intercept points   

    iCorner1 = iCorner(  Inside, 1 )
    iCorner2 = iCorner(  Inside, 2 ) 
    iCorner3 = iCorner( Outside, 1 )
    iCorner4 = iCorner( Outside, 2 ) 

    ! Interpolate to get modal values at source                     
    R = 0.0 

    IF ( RIn /= ROut ) THEN 
       alpha = ( R - RIn ) / ( ROut - RIn )
       alpha = MIN( MAX( alpha,  0.0D0 ), 1.0D0 )
    ELSE 
       alpha = 0.0 
    ENDIF

    ! following could be updated to do as a single vector operation
    Mode: DO L = 1, MProp 
       PhiIn(  L ) = CMPLX( PhiS( L, iCorner1 ) + sIn   * ( PhiS( L, iCorner2 ) - PhiS( L, iCorner1 ) ) )  
       PhiOut( L ) = CMPLX( PhiS( L, iCorner3 ) + SOut  * ( PhiS( L, iCorner4 ) - PhiS( L, iCorner3 ) ) )
       const(  L ) = CMPLX( PhiIn( L )          + alpha * ( PhiOut( L )         - PhiIn( L ) ) )
    END DO Mode

    ! Obtain values of modes at z = Rd and (x, y)=intercept points      

    CALL InterpolateModes( iElement, Inside,  sIn,  MProp, M, maxM, k, PhiR, kIn,  PhiIn  )          
    CALL InterpolateModes( iElement, Outside, sOut, MProp, M, maxM, k, PhiR, kOut, PhiOut )          

  END SUBROUTINE SourceElement

  !**********************************************************************!

  SUBROUTINE OUT( iElement, NewElement, Outside, ROut, xs, ys, MProp, M, maxM, k, PhiR, kOut, PhiOut )

    ! Given an element number and a side through which prop path enters
    ! computes mode info at exit point

    INTEGER, INTENT( IN  ) :: NewElement
    INTEGER, INTENT( IN  ) :: iElement         ! Current element number
    INTEGER, INTENT( IN  ) :: M( * ), maxM     ! number of modes, max number allowed
    REAL (KIND=8), INTENT( IN  ) :: xs, ys     ! source coordinate
    COMPLEX, INTENT( IN  ) :: PhiR( maxM, * )  ! mode shapes at each node
    COMPLEX, INTENT( IN  ) :: k(   maxM, * )   ! wavenumbers at each node
    INTEGER, INTENT( OUT ) :: Outside          ! side through which path exits
    INTEGER, INTENT( OUT ) :: MProp            ! number of propagating modes
    REAL (KIND=8), INTENT( OUT ) :: ROut             ! range at which path exits
    COMPLEX, INTENT( OUT ) :: kOut( maxM ), PhiOut( * )   ! modes at exit point
    INTEGER                :: iSide, node1, node2
    REAL (KIND=8)          :: ROutT, sOut, sT, x1S, y1S, Delta, Tx, Ty

    ! WRITE( *, * ) '   Crossing into new element = ', NewElement           

    ! If no adj elt then exit. Previous values of KOut, PhiOut retained 

    IF ( NewElement == 0 ) THEN 
       ROut = HUGE( ROut )
       RETURN 
    ENDIF

    ! *** Loop over sides to find outside ***                           

    sOut = HUGE( sOut )

    Side: DO iSide = 1, 3 

       ! Don't try an outside the same as the inside             
       IF ( AdjacentElement( iSide, NewElement ) /= iElement ) THEN 

          node1 = node( iCorner( iSide, 1 ), NewElement ) 
          node2 = node( iCorner( iSide, 2 ), NewElement ) 

          x1S = x( node1 ) - xs
          y1S = y( node1 ) - ys                                               
          Tx  = x( node2 ) - x( node1 )
          Ty  = y( node2 ) - y( node1 ) 

          Delta = tsx * Ty - tsy * Tx 

          ! *** Radial parallel to side? ***                            

          IF ( Delta == 0.0 ) THEN
             ROutT = HUGE( ROutT )
             sT    = HUGE( sT )
          ELSE 
             ROutT = ( x1S * Ty  - y1S * Tx  ) / Delta 
             sT    = ( x1S * tsy - y1S * tsx ) / Delta 
          ENDIF
          ! If intercept is in segment, store the side number       
          IF ( ABS( sT - 0.5 ) < ABS( sOut - 0.5) ) THEN 
             Outside = iSide 
             sOut    = sT 
             ROut    = ROutT
             ! write( *, * ) 'nodes', ModeFileName( Node1 ), ModeFileName( Node2 )
          ENDIF
       ENDIF
    END DO Side

    ! Obtain values of modes at intercept points                        

    CALL InterpolateModes( NewElement, Outside, sOut, MProp, M, maxM, k, PhiR, kOut, PhiOut )          

  END SUBROUTINE OUT

  !**********************************************************************!

  SUBROUTINE InterpolateModes( iElement, iSide, s, MProp, M, maxM, k, PhiR, kInt, PhiInt )          

    ! Given an element, a side, and a proportional distance along the side                                       
    ! Returns interpolated modal values

    INTEGER, INTENT( IN  ) :: iElement, iSide  ! index of the element and side
    INTEGER, INTENT( IN  ) :: M( * ), maxM     ! number of modes, max number allowed
    REAL (KIND=8), INTENT( IN  ) :: s                ! proportional distance along the side
    COMPLEX, INTENT( IN  ) :: PhiR( maxM, * )  ! mode shapes at each node
    COMPLEX, INTENT( IN  ) :: k(   maxM, * )   ! wavenumbers at each node
    INTEGER, INTENT( OUT ) :: MProp            ! number of propagating modes
    COMPLEX, INTENT( OUT ) :: PhiInt( * )      ! interpolated mode shapes
    COMPLEX, INTENT( OUT ) :: kInt( maxM )     ! interpolated wavenumbers
    INTEGER                :: iSet1, iSet2, I, node1, node2 
    REAL (KIND=8)          :: st

    node1 = node( iCorner( iSide, 1 ), iElement )
    node2 = node( iCorner( iSide, 2 ), iElement )

    iSet1 = iSet( node1 )
    iSet2 = iSet( node2 ) 

    MProp = MIN( MProp, M( iSet1 ), M( iSet2 ) ) 

    st = MIN( MAX( s, 0.0D0 ), 1.0D0 ) ! Extrapolation blocked by making sure s in [0, 1]

    ! could be updated to use a single vector operation
    Mode: DO I = 1, MProp 
       kInt(   I ) = CMPLX(    k( I, iSet1 ) + st * (    k( I, iSet2 ) -    k( I, iSet1 ) ) )
       PhiInt( I ) = CMPLX( PhiR( I, iSet1 ) + st * ( PhiR( I, iSet2 ) - PhiR( I, iSet1 ) ) )
    END DO Mode

  END SUBROUTINE InterpolateModes
END MODULE Evaluate3DMod
