MODULE EvaluateGBMod
  IMPLICIT NONE
  REAL,    PARAMETER :: pi = 3.141592, DegRad = pi / 180.0, RadDeg = 180.0 / pi, c0 = 1500   ! reference sound speed, should be speed at the source depth
  COMPLEX, PARAMETER :: i = ( 0.0, 1.0 )

CONTAINS
  SUBROUTINE EvaluateGB( k, phi, phiR, M, maxM, IelementSource, xS, yS, theta, Ntheta, Rmin, Rmax, NR, MLimit, Option, P, freq )

    ! Computes 3-D pressure field using adiabatic mode theory.
    ! Uses gaussian beam tracing for horizontal refraction.           
    ! Beam curvature change at interfaces still needed.          

    ! Phase may be off by e( i pi / 4 ) but TL is correct.

    ! Rays here should be complex; they are projected here onto the real plane
    ! This happens wherever we use DBLE( cA ) instead of the complex cA

    USE ElementMod
    USE FatalError

    INTEGER, PARAMETER    :: FLPFile = 5, RAYFile = 21
    INTEGER, INTENT( IN ) :: M( * ), MLimit, maxM            ! number of modes, limit on number of modes to propagate
    INTEGER, INTENT( IN ) :: Nr, Ntheta, IElementSource      ! number of receiver ranges and bearing lines
    REAL (KIND=8), INTENT( IN ) :: Rmin, Rmax                ! minimum and maximum receiver ranges in m
    REAL (KIND=8), INTENT( IN ) :: xs, ys                    ! source coordinate
    REAL (KIND=8), INTENT( IN ) :: freq                      ! source frequency
    REAL (KIND=8), INTENT( IN ) :: theta( Ntheta )           ! bearing angles for receivers
    COMPLEX, INTENT( IN ) :: k( maxM, * )                    ! wavenumbers
    COMPLEX, INTENT( IN ) :: PhiR( maxM, * ), Phi( maxM, * ) ! source/receiver mode shapes
    COMPLEX, INTENT( OUT) :: P( Ntheta, Nr )                 ! pressure field
    LOGICAL            :: EXITED
    INTEGER            :: Ialpha, IElement, Iset1, Iset2, Iset3, Iside, Istep, KMAHA, KMAHB, mode, Mprop, Nsteps, Nalpha, N2
    REAL    (KIND=8), ALLOCATABLE :: RadVec( :, : )
    REAL    (KIND=8), ALLOCATABLE :: xV( : ), yV( : )
    REAL    (KIND=8)   :: alpha, Dalpha, A1, A2, A3, D12, D13, D23, DB1, DB2, DB3, delta, &
         EpsilonMult, qAT, qBT, tsx, tsy, xa, xb, ya, yb, x1, y1, x2, y2, x3, y3, &
         xia, xib, etaa, etab, alpha1, alpha2, step, stepG, Hwidth
    COMPLEX (KIND=8)   :: const, phiA, phiB, pA, pB, qA, qB, tauA, tauB, cA, cB, eps, EpsOpt, cx, cy
    CHARACTER          :: Option*( * )
    CHARACTER (LEN=50) :: FileRoot = 'RayFile', Title = 'KRAKEN3D'

    ! following should only be done on first call !!!
    READ( FLPFile, * ) alpha1, alpha2, Nalpha 
    READ( FLPFile, * ) stepG, Nsteps 
    READ( FLPFile, * ) EpsilonMult

    ! automatic calculation of number of beams
    IF ( Nalpha == 0 ) THEN
       Nalpha = MAX( INT( 0.3 * Rmax * freq / c0 * ( alpha2 - alpha1 ) / 360 ), 300 )
    END IF

    WRITE( PRTFile, * ) 'alpha1, alpha2    = ', alpha1, alpha2
    WRITE( PRTFile, * ) 'Number of beams   = ', Nalpha
    WRITE( PRTFile, * ) 'stepsize          = ', stepG
    WRITE( PRTFile, * ) 'Number of steps   = ', Nsteps
    WRITE( PRTFile, * ) 'Epsilon multplier = ', EpsilonMult

    alpha1 = DegRad * alpha1 
    alpha2 = DegRad * alpha2 
    Dalpha = ( alpha2 - alpha1 ) / ( Nalpha - 1 )

    ! open ray file, if option selected
    IF ( Option( 6 : 6 ) == 'R' ) THEN
       OPEN ( FILE = TRIM( FileRoot ) // '.ray', UNIT = RAYFile, FORM = 'FORMATTED' )
       WRITE( RAYFile, * ) '''', Title( 1 : 50 ), ''''
       WRITE( RAYFile, * ) freq
       WRITE( RAYFile, * ) 1, 1, 1  ! Pos%NSx, Pos%NSy, Pos%NSz
       WRITE( RAYFile, * ) Nalpha, MLimit
       WRITE( RAYFile, * ) 0        ! Bdry%Top%HS%Depth
       WRITE( RAYFile, * ) 200      ! Bdry%Bot%HS%Depth
       WRITE( RAYFile, * ) '''xyz'''
    END IF

    ! Radial vector for each bearing line
    ALLOCATE( RadVec( 2, Ntheta ), xV( Nsteps ), yV( Nsteps ) )
    RadVec( 1, : )          = COS( DegRad * theta( 1 : Ntheta ) )
    RadVec( 2, : )          = SIN( DegRad * theta( 1 : Ntheta ) )

    P( 1 : Ntheta, 1 : NR ) = 0.0

    ! Loop over beam takeoff angle
    Angle: DO Ialpha = 1, Nalpha
       !IF ( Ialpha /= 110 ) CYCLE Angle   ! run a single beam
       alpha = alpha1 + ( Ialpha - 1 ) * Dalpha 
       tsx   = COS( alpha ) 
       tsy   = SIN( alpha ) 
       WRITE( *, * ) 'Ialpha, tsx, tsy', Ialpha, tsx, tsy 

       ! Loop over modes
       ModeLoop: DO mode = 1, MLimit 

          !IF ( mode /= 4 ) CYCLE ModeLoop   ! run a single mode
          ! Get mode excitation coef and initialize
          IElement = IelementSource 
          CALL NewElement( IElement, k, mode, M, maxM, Iset1, Iset2, Iset3, &
               x1, y1, x2, y2, x3, y3, D12, D13, D23, delta, Cx, Cy, Mprop )
          IF ( mode > Mprop ) CYCLE ANGLE

          ! Evaluate modes at source depth                          
          DB1 = xS * y1 - yS * x1 
          DB2 = xS * y2 - yS * x2 
          DB3 = xS * y3 - yS * x3 

          ! Compute areas of triangles                              
          A1   = (  D23 - DB3 + DB2 ) / delta 
          A2   = (  DB3 - D13 - DB1 ) / delta 
          A3   = ( -DB2 + DB1 + D12 ) / delta 
          cA   = A1 / k( mode, Iset1 ) + A2 / k( mode, Iset2 ) + A3 / k( mode, Iset3 )
          phiA = A1 * phiR( mode, 1  ) + A2 * phiR( mode, 2  ) + A3 * phiR( mode, 3  )

          SELECT CASE ( Option( 5 : 5 ) )
          CASE( 'F' )   ! Space filling in far field
             Hwidth = 2.0 / ( ( 1.0 / DBLE( cA ) ) * Dalpha ) 
             EpsOpt = 0.5 * Hwidth ** 2 
          CASE( 'M' )   ! Minimum width at Rmax
             Hwidth = SQRT( 2.0 * DBLE( cA ) * Rmax ) 
             EpsOpt = 0.5 * Hwidth ** 2
          CASE DEFAULT
             EpsOpt = 0.0
             CALL ERROUT( 'EvaluateBG', 'Unknown option fpr beam type' ) 
          END SELECT

          eps   = EpsilonMult * i * EpsOpt 
          const = phiA * SQRT( eps / cA ) * Dalpha 
          xA    = xS 
          yA    = yS 
          xiA   = tsx / DBLE( cA )
          etaA  = tsy / DBLE( cA ) 
          pA    = 1.0 
          qA    = eps 
          tauA  = 0.0 
          KMAHA = 1 

          ! Evaluate modes at rcvr depth                            
          phiA = A1 * phi( mode, Iset1 ) + A2 * phi( mode, Iset2 ) + A3 * phi( mode, Iset3 )        

          ! March forward in range
          EXITED = .FALSE. 
          step = stepG 

          Stepping: DO Istep = 1, Nsteps 
             xB   = xA   + step * DBLE( cA ) * xiA 
             yB   = yA   + step * DBLE( cA ) * etaA 
             xiB  = xiA  - step * DBLE( cx / ( cA * cA ) )
             etaB = etaA - step * DBLE( cy / ( cA * cA ) ) 
             pB   = pA 
             !qB   = qA   + step * cA * pA
             qB   = qA   + step * DBLE( cA ) * pA
             tauB = tauA + step / cA 

             step = stepG   ! Update step size

             IF ( Option( 6 : 6 ) == 'R' ) THEN 
                xV( Istep ) = xA 
                yV( Istep ) = yA 
             END IF

             ! Compute KMAH index                                   
             KMAHB = KMAHA 
             IF ( REAL( qB ) < 0.0 ) THEN 
                qAT = AIMAG( qA ) 
                qBT = AIMAG( qB ) 
                IF ( ( qAT < 0.0 .AND. qBT >= 0.0 ) .OR. &
                     ( qAT > 0.0 .AND. qBT <= 0.0 ) ) KMAHB = -KMAHA
             END IF

             ! Mode interpolation coefs                             
             IF ( .NOT. EXITED ) THEN 
                DB1 = xB * y1 - yB * x1 
                DB2 = xB * y2 - yB * x2 
                DB3 = xB * y3 - yB * x3 

                ! Compute areas of triangles                        
                A1 = (  D23 - DB3 + DB2 ) / delta 
                A2 = (  DB3 - D13 - DB1 ) / delta 
                A3 = ( -DB2 + DB1 + D12 ) / delta 

                ! Crossing into new element?                            
                ! EltSearch: DO WHILE ( ABS( A1 ) + ABS( A2 ) + ABS( A3 ) > 1.0 + 1e4 * EPSILON( A1 ) )
                EltSearch: DO WHILE ( A1 < 0.0 .OR. A2 < 0.0 .OR. A3 < 0.0 )

                   ! Identify the side through which exitting
                   Iside = 1   ! just in case none of the following happens
                   IF ( A1 < 0.0 ) Iside = 2 
                   IF ( A2 < 0.0 ) Iside = 3 
                   IF ( A3 < 0.0 ) Iside = 1
                   IElement = AdjacentElement( Iside, IElement ) 

                   ! Normal elt transition or exit into free space? 
                   IF ( IElement == 0 ) THEN 
                      EXITED = .TRUE. 
                      cx = 0.0 
                      cy = 0.0
                      EXIT EltSearch
                   ELSE 
                      CALL NewElement( IElement, k, mode, M, maxM, Iset1, Iset2, Iset3, &
                           x1, y1, x2, y2, x3, y3, D12, D13, D23, delta, cx, cy, Mprop )
                      IF ( Mprop < mode ) EXIT EltSearch ! If mode cuts off, skip to next mode
                      DB1 = xB * y1 - yB * x1 
                      DB2 = xB * y2 - yB * x2 
                      DB3 = xB * y3 - yB * x3 

                      ! Compute areas of triangles                        
                      A1 = (  D23 - DB3 + DB2 ) / delta
                      A2 = (  DB3 - D13 - DB1 ) / delta
                      A3 = ( -DB2 + DB1 + D12 ) / delta
                   END IF
                END DO EltSearch

                IF ( Mprop < mode .OR. IElement == 0 ) EXIT Stepping   ! If mode cuts off, skip to next mode

                ! Evaluate modes at the rcvr depth                  
                cB   = A1 /   k( mode, Iset1 ) + A2 /   k( mode, Iset2 ) + A3 /   k( mode, Iset3 )          
                phiB = A1 * phi( mode, Iset1 ) + A2 * phi( mode, Iset2 ) + A3 * phi( mode, Iset3 )
             END IF

             ! Compute beam influence
             CALL InfluenceR( xA - xS, yA - yS, xiA, etaA, pA, qA, tauA, cA, KMAHA, phiA, &
                  xB - xS, yB - yS, xiB, etaB, pB, qB, tauB, cB, KMAHB, phiB, RadVec, Ntheta, Rmin, Rmax, NR, const, P )

             xA    = xB
             yA    = yB
             xiA   = xiB
             etaA  = etaB
             pA    = pB
             qA    = qB
             tauA  = tauB
             cA    = cB
             KMAHA = KMAHB
             phiA  = phiB

          END DO Stepping

          ! Optionally dump rays to disk
          IF ( Option( 6 : 6 ) == 'R' ) THEN 
             N2 = Istep - 1 
             WRITE( RAYFile, * ) alpha / DegRad
             WRITE( RAYFile, * ) N2, mode, 0
             DO istep = 1, N2
                WRITE( RAYFile, * ) SNGL( xV( istep ) ), SNGL( yV( istep ) ), 0.0
             END DO
          END IF

       END DO ModeLoop

    END DO Angle

    RETURN 
  END SUBROUTINE EvaluateGB

  !**********************************************************************

  SUBROUTINE NewElement( IElement, k, mode, M, maxM, Iset1, Iset2, Iset3, &
       x1, y1, x2, y2, x3, y3, D12, D13, D23, delta, cx, cy, Mprop )

    ! Given elt number, returns info which is constant in elt           

    USE ElementMod

    INTEGER,          INTENT( IN    ) :: M( * ), mode, maxM, IElement
    INTEGER,          INTENT( INOUT ) :: Iset1, Iset2, Iset3, Mprop
    REAL (KIND=8),    INTENT( INOUT ) :: delta, x1, y1, x2, y2, x3, y3
    REAL (KIND=8),    INTENT( OUT   ) :: D12, D13, D23
    COMPLEX,          INTENT( IN    ) :: k( maxM, * )
    COMPLEX (KIND=8), INTENT( OUT   ) :: cx, cy
    INTEGER                           :: node1, node2, node3
    REAL (KIND=8)                     :: A1x, A2x, A3x, A1y, A2y, A3y

    node1 = node( 1, IElement ) 
    node2 = node( 2, IElement ) 
    node3 = node( 3, IElement ) 

    Iset1 = Iset( node1 ) 
    Iset2 = Iset( node2 ) 
    Iset3 = Iset( node3 ) 

    ! If mode cuts off, return to origin and do next mode           
    Mprop = MIN( M( Iset1 ), M( Iset2 ), M( Iset3 ) ) 
    IF ( mode > Mprop ) RETURN 

    x1 = x( node1 ) 
    y1 = y( node1 ) 
    x2 = x( node2 ) 
    y2 = y( node2 ) 
    x3 = x( node3 ) 
    y3 = y( node3 ) 

    D12   = x1 * y2 - y1 * x2 
    D13   = x1 * y3 - y1 * x3 
    D23   = x2 * y3 - y2 * x3 
    delta = D23 - D13 + D12 

    ! Gradient                                                      
    A1x = -y3 + y2 
    A2x =  y3 - y1 
    A3x = -y2 + y1 
    A1y =  x3 - x2 
    A2y = -x3 + x1 
    A3y =  x2 - x1 

    cx = ( A1x / k( mode, Iset1 ) + A2x / k( mode, Iset2 ) + A3x / k( mode, Iset3 ) ) / delta
    cy = ( A1y / k( mode, Iset1 ) + A2y / k( mode, Iset2 ) + A3y / k( mode, Iset3 ) ) / delta

    RETURN 
  END SUBROUTINE NewElement

  !**********************************************************************

  SUBROUTINE InfluenceR( xA, yA, xiA, etaA, pA, qA, tauA, cA, KMAHA, phiA, &
       xB, yB, xiB, etaB, pB, qB, tauB, cB, KMAHB, phiB, &
       RadVec, Ntheta, Rmin, Rmax, NR, const, P )

    ! Computes contribution to receivers assuming beam cannot           
    ! contribute to a radial in an incoming sense                       

    REAL,             PARAMETER       :: BeamWindow = 5
    INTEGER,          INTENT( IN )    :: KMAHA, KMAHB, Ntheta, NR
    REAL    (KIND=8), INTENT( IN )    :: RadVec( 2, * )
    REAL    (KIND=8), INTENT( IN )    :: Rmin, Rmax
    REAL    (KIND=8), INTENT( IN )    :: xA, yA, xiA, etaA, xB, yB, xiB, etaB
    COMPLEX (KIND=8), INTENT( IN )    :: const, phiA, phiB, pA, pB, qA, qB, tauA, tauB, cA, cB
    COMPLEX,          INTENT( INOUT ) :: P( Ntheta, * )
    INTEGER            :: KMAH, Itheta, ir1, ir2, ir
    REAL    (KIND=8)   :: nA, nB, Nsquared, R, rA, rB, deltar, deltaA, deltaB, W, qAT, qBT 
    COMPLEX (KIND=8)   :: phiMid, pMid, qMid, tauMid, cMid, contrib

    deltaR = ( Rmax - Rmin ) / ( NR - 1 ) 

    ! Loop over radials of receiver line and compute contribution

    Bearing: DO Itheta = 1, Ntheta

       ! Compute intercept range, ra, & index preceding rcvr        
       deltaA = RadVec( 1, Itheta ) * xiA + RadVec( 2, Itheta ) * etaA 
       deltaB = RadVec( 1, Itheta ) * xiB + RadVec( 2, Itheta ) * etaB

       IF ( ABS( deltaA ) < TINY( deltaA ) ) CYCLE
       IF ( ABS( deltaB ) < TINY( deltaB ) ) CYCLE

       rA  = ( yA * etaA + xA * xiA ) / deltaA
       rB  = ( yB * etaB + xB * xiB ) / deltaB

       ir1 = MAX( MIN( INT( ( rA - Rmin ) / deltaR ) + 1, Nr ), 0 )
       ir2 = MAX( MIN( INT( ( rB - Rmin ) / deltaR ) + 1, Nr ), 1 )

       ! If a receiver is bracketted, compute influence
       IF ( ir2 > ir1 .AND. deltaA * deltaB > 0 ) THEN

          ! Normal distance                                         
          nA = ( xA * RadVec( 2, Itheta ) - yA * RadVec( 1, Itheta ) ) / ( DBLE( cA ) * deltaA )     
          nB = ( xB * RadVec( 2, Itheta ) - yB * RadVec( 1, Itheta ) ) / ( DBLE( cB ) * deltaB )     

          Range: DO ir = ir1 + 1, ir2 
             r        = Rmin + ( ir - 1 ) * deltaR 
             W        = ( r - rA ) / ( rB - rA ) 
             pMid     =   pA + W * ( pB - pA ) 
             qMid     =   qA + W * ( qB - qA ) 
             nsquared = ( nA + W * ( nB - nA ) ) ** 2

             ! Within beam window?                               
             IF ( -0.5 * AIMAG( pMid / qMid ) * Nsquared < BeamWindow ) THEN
                cMid   = cA   + W * ( cB   - cA   ) 
                tauMid = tauA + W * ( tauB - tauA ) 
                phiMid = phiA + W * ( phiB - phiA )

                ! Compute KMAH index                                
                KMAH = KMAHA 
                IF ( REAL( qMid ) < 0.0 ) THEN 
                   qAT = AIMAG( qA ) 
                   qBT = AIMAG( qMid ) 
                   IF ( ( qAT < 0.0 .AND. qBT >= 0.0 ) .OR. ( qAT > 0.0 .AND. qBT <= 0.0 ) ) KMAH = -KMAH                
                END IF

                contrib = const * phiMid * SQRT( cMid / qMid ) * EXP( -i * ( tauMid + 0.5 * pMid / qMid * nsquared ) )           
                IF ( KMAH < 0 ) contrib = -contrib
                P( Itheta, ir ) = P( Itheta, ir ) + CMPLX( contrib ) 
             END IF
          END DO Range
       END IF

    END DO Bearing

    RETURN 
  END SUBROUTINE InfluenceR
END MODULE EvaluateGBMod
