function [ B1, B2, B3, B4, rho, Mater ] = init_elmatrix( NptsAll )

% Initializes arrays defining the finite-difference equations

global Bdry SSP omega
global NFirstAcoustic NLastAcoustic Loc h

SSPType = Bdry.Top.Opt( 1 : 1 );

% pre-allocate space for efficiency
B1  = zeros( NptsAll, 1 );
B2  = zeros( NptsAll, 1 );
B3  = zeros( NptsAll, 1 );
B4  = zeros( NptsAll, 1 );
rho = zeros( NptsAll, 1 );
Mater( SSP.NMedia, : ) = 'ACOUSTIC';

cMin           = 1.0E6;
NFirstAcoustic = 0;
Loc( 1 )       = 0;
omega2 = omega^2;

% *** Loop over media ***

for Medium = 1 : SSP.NMedia
    if ( Medium ~= 1 )
        Loc( Medium ) = Loc( Medium - 1 ) + SSP.N( Medium - 1 ) + 1;
    end
    N1   = SSP.N( Medium ) + 1;
    ii   = Loc(   Medium ) + 1;

    % calculate indices in user specified SSP for that medium
    if Medium == 1
        I1 = 1;
    else
        I1 = 1 + sum( SSP.Npts( 1 : Medium - 1 ) );
    end
    I2 = I1 + SSP.Npts( Medium ) - 1;
    II = I1 : I2;
    ZT = linspace( SSP.z( I1 ), SSP.z( I2 ), N1 );
    jj = ii : ii + SSP.N( Medium );

    switch ( SSPType )
        case ( 'N' )   % n2-linear approximation to SSP
            cp        = 1./ sqrt( interp1( SSP.z( II ), 1./SSP.c(  II ).^2,  ZT ) );
            if ( any( SSP.cs( II ) ) )
                cs    = 1./ sqrt( interp1( SSP.z( II ), 1./SSP.cs( II ).^2,  ZT ) );
            else
                cs = zeros( size( SSP.cs ) );
            end
            rho( jj ) = interp1( SSP.z( II ), SSP.rho(   II ),  ZT );
        case ( 'C' )   % C-LINEAR approximation to SSP
            cp        = interp1( SSP.z( II ), SSP.c(   II ),  ZT );
            cs        = interp1( SSP.z( II ), SSP.cs(  II ),  ZT );
            rho( jj ) = interp1( SSP.z( II ), SSP.rho( II ),  ZT );
        case ( 'P' )   % monotone PCHIP ACS (almost a cubic spline)
            cp        = pchip_acs( SSP.z( II ), SSP.c(   II ),  ZT );
            cs        = pchip_acs( SSP.z( II ), SSP.cs(  II ),  ZT );
            rho( jj ) = pchip_acs( SSP.z( II ), SSP.rho( II ),  ZT );
        case ( 'S' )
            cp        = interp1( SSP.z( II ), SSP.c(   II ),  ZT, 'spline' );
            cs        = interp1( SSP.z( II ), SSP.cs(  II ),  ZT, 'spline' );
            rho( jj ) = interp1( SSP.z( II ), SSP.rho( II ),  ZT, 'spline' );
        case ( 'A' )
            error( '    ANALYTIC SSP option not available' )
    end

    if ( ~any( cs ) )   % *** Case of an acoustic Medium ***
        Mater( Medium, : )  =  'ACOUSTIC';
        
        if ( NFirstAcoustic == 0 )
            NFirstAcoustic = Medium;
        end
        NLastAcoustic  = Medium;

        cMinV = min( real( cp ) );
        cMin  = min( cMin, cMinV );
        B1( jj ) = -2 + h( Medium )^2 * omega2 ./ cp.^2;

    else                % *** Case of an elastic medium ***
        Mater( Medium, : ) = 'ELASTIC ';
        %Twoh = 2.0 * h( Medium );

        cMinV = min( real( cs ) );
        cMin  = min( cMin, cMinV );

        cp2 = cp.^2;
        cs2 = cs.^2;

        lambda = rho( jj )' .* ( cp2 - 2.0 * cs2 );
        mu     = rho( jj )' .* cs2;

        B1(  jj ) = 1.0 ./ mu;
        B2(  jj ) = 1.0 ./ ( lambda + 2.0 * mu );
        B3(  jj ) = 4.0 .* mu .* ( lambda + mu ) ./ ( lambda + 2.0 * mu );
        B4(  jj ) = lambda ./ ( lambda + 2.0 * mu );
        rho( jj ) = omega2 * rho( jj );

        %B1(  jj ) = Twoh ./ ( rho( jj )' .* cs2 );
        %B2(  jj ) = Twoh ./ ( rho( jj )' .* cp2 );
        %B3(  jj ) = 4.0 * Twoh * rho( jj )' .* cs2 .* ( cp2 - cs2 ) ./ cp2;
        %B4(  jj ) = Twoh * ( cp2 - 2.0 .* cs2 ) ./ cp2;
        %rho( jj ) = Twoh * omega2 * rho( jj )';
    end
end % next Medium
