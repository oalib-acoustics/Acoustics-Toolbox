function [ c2R, c2I, rho, crosst, cMin, cMax ] = initsparc( SSP, Pos )

% Initializes arrays defining difference equations

global omega Bdry
global Loc H
global WS isz WR irz Z

SSPType = Bdry.Top.Opt( 1:1 );

cMin           =  1.0E6;
cMax           = -1.0E6;
crosst         =  1.0E6;
Loc( 1 )       = 0;

% *** Loop over media ***

for Medium = 1:SSP.NMedia
    if ( Medium ~= 1 )
        Loc( Medium ) = Loc( Medium - 1 ) + SSP.N( Medium - 1 ) + 1;
    end
    N1   = SSP.N(   Medium ) + 1;
    ii   = Loc( Medium ) + 1;

    % calculate indices in user specified SSP for that medium
    if Medium == 1
        I1 = 1;
    else
        I1 = 1 + sum( SSP.Npts( 1:Medium-1) );
    end
    I2 = I1 + SSP.Npts( Medium ) - 1;
    II = I1 : I2;
    ZT  = linspace( SSP.z( I1 ), SSP.z( I2 ), N1 );
    jj = ii:ii + SSP.N( Medium );

    switch ( SSPType )
        case ( 'N' )   % n2-linear approximation to SSP
            cp        = 1./ sqrt( interp1( SSP.z( II ), 1./SSP.c(  II ).^2,  ZT ) );
            if ( any( SSP.cs( II ) ) )
               error( 'shear not allowed' )
            end
            rho( jj ) = interp1( SSP.z( II ), SSP.rho( II ),  ZT );
        case ( 'C' )   % C-LINEAR approximation to SSP
            cp        = interp1( SSP.z( II ), SSP.c(   II ),  ZT );
            rho( jj ) = interp1( SSP.z( II ), SSP.rho( II ),  ZT );
        case ( 'S' )
            cp        = interp1( SSP.z( II ), SSP.c(   II ),  ZT, 'spline' );
            rho( jj ) = interp1( SSP.z( II ), SSP.rho( II ),  ZT, 'spline' );
        case ( 'A' )
            error( '    ANALYTIC SSP option not available' )
    end

    cMinV = min( real( cp ) );
    cMin  = min( cMin, cMinV );
    cMaxV = max( real( cp ) );
    cMax  = max( cMax, cMaxV );
    
    c2R( jj ) = real( cp.^2 );
    c2I( jj ) = imag( cp.^2 ) ./ real( cp.^2) / omega ;
    
    crosst1 = min( real( H( Medium ) ./ cp ) );
    crosst  = min( crosst1, crosst );

end % next Medium

% Calculate weights for interpolation of source/receiver data

% Tabulate z coordinates
Z( 1 ) = SSP.depth( 1 );
j      = 2;
for Medium = 1 : SSP.NMedia
    Z( j: j + SSP.N(Medium) - 1 ) = linspace( SSP.depth( Medium ) + H( Medium ), SSP.depth( Medium + 1 ), SSP.N( Medium ) );
    j = j + SSP.N( Medium );
end

% Compute weights for source/rcvr depth interpolation

[ WS, isz ] = weight( Z, Pos.s.z );
[ WR, irz ] = weight( Z, Pos.r.z );
