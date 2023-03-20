function reflect( I, BeamType, BC, cHS, rhoHS, BotTop, TBdry, NBdry, kappa, SSP )
% Reflect a ray/beam off a boundary

global ray omega Layer

RADDEG = 180 / pi;   % used to convert radians to degrees

Tg = dot( ray( I ).Tray', TBdry' );  % component of ray tangent, along boundary
Th = dot( ray( I ).Tray', NBdry' );  % component of ray tangent, normal to boundary

% *** calculate the change in curvature ***
% Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).

[ c, c_imag, gradc, ~, ~, ~, Layer ] = ssp( ray( I ).x( : ), SSP, Layer );

% incident unit ray tangent and normal
rayt = ray( I ).Tray;
rayn = [ -rayt( 2 ) rayt( 1 ) ];

%cnjump = 2 * gradc( 2 ) .* ray( I ).Tray( 1 );
%csjump = 2 * gradc( 2 ) .* ray( I ).Tray( 2 );   % assumes gradc( 1 ) = cr = 0
cnjump = 2 * gradc * rayn.';
csjump = 2 * gradc * rayt.';

RN = 2 * kappa / c^2 / Th;    % boundary curvature correction

if ( strcmp( BotTop, 'TOP' ) )
    cnjump = -cnjump;    % flip sign for top reflection
    RN     = -RN;
end

RM = Tg ./ Th;
RN = RN + RM' .* ( 2 * cnjump - RM' .* csjump ) ./ c;

switch ( BeamType(2:2) )
    case ( 'D' )
        RN = 2.0 * RN;
    case ( 'Z' )
        RN = 0.0;
end

ray( I ).Tray( 1 ) =  ray( I ).Tray( 1 ) - 2.0 * ( Th .* NBdry( 1 )' )';
ray( I ).Tray( 2 ) =  ray( I ).Tray( 2 ) - 2.0 * ( Th .* NBdry( 2 )' )';

ray( I ).p( 1 ) = ray( I ).p( 1 ) + ray( I ).q( 1 ) * RN;
ray( I ).p( 2 ) = ray( I ).p( 2 ) + ray( I ).q( 2 ) * RN;

% *** account for phase change ***

switch ( BC )
    case ( 'R' )                 % rigid
    case ( 'V' )                 % vacuum
        ray( I ).Rfa = -ray( I ).Rfa;
    case ( 'F' )                 % file
        theInt = RADDEG * abs( atan2( Th, Tg ) );   % angle of incidence (relative to normal to bathymetry)
        %[ theta, RefC, phi ] = RefCO( theInt, rInt, phiInt, Npts  )
        ray( I ).Rfa = ray( I ).Rfa * rInt * exp( 1i * phiInt );
    case ( 'A' )                 % half-space
        GK       = omega * Tg;   % wavenumber in direction parallel to bathymetry
        gamma1SQ = ( omega / c'  ).^ 2 - GK^ 2;
        gamma2SQ = ( omega / cHS ) ^ 2 - GK^ 2;
        gamma1   = sqrt( -gamma1SQ );
        gamma2   = sqrt( -gamma2SQ );
        Refl = ( rhoHS * gamma1 - gamma2 ) / ( rhoHS * gamma1 + gamma2 );
        % above should really be:
        % Refl = ( rhoHs * gamma1 - rho * gamma2 ) / ( rhoHS * gamma1 + rho * gamma2 )
        % However, BELLHOP assumes rho=1 in water
        % Need to modify this as was done in the Fortran version to have
        % ssp return rho
        ray( I ).Rfa = Refl.' .* ray( I ).Rfa;
end
