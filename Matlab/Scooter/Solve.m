function Green = Solve( B1, B2, B3, B4, rho, NptsAcoustic, X, DF, EF, rhoElement, Pos ) 

% Set up the linear system and solve

global SSP Bdry
global NFirstAcoustic NLastAcoustic Loc H
global WS Isz WR Irz

Nsz        = Pos.Nsz;   % # of source depths
Nrz        = Pos.Nrz;   % # of receiver ranges

Green = zeros( Nsz, Nrz );
d     = zeros( NptsAcoustic, 1 );
e     = zeros( NptsAcoustic, 1 );

% Complete assembly of matrix by adding in X ***

j      = 1;
L      = Loc( NFirstAcoustic ) + 1;
d( 1 ) = DF( 1 );

% the following should be rewritten to pre-calculate the mass matrix
% should not be recalculated for every X

for Medium = NFirstAcoustic : NLastAcoustic
    XT = -H( Medium ) * X / 12.0;
    BElement = XT ./ rhoElement( L : L + SSP.N( Medium ) - 1 );
      
    % form global matrix as sum of local matrices
    j1    = j + 1;
    jend  = j + SSP.N( Medium ) -1;
    jend1 = jend + 1;
    d( j1 : jend1 ) = DF( j1 : jend1 ) + 5 * BElement;
    d( j  : jend  ) = d(  j  : jend  ) + 5 * BElement;
    e( j1 : jend1 ) = EF( j1 : jend1 ) +     BElement;
   
    j = j + SSP.N( Medium );
    L = L + SSP.N( Medium ) + 1;
end    % next element

% Corner elt requires top impedance

BCType( 1 : 1 ) = Bdry.Top.Opt( 2 : 2 );
[ F, G, IPow ]  = bcimp( B1, B2, B3, B4, rho, X, 'TOP', Bdry.Top );

if ( G == 0.0 )
    d( 1 ) = 1.0;
    e( 2 ) = 0.0;
else
    d( 1 ) = d( 1 ) + F / G;
end

% Corner elt requires bottom impedance

BCType( 1 : 1 ) = Bdry.Bot.Opt( 1 : 1 );
[ F, G, IPow ]  = bcimp( B1, B2, B3, B4, rho, X, 'BOT', Bdry.Bot );

if ( G == 0.0 )
    d( NptsAcoustic ) = 1.0;
    e( NptsAcoustic ) = 0.0;
else
    d( NptsAcoustic ) =  d( NptsAcoustic ) - F / G;
end

[ mults, dt, et ] = factortri( NptsAcoustic, d, e );

for iS = 1 : Nsz   	% Loop over all source positions
    
    % Set up RHS in b (previously used for diagonal)
    
    b          = zeros( NptsAcoustic, 1 );
    I          = Isz( iS );
    
    elt        = Loc( NFirstAcoustic ) + I;
    rhosd      = rhoElement( iS );    % FIXME: this may be mis-indexed !!!

    b( I     ) = 2.0 * ( 1.0 - WS( iS ) ) / rhosd;
    b( I + 1 ) = 2.0 *     WS( iS )       / rhosd;
    b = complex( b );   % if you use the mex version of backsub, it requires b is complex

    % Gaussian field via complex source point
%     a  = 1000;
%     c0 = 1500;
%     freq = 50;
%     k  = 2 * pi * freq / c0;
%     %const = exp( -k * a ) / ( 1i * a );
%     %const = 1 / ( 1i * a );
%     const = 1;
%     z  = linspace( 0, 5000, NptsAcoustic );
%     x  = 100.0;
%     
%     r = sqrt( ( x - 1i * a )^2 + ( z - Pos.s.z( 1 ) ).^2 );
% 
%     pscaled = exp( k * ( 1i * r - a ) ) ./ r;
%     pscaled = pscaled / const;
%     b = pscaled;

    %save
    %figure; plot( z, abs( b ) ); drawnow; pause
    % solve the linear system
    x = backsub( NptsAcoustic, mults, dt, et, b );
    Green( iS, : ) = ( x( Irz ) + WR .* ( x( Irz + 1 ) - x( Irz ) ) ); % extract the solution at the rcvr depths
end
