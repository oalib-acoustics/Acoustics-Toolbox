function [ psi, k ] = Solve2( B1, B2, B3, B4, rho, NptsAcoustic, DF, EF, rhoElement, Pos ) 

% Set up the linear system and solve

global Bdry
global N NFirstAcoustic NLastAcoustic Loc H
global WS Isd WR Ird

Nsd = length( Pos.s.z );    % # of source depths
Nrd = length( Pos.r.z );    % # of receiver ranges

Green = zeros( Nsd, Nrd );
d     = zeros( NptsAcoustic, 1 );
e     = zeros( NptsAcoustic, 1 );
dR    = zeros( NptsAcoustic, 1 );
eR    = zeros( NptsAcoustic, 1 );
X = 0   % dummy value used in BCIMP

% Assemble global matrices from local matrices

j       = 1;
L       = Loc( NFirstAcoustic ) + 1;
d(  1 ) = DF( 1 );
dR( 1 ) = 0;

for Medium = NFirstAcoustic : NLastAcoustic
    XT = H( Medium ) / 12.0;
    BElement = XT ./ rhoElement( L : L + N( Medium ) - 1 );
      
    % form global matrix as sum of local matrices
    d( j+1 : j + N( Medium )     ) = DF( j+1 : j + N( Medium )     );
    d( j   : j + N( Medium ) - 1 ) = d(  j   : j + N( Medium ) - 1 );
    e( j+1 : j + N( Medium )     ) = EF( j+1 : j + N( Medium )     );

    dR( j+1 : j + N( Medium )     ) = 5 * BElement( 1: N( Medium ) );
    dR( j   : j + N( Medium ) - 1 ) = dR(  j   : j + N( Medium ) - 1 ) + 5 * BElement( 1: N( Medium ) );
    eR( j+1 : j + N( Medium )     ) = BElement( 1: N( Medium ) );
   
    j = j + N( Medium );
    L = L + N( Medium ) + 1;
end    % next element

% Corner elt requires top impedance

BCType(1:1) = Bdry.Top.Opt(2:2);
[ F, G, IPow ] = bcimp( B1, B2, B3, B4, rho, X, 'TOP', Bdry.Top );

if ( G == 0.0 )
    d( 1 ) = 1.0;
    e( 2 ) = 0.0;
else
    d( 1 ) = d( 1 ) + F / G;
end

% Corner elt requires bottom impedance

BCType(1:1) = Bdry.Bot.Opt(1:1);
[ F, G, IPow ] = bcimp( B1, B2, B3, B4, rho, X, 'BOT', Bdry.Bot );

if ( G == 0.0 )
    d( NptsAcoustic ) = 1.0;
    e( NptsAcoustic ) = 0.0;
else
    d( NptsAcoustic ) =  d( NptsAcoustic ) - F / G;
end

%[ mults, dt, et ] = factortri( NptsAcoustic, d, e );

nz = NptsAcoustic;

D = sparse( 1:nz, 1:nz  , d, nz, nz, nz   );
E = sparse( 2:nz, 1:nz-1, e( 2 : end ),     nz, nz, nz-1 );
A = D + E + E';

D = sparse( 1:nz, 1:nz  , dR, nz, nz, nz   );
E = sparse( 2:nz, 1:nz-1, eR( 2 : end ),     nz, nz, nz-1 );
B = D + E + E';
%[ D; E ]

%[ psi, x ] = eigs( A, B, nz );
[ psi, x ] = eig( full( A ), full( B ) );
save
k = sqrt( x );

%Green( IS, : ) = ( x( Ird ) + WR .* ( x( Ird + 1 ) - x( Ird ) ) ); % extract the solution at the rcvr depths
