function Green = kernel( B1, B2, B3, B4, rho, rk, atten, NptsAcoustic, Pos )

% Solve system for a sequence of k-values
%
global SSP NFirstAcoustic NLastAcoustic Loc H
global WS Isz WR Irz

% pre-allocate for efficiency
Nsz        = Pos.Nsz;   % # of source depths
Nrz        = Pos.Nrz;   % # of receiver ranges

Green      = zeros( Nsz, Nrz, length( rk ) );
DF         = zeros( NptsAcoustic, 1 );
EF         = zeros( NptsAcoustic, 1 );
rhoElement = zeros( NptsAcoustic, 1 );

% Tabulate z coordinates
Z( 1 ) = SSP.depth( NFirstAcoustic );
j      = 2;

for Medium = NFirstAcoustic : NLastAcoustic
    Z( j: j + SSP.N(Medium) - 1 ) = linspace( SSP.depth( Medium ) + H( Medium ), SSP.depth( Medium + 1 ), SSP.N( Medium ) );
    j = j + SSP.N( Medium );
end

[ WS, Isz ] = weight( Z, Pos.s.z ); % Compute weights for source depth interpolation
[ WR, Irz ] = weight( Z, Pos.r.z ); % Compute weights for rcvr   depth interpolation

% Assemble matrix

j       = 1;
l       = Loc( NFirstAcoustic ) + 1;
DF( 1 ) = 0.0;

for Medium = NFirstAcoustic : NLastAcoustic
    for I = 1 : SSP.N( Medium )
        rhoElement( l ) = ( rho( l ) + rho( l + 1 ) ) / 2.0;
        rhoH = rhoElement( l ) * H( Medium );
        BElement = H( Medium ) * ( ( B1( l ) + B1( l + 1 ) ) / 2.0 ) / ( 12.0 * rhoElement( l ) );
        
        DF( j     ) = DF( j ) - 1.0 / rhoH + 5.0 * BElement;
        DF( j + 1 ) =         - 1.0 / rhoH + 5.0 * BElement;
        EF( j + 1 ) =           1.0 / rhoH +       BElement;
        
        j = j + 1;
        l = l + 1;
    end
    l = l + 1;
end

for Ik = 1 : length( rk )   % Step through each point in k-space
    X = ( rk( Ik ) + 1i * atten ) ^ 2;
    Green( :, :, Ik ) = Solve( B1, B2, B3, B4, rho, NptsAcoustic, X, DF, EF, rhoElement, Pos );  % Solve for G(k)
end