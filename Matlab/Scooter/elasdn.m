function [ YVout, IPowout ] = elasdn( B1, B2, B3, B4, rho, X, YV, IPow, Medium )

%     Propagates through an elastic layer using
%     compound matrix formulation
global SSP H Loc

Roof  = 1.0E5;
Floor = 1.0E-5;
IPowR = 5;
IPowF = -5;

% Euler's method for first step

TwoX   = 2.0 * X;
TwoH   = 2.0 * H( Medium );
FourHX = 4.0 * H( Medium ) * X;
j = Loc( Medium ) + 1;

XB3 = X * B3( j ) - rho( j );

ZV(1) = YV(1) + 0.5*(   B1( j ) * YV( 4 ) - B2( j ) * YV( 5 ) );
ZV(2) = YV(2) + 0.5*( -rho( j ) * YV( 4 ) -     XB3 * YV( 5 ) );
ZV(3) = YV(3) + 0.5*(      TwoH * YV( 4 ) + B4( j ) * YV( 5 ) );
ZV(4) = YV(4) + 0.5*(    XB3 * YV(1) + B2(j) * YV(2) -TwoX * B4(j) * YV(3) );
ZV(5) = YV(5) + 0.5*( rho(j) * YV(1) - B1(j) * YV(2) -      FourHX * YV(3) );

%     Modified midpoint method

for I = 1 : SSP.N( Medium )
    j = j+1;

    XV = YV   ;   YV = ZV;

    XB3 = X * B3( j ) - rho( j );

    ZV(1) = XV(1) + (   B1( j ) * YV( 4 ) - B2( j ) * YV( 5 ) );
    ZV(2) = XV(2) + ( -rho( j ) * YV( 4 ) -     XB3 * YV( 5 ) );
    ZV(3) = XV(3) + (      TwoH * YV( 4 ) + B4( j ) * YV( 5 ) );
    ZV(4) = XV(4) + (    XB3 * YV(1) + B2(j) * YV(2) -TwoX * B4(j) * YV(3) );
    ZV(5) = XV(5) + ( rho(j) * YV(1) - B1(j) * YV(2) -      FourHX * YV(3) );

    % Scale if necessary

    if ( I ~= SSP.N( Medium ) )

        if     ( abs( real( ZV( 2 ) ) ) < Floor )
            ZV = Roof  * ZV   ;   YV = Roof * YV;
            IPow = IPow - IPowR;
        elseif ( abs( real( ZV( 2 ) ) ) > Roof )
            ZV = Floor * ZV   ;   YV = Floor * YV;
            IPow = IPow - IPowF;
        end
    end

end

%     Apply the standard filter at the terminal point

YVout   = ( XV + 2.0 * YV + ZV ) / 4.0;
IPowout = IPow;