function [ c, gradc, crr, crz, czz, Layer ]= ssp( x, ~, ~, ~, ~ ) 

% Range-dependent Munk profile given as an analytic function

global Layer 

% ssp
% tabulates the sound speed profile

Npts = length( x( :, 1 ) );
Layer = ones( Npts, 1 );

r = x( :, 1 );
z = x( :, 2 );

c0   = 1500.0;
DcoefDr = 0.03 / 100000;
coef = 0.00737 + DcoefDr * r;
coef( 1 )

zt   = 2.0 * ( z - 1300.0 ) / 1300.0;  % used to be called x
DxDz = 2.0 / 1300.0;
c    = c0 * ( 1.0 + coef .* ( zt - 1.0 + exp( -zt ) ) );
cz   = c0 * 0.00737 * ( 1.0 - exp( -zt ) ) .* DxDz;
czz  = c0 * 0.00737 * exp( -zt ) .* DxDz ^2;

cr = c0 * DcoefDr * ( zt - 1 + exp( -zt ) );
%c = c';
gradc = [ cr cz ];

crr = zeros( Npts, 1 ) ;
crz = c0 * DcoefDr * ( 1.0 - exp( -zt ) ) * DxDz;
     