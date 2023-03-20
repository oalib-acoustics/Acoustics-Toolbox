function [ c, gradc, crr, crz, czz, Layer ]= Munk( x, ~, ~, ~, ~ ) 

% Munk profile given an an analytic function

global Layer 

% ssp
% tabulates the sound speed profile

Npts  = length( x( :, 1 ) );
Layer = ones( Npts, 1 );

z = x( :, 2 );

c0   = 1500.0;
zt   = 2.0 * ( z - 1300.0 ) / 1300.0;  % used to be called x
DxDz = 2.0 / 1300.0;
c    = c0 * ( 1.0 + 0.00737 * ( zt - 1.0 + exp( -zt ) ) );
cz   = c0 * 0.00737 * ( 1.0 - exp( -zt ) ) .* DxDz;
czz  = c0 * 0.00737 * exp( -zt ) .* DxDz ^2;

%c = c';
gradc = [ zeros(Npts, 1 ) cz ];

crr = zeros( Npts, 1 ) ;
crz = zeros( Npts, 1 );
