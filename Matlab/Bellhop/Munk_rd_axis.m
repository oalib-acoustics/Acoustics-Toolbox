function [ c, gradc, crr, crz, czz, Layer ]= ssp( x, zSSPV, cSSPV, czV, TopOpt ) 

global Layer 

% ssp
% tabulates the sound speed profile

Npts = length( x( :, 1 ) );
Layer = ones( Npts, 1 );

r = x( :, 1 );
z = x( :, 2 );

c0   = 1500.0;
coef = 0.00737;

DsDr = -500/100000;
sspaxis = 1300.0 + DsDr * r;

xt   = 2.0 * ( z - sspaxis ) ./ sspaxis;
DxDz = 2.0 ./ sspaxis;
DxDr = -2.0 * z ./ sspaxis.^2 .* DsDr;
DxDrDz = -2 ./ sspaxis^2 .* DsDr;
DxDrDr =  4 * z ./ sspaxis^3 .* DsDr^2;
  
c    = c0 * ( 1.0 + coef * ( xt - 1.0 + exp( -xt ) ) );

cz   = c0 * 0.00737 * ( 1.0 - exp( -xt ) ) .* DxDz;
czz  = c0 * 0.00737 * exp( -xt ) .* DxDz .^2;

cr = c0 * coef * DxDr .* ( 1 - exp( -xt ) );
%c = c';
gradc = [ cr cz ];

crz = c0 * coef * ( exp( -xt ) .* DxDz .* DxDr + ( 1 - exp( -xt ) ) .* DxDrDz );
crr = c0 * coef * ( exp( -xt ) .* DxDr .^2     + ( 1 - exp( -xt ) ) .* DxDrDr );
