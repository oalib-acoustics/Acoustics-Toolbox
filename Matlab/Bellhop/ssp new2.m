function [ c, gradc, crr, crz, czz, Layer ]= ssp( x, zSSPV, cSSPV, czV, TopOpt ) 

global Layer 

% ssp
% tabulates the sound speed profile and its derivatives
% also returns a vector Layer indicating the layer a depth point is in

% Layer is a vector of indices indicating the layer that each ray is in
% [zSSPV, sSSPV] contains the depth/sound speed values
% cSSPV is a vector containing the 
% all vectors are set up as column vectors ...

Npts    = length( x( :, 1 ) );
Layer   = ones( Npts, 1 );

c  = interp1( zSSPV,          cSSPV, x( :, 2 ), 'cubic', 'extrap' );
cz = interp1( zSSPV(1:end-1), czV,   x( :, 2 ), 'cubic', 'extrap' );

gradc = [ zeros(Npts, 1 ) cz ];
crr   = zeros( Npts, 1 );
crz   = zeros( Npts, 1 );
czz   = zeros( Npts, 1 );
