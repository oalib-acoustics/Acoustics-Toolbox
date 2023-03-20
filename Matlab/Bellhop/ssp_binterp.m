function [ C, gradc, Crr, Crz, Czz, Layer ]= ssp( x, zSSPV, cSSPV, czV, TopOpt ) 

% tabulates the sound speed profile and its derivatives
% also returns a vector Layer indicating the layer a depth point is in

global Layer 
persistent firstcall r z c cz cr crr crz czz


% Layer is a vector of indices indicating the layer that each ray is in
% [zSSPV, sSSPV] contains the depth/sound speed values
% cSSPV is a vector containing the 
% all vectors are set up as column vectors ...

'this version doesnt work; interp2 expects a matrix of points on which to interpolate'
if ( ~exist( 'firstcall' ) )
   load sspgrid r z c cz cr crr crz czz
   firstcall = 0;
   r = r';  % make it a row vector
end

Npts  = length( x( :, 1 ) );
Layer = ones( Npts, 1 );

C  = interp2( r, z, c,  x( :, 1 ), x( :, 2 ), 'cubic' );
Cz = interp2( r, z, cz, x( :, 1 ), x( :, 2 ), 'cubic' );

gradc = [ zeros(Npts, 1 ) cz ];
Crr   = zeros( Npts, 1 );
Crz   = zeros( Npts, 1 );
Czz   = zeros( Npts, 1 );
