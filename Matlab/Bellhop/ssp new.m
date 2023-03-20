function [ C, gradc, Crr, Crz, Czz, Layer ]= ssp( x, zSSPV, cSSPV, czV, TopOpt ) 

global Layer 
persistent firstcall r z c cz cr crr crz czz

% ssp
% tabulates the sound speed profile and its derivatives
% also returns a vector Layer indicating the layer a depth point is in

% Layer is a vector of indices indicating the layer that each ray is in
% [zSSPV, sSSPV] contains the depth/sound speed values
% cSSPV is a vector containing the 
% all vectors are set up as column vectors ...

if ( ~exist( 'firstcall' ) )
   load sspgrid r z c cz cr crr crz czz
   firstcall = 0;
   r = r';  % make it a row vector
end

Npts  = length( x( :, 1 ) );
Layer = ones( Npts, 1 );

C  = interp1( z, c,  x( :, 2 ), 'cubic', 'extrap' );
Cz = interp1( z, cz, x( :, 2 ), 'cubic', 'extrap' );

gradc = [ zeros(Npts, 1 ) cz ];
Crr   = zeros( Npts, 1 );
Crz   = zeros( Npts, 1 );
Czz   = zeros( Npts, 1 );


