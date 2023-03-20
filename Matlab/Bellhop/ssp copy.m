function [ c, gradc, crr, crz, czz, Layer ]= ssp( x, SSP, Layer )

% ssp
% tabulates the sound speed profile and its derivatives
% also returns a vector Layer indicating the layer a depth point is in

% Layer is the index of the layer that each ray is in
% [SSP.z, SSP.c] contains the depth/sound speed values
% all vectors are set up as column vectors ...
% x can be a vector or a scalar

%Layer = find( x( :, 2 ) > zSSPV( : ) );

NSSPLayers = length( SSP.z ) - 1;

% search through deeper layers
while ( x( 2 ) >= SSP.z( Layer + 1 ) && Layer < NSSPLayers )
   Layer = Layer + 1;
end

% search through shallower layers
while ( x( 2 ) < SSP.z( Layer ) && Layer > 1 )
   Layer = Layer - 1;
end

w  = x( 2 ) - SSP.z( Layer ); % distance into layer
cz = SSP.cz( Layer );
c  = SSP.c(  Layer ) + w * cz;

%c  = interp1( zSSPV, cSSPV, x( :, 2 ), 'spline', 'extrap' );
%cz = interp1( zSSPV( 1:end-1) + HV, czV, x( :, 2 ), 'spline', 'extrap' );

gradc = [ 0 cz ];
crr   = 0;
crz   = 0;
czz   = 0;
