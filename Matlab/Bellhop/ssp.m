function [ c, c_imag, gradc, crr, crz, czz, Layer ]= ssp( x, SSP, Layer )

% tabulates the sound speed profile and its derivatives
% also returns a vector Layer indicating the layer a depth point is in

% Layer is the index of the layer that each ray is in
% [SSP.z, SSP.c] contains the depth/sound speed values
% all vectors are set up as column vectors ...
% x can be a vector or a scalar

NSSPLayers = length( SSP.z ) - 1;

% following lines are equivalent to
% Layer = find( x( :, 2 ) > zSSPV( : ) );

% search through deeper layers

while ( ( x( 2 ) >= SSP.z( Layer + 1 ) ) && ( Layer < NSSPLayers ) )
   Layer = Layer + 1;
end

% search through shallower layers
while ( ( x( 2 )  < SSP.z( Layer     ) ) && ( Layer > 1          ) )
   Layer = Layer - 1;
end

w = x( 2 ) - SSP.z( Layer ); % distance into layer

switch ( SSP.Type )
   case 'N'   % n^2 Linear interpolation
      n2z    = SSP.n2z( Layer );
      c_eval = 1.0 / sqrt( SSP.n2( Layer ) + w * n2z );
      c      = real( c_eval );
      c_imag = imag( c_eval );

      gradc  = [ 0.0 real( -0.5 * c_eval * c_eval * c_eval * n2z ) ];
      czz    = 3 * gradc( 2 ) * gradc( 2 ) / c;

   case { 'P', 'S' }   % monotone PCHIP or natural Cublic Spline interpolation
      c_eval = ppval( SSP.ppc, x( 2 ) );
      c      = real( c_eval );
      c_imag = imag( c_eval );

      gradc  = [ 0.0 real( ppval( SSP.ppcz,  x( 2 ) ) ) ];
      czz    =       real( ppval( SSP.ppczz, x( 2 ) ) );
   otherwise   % c Linear interpolation
      cz     = SSP.cz( Layer );
      c_eval = SSP.c( Layer ) + w * cz;

      c      = real( c_eval );
      c_imag = imag( c_eval );

      gradc  = [ 0.0 real( cz ) ];
      czz    = 0.0;
end

crr = 0.0;
crz = 0.0;

return


