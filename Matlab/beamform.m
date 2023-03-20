function beamform( filename )
%BEAMFORM simple beamformer
% mpb 2 March 2001

isd = 1      % select the index of the source depth
ird = 1:101  % select the indices of the pressure field

SL  = 150    % source level (dB)
NL  = 0      % noise  level (dB)

% read in the field data
[ pltitl, freq, ~, ~, ~, ~, rd, rr, pressure ] = read_shd( filename );

% set up replica vectors
phone_coords = rd( ird );
angles       = -90:1:90;

e = planewave_rep( phone_coords, angles, freq );   % construct steering vectors

power = 20 * log10 ( abs( e * squeeze( pressure( ird, : ) ) ) ) + SL - NL;
peak  = max( max( power ) );

figure
pcolor( rr, angles, power ); shading interp; colorbar
%imagesc( time, flipud( theta_con ), flipud( power ) ); colorbar
xlabel( 'Range (m)' )
ylabel( 'conical angle' )
title( pltitl )
%axis( [ 0 5000 -90 90 ] )


