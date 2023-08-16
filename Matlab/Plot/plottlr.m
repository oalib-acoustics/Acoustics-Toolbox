function plottlr( filename, rdt )

% plot a single TL slice from the shade file
%
% usage:
% plottlr( filename, rdt )
% where
%   filename is the shadefile (with extension)
%   rdt is the receiver depth in m
%   if rdt is a vector then one plot is generated for each element
% mbp

global units

disp( 'PlotTLr uses the first bearing and source depth in the shade file; check OK' )
itheta = 1;
isz    = 1;
ifreq  = 2
freq   = 100   % 

% read

[ PlotTitle, ~, freqvec, ~, ~, Pos, pressure ] = read_shd( filename );
[ PlotTitle, ~, freqvec, ~, ~, Pos, pressure ] = read_shd( filename, freq );
%freq = freqvec( ifreq );
rkm = Pos.r.r / 1000.0;         % convert to km

pressure = pressure( itheta, isz, :, : );

tlt = abs( pressure );	            % this is really the negative of TL
tlt( tlt == 0 ) = max( max( tlt ) ) / 1e10;      % replaces zero by a small number
tlt = -20.0 * log10( tlt );          % so there's no error when we take the log

% interpolate the TL field at the receiver depth
% note: interp1 won't interpolate a vector with only 1 element

if ( length( Pos.r.z ) == 1 )
  tlslice = squeeze( tlt( 1, : ) );
  rdt     = Pos.r.z( 1 );
else
  TLtemp  = squeeze( tlt );   % need to avoid dimensional problems when TLT has 1, 2, or 3 dimensions
  tlslice = interp1( Pos.r.z, TLtemp, rdt );
end

hh = plot( rkm, tlslice', 'b' );

set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
xlabel( 'Range (km)' )
ylabel( 'TL (dB)' )
title( { deblank( PlotTitle ); [ 'Freq = ' num2str( freq ) ' Hz    Sz = ' num2str( Pos.s.z( isz ) ) ' m' ] } )

set( hh, 'LineWidth', 2 )

% generate legend
for irz = 1: length( rdt )
    legendstr( irz, : ) = [ 'Rz = ', num2str( rdt( irz ) ), ' m' ];
end

legend( legendstr, 'Location', 'Best' )
legend( 'boxoff' )
drawnow

% %figure; plot( rkm, abs( interp1( Pos.r.z, squeeze( pressure ), rdt ) ) );
