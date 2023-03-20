function plottl_zf( filename, rrt )

% plot TL from the shade file as a function of depth, z, and freq
%
% usage:
% plottl_zf( filename, rrt )
% where
%   filename is the shadefile (with extension)
%   rrt is the receiver range in km
%   if rrt is a vector then one plot is generated for each element
% mbp 1/2018

disp( 'PlotTL_zf uses the first bearing and source depth in the shade file; check OK' )
itheta = 1;
isz    = 1;

% read the frequency vector

[ PlotTitle, ~, freqVec, ~, Pos, ~ ] = read_shd( filename );

Nfreq = length( freqVec );
Nz    = length( Pos.r.z );

% loop to read the pressure for each frequency

TLmat = zeros( Nfreq, Nz );

for ifreq = 1 : length( freqVec )
   [ ~, ~, ~, ~, ~, pressure ] = read_shd( filename, freqVec( ifreq ) );
   
   pressure = pressure( itheta, isz, :, : );

   tlt = abs( pressure );	            % this is really the negative of TL
   tlt( tlt == 0 ) = max( max( tlt ) ) / 1e10;      % replaces zero by a small number
   tlt = -20.0 * log10( tlt );         % so there's no error when we take the log
   
   %figure; imagesc( squeeze( tlt ) ); colormap( flipud( jet ) ); caxisrev( [ 60 100 ] )

   % interpolate the TL field at the receiver range
   % note: interp1 won't interpolate a vector with only 1 element
   
   if ( length( Pos.r.r ) == 1 )
      tlslice = tlt;
   else
      TLtemp  = squeeze( tlt )';   % need to avoid dimensional problems when TLT has 1, 2, or 3 dimensions
      tlslice = interp1( Pos.r.r, TLtemp, 1000.0 * rrt );
   end
   TLmat( ifreq, : ) = tlslice;
   %figure; plot( tlslice ); axis( [ 0 65 60 100 ] )
end

%%
% plot TL vs frequency and depth

pcolor( freqVec, Pos.r.z, TLmat' );  ...
shading flat

tej = flipud( jet( 256 ) );  % 'jet' colormap reversed
colormap( tej )
colorbar
caxis( [ 60 100 ] )

%caxisrev( [ tlmin, tlmax ] )
set( gca, 'YDir', 'Reverse' )
xlabel( 'Frequency (Hz)' )
ylabel( 'Depth (m)' )

set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
title( { deblank( PlotTitle ); [ 'z_{src} = ' num2str( Pos.s.z( isz ) ) ' m    r_{rcvr} = ' num2str( rrt ) ' km' ] } )

