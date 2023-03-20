function plotmode( filename, freq, modes )

% plot the modes produced by KRAKEN
%
% usage: plotmode( filename, freq, modes )
% where
%    filename is the mode file (extension is optional to indicate binary or ascii mode file)
%    freq is the frequency (needed if there are multiple modesets
%       corresponding to different frequencies)
%    modes is a vector of indices
%      can be specified in any order and with duplicates
%
% Examples:
%    plotmode( 'duct' )      % plots all the modes from duct.mod
%    plotmode( 'duct.mod' )  %   "
%    plotmode( 'duct', 1 )   % plots mode 1
%    plotmode( 'duct', [ 1 5 9 11 ] ); % plots the selected modes
%    plotmode duct           % a sort of sloppy syntax that Matlab accepts
%
% mbp

clear read_modes_bin % to force rewind to beginning of mode file

% add a default extension '.mod' if it's missing
[ ~, ~, ext ] = fileparts( filename );

if ( isempty( ext ) )
   filename = [ filename '.mod' ];
end

%for ifreq = 1 : 1   % if broadband run, do a frequency loop
if nargin == 2
   [ Modes ] = read_modes( filename, freq );
   modes = 1 : length( Modes.k );
else
   [ Modes ] = read_modes( filename, freq, modes );
end
%end

if ( Modes.M == 0 )
   error( 'No modes in mode file' )
end
   
% identify the index of the frequency closest to the user-specified value
freqdiff = abs( Modes.freqVec - freq );
[ ~, freq_index ] = min( freqdiff );

% extract the specfied component from the stress-displacement vector
phi = get_component( Modes, 'N' );
% Modes.z = Modes.z( 1 : 101 );

% Extract stress-displacement components

% phiV = get_component( Modes, 'V' );
% phiH = get_component( Modes, 'H' );
% phiN = get_component( Modes, 'N' );
% phiT = get_component( Modes, 'T' );
%
% phiV = phiV( 602 : end, : );
% phiH = phiH( 602 : end, : );
% phiN = phiN( 602 : end, : );
% phiT = phiT( 602 : end, : );
% Modes.z = Modes.z( 602 : end );

%figure; plot( z, psi );
%xlabel( 'Depth (m)' )
%view( 90, 90 )	% rotate

% plots of the wavenumbers in the complex plane

figure
plot( real( Modes.k ), imag( Modes.k ), 'o' );
title( { deblank( Modes.title ); ...
   [ 'Freq = ' num2str( Modes.freqVec( freq_index ) ) ' Hz' ] } )
xlabel( 'real( k )' );
ylabel( 'imag( k )' );

%%

% plots of the modes as a single color image
% skip if there's only one mode selected

nx = size( Modes.phi, 2 );

if nx > 1
   x = 1 : nx;
   
   figure
   doo = double( real( phi ) );
   pcolor( x, Modes.z, doo ); shading flat
   set( gca, 'YDir', 'Reverse' )
   
   colormap( jet( 256 ) )
   colorbar
   caxis_lim = max( abs( caxis ) );      % get extrema of colorbar limits
   caxis( [ -caxis_lim, caxis_lim ] );   % make a colorbar symmetric about 0
   
   xlabel( 'Mode index' )
   ylabel( 'Depth (m)' )
   title( { deblank( Modes.title ); ...
      [ 'Freq = ' num2str( Modes.freqVec( freq_index ) ) ' Hz' ] } )
   
end

%%
% line plots of the modes

Nplots = min( length( modes ), 10 );
iskip  = floor( length( modes ) / Nplots );

figure

for iplot = 1 : Nplots
   subplot( 1, Nplots, iplot );
   imode = 1 + ( iplot - 1 ) * iskip;
   plot( real( phi( :, imode ) ), Modes.z )
   hold on
   plot( imag( phi( :, imode ) ), Modes.z, 'b--' )
   
   xlabel( [ 'Mode ' num2str( modes( imode ) ) ] )
   set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
   if ( iplot == 1 )
      ylabel( 'Depth (m)')
      title( { deblank( Modes.title ); ...
         [ 'Freq = ' num2str( Modes.freqVec( freq_index ) ) ' Hz' ] } )
   else
      set( gca, 'YTickLabel', [' ';' '] ) % no tick lables
      %set( gca, 'YTickMode', 'manual')
   end
   
   %    temp = axis;
   %    temp( 3 : 4 ) = [ 0 300 ];
   %    axis( temp )
end

%%
%
% % vertical component
% subplot( 1, 4, 1 );
% imode = 1;
% plot( real( phiV( :, imode ) ), Modes.z )
% hold on
% plot( imag( phiV( :, imode ) ), Modes.z, 'b--' )
%
% xlabel( [ 'u' ] )
% set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
% ylabel( 'Depth (m)')
%
% temp = axis;
% temp( 3 : 4 ) = [ 150 225 ];
% axis( temp )
%
% % horizontal component
% subplot( 1, 4, 2 );
%
% plot( real( phiH( :, imode ) ), Modes.z )
% hold on
% plot( imag( phiH( :, imode ) ), Modes.z, 'b--' )
%
% xlabel( [ 'w' ] )
% set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
%
% set( gca, 'YTickLabel', [' ';' '] ) % no tick lables
% %set( gca, 'YTickMode', 'manual')
% temp = axis;
% temp( 3 : 4 ) = [ 150 225 ];
% axis( temp )
%
%
% % normal stress
% subplot( 1, 4, 3 )
% plot( real( phiN( :, imode ) ), Modes.z )
% hold on
% plot( imag( phiN( :, imode ) ), Modes.z, 'b--' )
%
% xlabel( [ '\tau_{zz}' ] )
% set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
%
% set( gca, 'YTickLabel', [' ';' '] ) % no tick lables
% %set( gca, 'YTickMode', 'manual')
% temp = axis;
% temp( 3 : 4 ) = [ 150 225 ];
% axis( temp )
%
% % tangential stress
% subplot( 1, 4, 4 );
%
% plot( real( phiT( :, imode ) ), Modes.z )
% hold on
% plot( imag( phiT( :, imode ) ), Modes.z, 'b--' )
%
% xlabel( [ '\tau_{zx}' ] )
% set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
%
% set( gca, 'YTickLabel', [' ';' '] ) % no tick lables
% %set( gca, 'YTickMode', 'manual')
% temp = axis;
% temp( 3 : 4 ) = [ 150 225 ];
% axis( temp )
%
% %    temp = axis;
% %    temp( 3 : 4 ) = [ 0 300 ];
% %    axis( temp )
%
% %title( { deblank( Modes.title ); ...
% %       [ 'Freq = ' num2str( Modes.freq ) ' Hz' ] } )
% %figure; surf( modes, z, real( psi ) ); colorbar
