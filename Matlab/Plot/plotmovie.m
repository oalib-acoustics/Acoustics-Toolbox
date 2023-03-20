function plotmovie( filename )

% plot a single TL surface in dB
% usage: plotmovie filename
% where filename is a .mat file without the extension

dbflag = 0;

[ PlotTitle, ~, freqVec, ~, ~, Pos, p ] = read_shd( filename );

p = squeeze( p );

Nfreq = length( freqVec );
nsd = length( Pos.s.z );

% figure( 'units', 'normalized', 'outerposition', [ 0 0 0.5 0.5 ] );  % create full screen figure
% newmap = colormap( jet( 256 ) );
% newshortmap = newmap( 1 : end, : );
tej = flipud( jet( 256 ) );  % 'jet' colormap reversed

if ( ~dbflag )
   tej( 110 : 146, : ) = ones( 37, 3 );   % no color for points near zero
end

figure

% set( gca, 'ActivePositionProperty', 'Position', 'Units', 'centimeters' )
% set( gcf, 'Units', 'centimeters' )
% set( gcf, 'PaperPositionMode', 'auto');   % this is important; default is 6x8 inch page
% 
% 
% set( gca, 'Position', [ 2 0.5   14.0 7.0 ] )
% set( gcf, 'Units', 'centimeters' )
% set( gcf, 'Position', [ 3 15 17.0 8.0 ] )

%%
for ifreq = 1 : Nfreq
   tlt = squeeze( real( p( ifreq, :, : ) ) );
   
   if dbflag
      tlt = -20 * log10( abs( tlt ) );
   else
      tlt = 1e6 * tlt;   % pcolor routine has problems when the values are too low
      
      %tlt( :, 1 ) = zeros( nrd, 1 );   % zero out first column for SPARC run
      tlmax = max( max( abs( tlt ) ) );
      tlmax = 0.2 * max( tlmax, 0.000001 );
      % tlmax = maxtlt/10;
      % tlmax = 0.02/i;
   end
   
   pcolor( Pos.r.r, Pos.r.z, tlt )
   % imagesc( Pos.r.r, Pos.r.z, tlt )
   shading interp
   colormap( tej )
   axis image
   xlabel( 'Range (m)' )
   ylabel( 'Depth (m)' )
   title( { deblank( PlotTitle ); [ 'Time = ' num2str( freqVec( ifreq ) ) ' s' ] } )
   set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
   
    if dbflag
       caxisrev( [ 50 110 ] )
    else
       caxis( [ -tlmax, tlmax ] )
       % rectangles:
%        hold on
%        
%        % air
%        skyblue = [ 0.7 0.9 1.0 ];
%        
%        r = [    0 1000 1000 0    0 ]';
%        z = [ -100 -100    0 0 -100 ]';
%        
%        h = fill( r, z, skyblue' );
%        alpha( h, .2 )
%        
%        % water
%        seablue = [ 0.5 0.5 1.0 ];
%       r = [    0 1000 1000 0    0 ]';
%       z = [  200  200    0 0  200 ]';
%       
%       h = fill( r, z, seablue );
%       alpha( h, .2 )
    end
   colorbar
   
   drawnow
   
%    if ( mod( isd, 10 ) == 0 )
%       eval( [ 'print -dpng ' num2str( isd ) ] )
%    end
   
   %    if isd == 1
   %       M = moviein( nsd ); % initialize movie storage
   %    end
   %
   %    A = getframe( gcf );
   %    M( :, isd ) = A;
   
end

% save movie M
% movie2avi( M, 'movie', 'Quality', 100 );
