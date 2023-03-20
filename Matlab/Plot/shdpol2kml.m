function plotshdpol( filename )

% convert a shd file to a single TL surface in dB (polar coordinates),
% written in kml format for Google Earth
%
% usage: shdpol2kml( filename )
% mike porter, October 2011

'this is just barely started ...'

global units jkpsflag

% open the file and read data

xs =  0;   % km !!!
ys =  0;

for xs = 950 : 100 : 1050 % 0 : 100 : 100
   for ys = 750 : 200 : 1050
      [ PlotTitle, ~, ~, ~, Pos, pressure ] = read_shd( filename, xs, ys );
      
      pressure = squeeze( pressure( :, 1, 9, : ) );   % take first source and receiver depth
      
      % make plot polar
      
      [ th, r ] = meshgrid( Pos.theta, Pos.r.r );
      th        = ( 2 * pi / 360. ) * th;   % convert to radians
      [ x, y ]  = pol2cart( th, r );
      
      x = x + 1000. * xs * ones( size( x ) );
      y = y + 1000. * ys * ones( size( x ) );
      
      if ( strcmp( units, 'km' ) )
         x = x / 1000;   % convert to km
         y = y / 1000;
      end
     
      tlt = abs( pressure );

      tlt( isnan( tlt ) ) = 1e-6;   % remove NaNs
      tlt( isinf( tlt ) ) = 1e-6;   % remove infinities
      tlt( tlt < 1e-37 ) = 1e-37;   % remove zeros

      tlt = -20.0 * log10( tlt );
      tlt = tlt';
      
      % *** plot ***
      
      tej = flipud( colormap( jet( 256 ) ) );
      surfc( x, y, tlt ); shading interp
      colormap( tej );
      colorbar
      % caxisrev( [ tlmin, tlmax ] )
      
      view( 2 )
      xlabel( 'Range, x (m)' )
      ylabel( 'Range, y (m)' )
      if ( strcmp( units, 'km' ) )
         xlabel( 'Range, x (km)' )
         ylabel( 'Range, y (km)' )
      end
      
      zlabel( 'Depth (m)' )
      
      title( deblank( PlotTitle ) )
      axis( 'equal' )
      drawnow
      pause( 5 )
      hold on
   end
end
