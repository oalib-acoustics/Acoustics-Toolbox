function plotshdpol( varargin )

% plot a TL surface in dB (polar coordinates)
% usage:
% plotshdpol( filename )
%    or
% plotshdpol( filename, xs, ys, rd )
% where
%   xs, ys is the source coordinate in km
%   rd is the receiver depth in m

global units

filename = varargin{1};

if ( nargin == 2  || nargin == 3 )
   % Generate a warning.
   Message = 'Call plotshdpol with 1 or 4 inputs.';
   warning( Message );
end

if nargin > 1
   xsvec = varargin{ 2 };
   ysvec = varargin{ 3 };
   rd    = varargin{ 4 };
else
   xsvec = NaN;
   ysvec = NaN;
   rd    = 0.0;
end

% open the file and read data

isz = 1;   % select a source depth

for xs =  xsvec
   for ys = ysvec
      [ PlotTitle, ~, freqVec, ~, ~, Pos, pressure ] = read_shd( filename, xs, ys );
      PlotTitle = replace( PlotTitle, '_', ' ' );  % Get rid of underlines (for Laurel who uses them in titles)

      % get nrz so that tlt doesn't loose the singleton dimension when nrz = 1
      nrz = length( Pos.r.z );
      clear tlt
      tlt( :, 1 : nrz, : ) = abs( pressure( :, isz, :, : ) );   % taken chosen source depth
      %tlt( :, 1 : nrd, : ) = abs( pressure( 1, 1, :, isz, :, : ) );   % take chosen source depth
      
      tlt = permute( tlt, [ 2 1 3 ] );   % order so that TL( rd, x, y )
      
      % interpolate the TL field at the receiver depth
      % note: interp1 won't interpolate a vector with only 1 element
      
      if ( length( Pos.r.z ) == 1 )
         tl = squeeze( tlt( 1, :, : ) );
      else
         tl = squeeze( interp1( Pos.r.z, tlt, rd ) );
      end
      
      tl( isnan( tl ) ) = 1e-6;   % remove NaNs
      tl( isinf( tl ) ) = 1e-6;   % remove infinities
      tl( tl < 1e-37  ) = 1e-37;  % remove zeros
      
      tl = -20.0 * log10( tl );
      
      % if full circle, duplicate the first bearing
      
      ntheta = length( Pos.theta );
      d_theta = ( Pos.theta( end ) - Pos.theta( 1 ) ) / ( ntheta - 1 );
      
      if ( mod( Pos.theta( end ) + d_theta - Pos.theta( 1 ) + .001, 360.0 ) < .002 )
         Pos.theta( end + 1 ) = Pos.theta( end ) + d_theta;
         tl( end + 1, : ) = tl( 1, : );
      end
      tl = tl';
      
      % make plot polar
      
      [ th, r ] = meshgrid( Pos.theta, Pos.r.r );
      
      th        = ( 2 * pi / 360. ) * th;   % convert to radians
      [ x, y ]  = pol2cart( th, r );
      
      % offset grid by the source coordinate
      if ( isnan( xs ) && isnan( ys ) )
         xst = Pos.s.x / 1000;
         yst = Pos.s.y / 1000;
      else
         xst = xs;
         yst = ys;
      end
      
      xt = x + 1000. * xst * ones( size( x ) );
      yt = y + 1000. * yst * ones( size( x ) );
         
      if ( strcmp( units, 'km' ) )
         xt = xt / 1000;   % convert to km
         yt = yt / 1000;
      end

      % *** plot ***
      
      tej = flipud( jet( 256 ) );  % 'jet' colormap reversed
      %tej = flipud( parula( 256 ) );  % 'parula' colormap reversed
      
      surf( xt, yt, tl ); shading interp
      % surfc( x, y, tl ); shading interp
      % pcolor( x, y, tl ); shading flat
      
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

      title( { deblank( PlotTitle ); [ ...
         'Freq = '     num2str( freqVec( 1 )   ) ' Hz   ' ...
         'x_{src} = '  num2str( xst            ) ' km   ' ...
         'y_{src} = '  num2str( yst            ) ' km   ' ...
         'z_{src} = '  num2str( Pos.s.z( isz ) ) ' m     ' ...
         'z_{rcvr} = ' num2str( rd             ) ' m' ] } )
       
% axis( 'image' )
      drawnow
      hold on
   end
end

%%

% fixed size for publications
% set( gca, 'ActivePositionProperty', 'Position', 'Units', 'centimeters' )
% set( gcf, 'Units', 'centimeters' )
% set( gcf, 'PaperPositionMode', 'auto');   % this is important; default is 6x8 inch page
% set( gca, 'Position', [ 1    0                       14.0       7.0 ] )
% set( gcf, 'Units', 'centimeters' )
% set( gcf, 'Position', [ 1 15 17.0 7.5 ] )

%     set( gcf, 'Units', 'centimeters' )
%     set( gcf, 'PaperPositionMode', 'manual' );
%     set( gcf, 'PaperPosition', [ 3 3 15.0 10.0 ] )


