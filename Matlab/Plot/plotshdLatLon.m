function plotshdLatLon( varargin )

% plot a TL surface in dB (lat/long coordinates)
% usage:
% plotshdpol( filename )
%
% plotshdpol( filename, xs, ys, rd )

global units jkpsflag

filename = varargin{1};

if nargin == 2
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
end

% initialize the projection
m_proj( 'Gnomonic','longitude', 60, ...
   'latitude', -25, 'radius', 30 );
R_earth = 6370.; % approximate radius of the earth in km

%m_coast( 'linewidth', 4,'color','k');
%m_coast( 'patch', [.7 .7 .7],'edgecolor','none');
%m_grid % ( 'box', 'on' );

hold on

% open the file and read data

for xs =  xsvec
   for ys = ysvec
      [ PlotTitle, ~, ~, ~, Pos, pressure ] = read_shd( filename, xs, ys );
      
      % get nrd so that tlt doesn't loose the singleton dimension when nrd = 1
      nrd = length( Pos.r.z );
      clear tlt
      tlt( :, 1 : nrd, : ) = abs( pressure( :, 1, :, : ) );   % take first source depth
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
      tl( tl < 1e-37  ) = 1e-37;   % remove zeros
      
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
      
      x = x + 1000. * xs * ones( size( x ) );
      y = y + 1000. * ys * ones( size( x ) );
      
      % map to lat-lon
      [ lon, lat ] = m_xy2ll( x / 1000 / R_earth, y / 1000 / R_earth );
          
      % *** plot ***
      
      tej = flipud( jet( 256 ) );
      surfc( lon, lat, tl ); shading interp
      colormap( tej );
      colorbar
      % caxisrev( [ tlmin, tlmax ] )
      
      view( 2 )
      xlabel( 'Longitude' )
      ylabel( 'Latitude' )      
      zlabel( 'Depth (m)' )
      
      title( deblank( PlotTitle ) )
      % axis( 'image' )
      drawnow
      hold on
   end
end


