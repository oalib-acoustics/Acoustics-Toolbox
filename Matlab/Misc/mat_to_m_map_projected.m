function [ Bathy ] = mat_to_m_map_projected( xyzfile )
% Loads the bathymetry file
% assumed to be a rectangular grid of data written in the standard GEODAS
% xyz form
% xyzfile should include the extension
%
% Returns a structure:
% Bathy.depth depths
% Bathy.Lat   latitudes
% Bathy.Lon   longitudes
% Bathy.Latkm mapping of lat to km assuming (x,y)=(0,0) corresponds to ( Lon(1), Lat(1) )
% Bathy.Lonkm
%
% mike porter 10/2011 based on earlier LoadBathymetry.m

load ../SRTM/SRTM_6min;

% plot the site

figure
% Indian Ocean
m_proj( 'Gnomonic','longitude', 60, ...
                   'latitude', -25, 'radius', 45 );
R_earth = 6370.; % approximate radius of the earth in km

%m_proj( 'UTM','longitudes', [  30 90 ], ...
%              'latitudes',  [ -70 20 ] );

                % Taiwan
% m_proj( 'UTM','longitudes', [ -71 -70 ], ...
%               'latitudes', [  42  43 ] );

% convert the source coordinate (Heard Island)
% [ src.x, src.y ] = m_ll2xy( 74.5, -53.4 );

% convert the source coordinate (HA04 Crochet Islands)
[ src.x, src.y ] = m_ll2xy( 51.8978, -46.8433 );

src.x = src.x * R_earth
src.y = src.y * R_earth

m_coast( 'linewidth', 4,'color','k');
m_coast( 'patch', [.7 .7 .7],'edgecolor','none');

m_grid % ( 'box', 'on' );

% convert lat/longs to x-y

[ Lon, Lat ] = meshgrid( Bathy.Lon, Bathy.Lat );
[ Bathy.X, Bathy.Y ] = m_ll2xy( Lon, Lat );
Bathy.X = R_earth * Bathy.X;
Bathy.Y = R_earth * Bathy.Y;

% convert grid coordinates to lat/longs, using the box limits suggested by
% the UTM mapping

grid.nx = 101;
grid.ny = 101;

grid.x.min = min( min( Bathy.X ) );
grid.x.max = max( max( Bathy.X ) );
grid.y.min = min( min( Bathy.Y ) );
grid.y.max = max( max( Bathy.Y ) );

grid.x = linspace( grid.x.min, grid.x.max, grid.nx );
grid.y = linspace( grid.y.min, grid.y.max, grid.ny );

[ grid.X, grid.Y ] = meshgrid( grid.x, grid.y );

% plot in km

% convert the x-y grid to lat/longs
[ grid.lon, grid.lat ] = m_xy2ll( grid.X / R_earth, grid.Y / R_earth );

% interpolate the bathymetry at those points
Bathy.zI = interp2( Lon, Lat, Bathy.depth, grid.lon, grid.lat );

tej = flipud( colormap( jet( 256 ) ) );

figure
surf( Bathy.Lon, Bathy.Lat, Bathy.depth )
shading flat
colormap( tej )
colorbar
xlabel( 'Longitude' )
ylabel( 'Latitude' )
title( 'Bathymetry' )

%surf( grid.x, grid.y, Bathy.zI )
%shading interp

%
figure
pcolor( grid.x, grid.y, Bathy.zI )
shading flat
colormap( tej )
colorbar
xlabel( 'Eastings (km)' )
ylabel( 'Northings (km)' )
title( 'Bathymetry' )
axis equal

surf( grid.x, grid.y, Bathy.zI )
shading interp
colormap( tej )
colorbar
hold on
xlabel( 'Eastings (km)' )
ylabel( 'Northings (km)' )
zlabel( 'Depth (m)' )
axis image

earthbrown = [ 0.5 0.3 0.1 ];
% h          = fill( r, z, earthbrown );

set( gca, 'ZDir', 'Reverse' )   % plot with depth-axis positive down

%if ( nargout == 1 )
%   varargout( 1 ) = { h };   % return a handle to the figure
%end

% write it out in the standard BELLHOP3D format

Bathy.X     = grid.X;
Bathy.Y     = grid.Y;
Bathy.depth = Bathy.zI;

writebty3d( 'foo.bty', Bathy )

