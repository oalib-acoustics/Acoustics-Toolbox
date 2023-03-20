function varargout = plotbdry3d( bdryfil )

% plot the boundary file
% usage: plotbdry3D( bdryfil )
% where bdryfil is the boundary file
% e.g. plotbdry( 'foofoo.bty' )
%
% plots the boundary file used by Bellhop3D
% MBP April 2011

global units

[ x, y, z, ~, ~ ] = readbdry3d( bdryfil );

% set labels in m or km
xlab = 'Range-x (m)';
ylab = 'Range-y (m)';
x    = x * 1000.0;
y    = y * 1000.0;

if ( strcmp( units, 'km' ) )
   x    = x / 1000.0;
   y    = y / 1000.0;
   xlab = 'Range-x (km)';
   ylab = 'Range-y (km)';
end

[ X, Y ] = meshgrid( x, y );
surf( X, Y, z )
shading faceted
shading interp
colormap( flipud( jet ) )
colorbar

%%

hold on
xlabel( xlab )
ylabel( ylab )
zlabel( 'Depth (m)' )

set( gca, 'ZDir', 'Reverse' )   % plot with depth-axis positive down

if ( nargout == 1 )
   varargout( 1 ) = { h };   % return a handle to the figure
end

