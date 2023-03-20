function varargout = plotbty( btyfil )

% plot the bathymetry file
% usage: plotbty( btyfil )
% where btyfil is the BaThYmetry file (the extension is optional)
% e.g. plotbty( 'foofoo' )
%
% MBP July 1999

global xBot NbtyPts
global units

BotBTY = '*';   % flag to read from bty file
depthB = inf;
rBox   = inf;
readbty( btyfil, BotBTY, depthB, rBox ) % read the bathymetry data

% copy, removing +/- inf values
NbtyPts = NbtyPts - 2;
r = xBot( 1, 2 : end - 1 );
z = xBot( 2, 2 : end - 1 );

% set labels in m or km
xlab     = 'Range (m)';
if ( strcmp( units, 'km' ) )
  r      = r / 1000.0;
  xlab   = 'Range (km)';
end
%%

hold on
xlabel( xlab )
ylabel( 'Depth (m)' )

% close the polygon
zmax      = max( z );
thickness = ( max( z ) - min( z ) );
thickness = max( thickness, 1 );   % make sure there's some thickness

r( NbtyPts + 1 ) = r( NbtyPts );
z( NbtyPts + 1 ) = zmax + thickness;
r( NbtyPts + 2 ) = r( 1 );
z( NbtyPts + 2 ) = zmax + thickness;

earthbrown = [ 0.5 0.3 0.1 ];
h          = fill( r, z, earthbrown );

set( gca, 'YDir', 'Reverse' )   % plot with depth-axis positive down

if ( nargout == 1 )
   varargout( 1 ) = { h };   % return a handle to the figure
end

