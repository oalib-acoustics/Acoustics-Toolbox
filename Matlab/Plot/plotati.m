function varargout = plotati( atifil )

% plot the altimetry file
% usage: plotati( atifil )
% where atifil is the Altimetry file (the extension is optional)
% e.g. plotati( 'foofoo' )
%
% plots the AlTImetry file used by Bellhop
% MBP July 1999

global xTop NatiPts
global units

TopATI = '*';   % flag to read from bty file
depthT = -inf;
rBox   = inf;
readati( atifil, TopATI, depthT, rBox ) % read the bathymetry data

% copy, removing +/- inf values
NatiPts = NatiPts - 2;
r = xTop( 1, 2 : end - 1 );
z = xTop( 2, 2 : end - 1 );

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
zmin = min( z );
thickness = ( max( z ) - min( z ) );
thickness = max( thickness, 1 );   % make sure there's some thickness

r( NatiPts + 1 ) = r( NatiPts );
z( NatiPts + 1 ) = zmin - thickness;
r( NatiPts + 2 ) = r( 1 );
z( NatiPts + 2 ) = zmin - thickness;

skyblue = [ 0.7 0.7 1.0 ];
h       = fill( r, z, skyblue );

set( gca, 'YDir', 'Reverse' )   % plot with depth-axis positive down

if ( nargout == 1 )
   varargout( 1 ) = { h };   % return a handle to the figure
end

