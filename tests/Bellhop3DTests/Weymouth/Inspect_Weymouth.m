load WeymouthForeRiver

% from
%    http://www.nae.usace.army.mil/Missions/Navigation/Massachusetts-Projects/

gr = 10; % make 10-meter grid
x1 = 743825;
y1 = 452474;

y1m = -800;
y2m =  500;

x = ( x - x1 ) * 12 / 39.37;
y = ( y - y1 ) * 12 / 39.37; 
% state plane 1983 coordinate system, feet, add offset, convert to meters

z = z * 12 / 39.37; % convert feet to meters

ii = find (y > y1m & y < y2m);

x = x( ii );
y = y( ii );
z = z( ii );

% change small heights to a small number to make BELLHOP3D happy
ii = find( z > -.1 );
z( ii ) = -.1;

F = TriScatteredInterp( x, y, z ); % or scatteredinterpolator

xg = round( min( x ) ) : gr : round( max( x ) ); 
yg = round( min( y ) ) : gr : round( max( y ) );
[ X, Y ] = meshgrid( xg, yg );

bathy = F( X, Y );
bathy( isnan( bathy ) ) = -.1;   % remove NaNs

figure( 1 ); clf
surf( X, Y, bathy ); shading flat
set( gca, 'dataaspectratio', [ 1 1 .03 ])
grid on

xlabel( 'x (m)' )
ylabel( 'y (m)' )
zlabel( 'z (m)' )

figure( 2 ); clf
contourf( X, Y, bathy, -[0:2:16] ); 
set(gca,'dataaspectratio', [ 1 1 .03 ] )
grid on

xlabel( 'x (m)' )
ylabel( 'y (m)' )

%%
% mbp: added this section to write a bathymetry file in the format
% BELLHOP3D uses
btyfil = 'Weymouth.bty';
Bathy.X = xg / 1000.;
Bathy.Y = yg / 1000.;
Bathy.depth = -bathy;

writebdry3d( btyfil, Bathy )