btyfil = 'wedge2d.bty';
interp_type = 'R';

th0 = 1.2;  % Wedge angle (deg)

% Write a bathymetry from the workspace variables

% write bty file for BELLHOP3D
% xaz  1/2012

xctr = 0;
yctr = 0;

xmin = -120;
xmax = +120;
nx   = 21;

ymin = -120;
ymax = 120;
ny   = 21;

x = linspace( xmin, xmax, nx );
y = linspace( ymin, ymax, ny );

z = zeros( nx, ny );

% ydepth = [0 range(y)*1000*tand(th0)];
ydepth = -y * 1000 * tand( th0 );

z = repmat( ydepth(:), 1, length(x) );

Bathy.X = x;
Bathy.Y = y;
Bathy.depth = z;

writebdry3d( btyfil, interp_type, Bathy )

figure
plotbdry3d( btyfil )
axis( [ -100 100 -40 0 0 900 ] )
caxis( [ 0 900 ] )
print -dpng Bathy

copyfile( 'wedge2d.bty', 'wedge3dHat.bty' )
copyfile( 'wedge2d.bty', 'wedge3dGaussian.bty' )
