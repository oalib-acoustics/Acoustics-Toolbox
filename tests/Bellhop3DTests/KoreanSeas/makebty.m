btyfil = 'KoreanSea.bty';
interp_type = 'R';
interp_type = 'C';

% Write a bathymetry from the workspace variables

% write bty file for BELLHOP3D
% mbp 3/2011

xctr = 30;
yctr = 30;

xmin = -100;
xmax = +100;
nx   = 41;

ymin = -100;
ymax = +100;
ny   = 41;

x = linspace( xmin, xmax, nx );
y = linspace( ymin, ymax, ny );

z = zeros( nx, ny );

for ix = 1 : nx
   for iy = 1 : ny
      r2 = ( x( ix ) - xctr )^ 2 + ( y( iy ) - yctr )^ 2;
      z( ix, iy ) = 5000.0 - 2000 / ( 1 + r2 / 20^2 );
      % z( ix, iy ) = 1000.0;
   end
end

Bathy.X = x;
Bathy.Y = y;
Bathy.depth = z;

writebdry3d( btyfil, Bathy )

