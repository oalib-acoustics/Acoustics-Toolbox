clear
btyfil = 'Wedge3d.bty';
interp_type = 'R';

zmax = 380;
zmin = 20;
slope_range = 7200;
wedge_angle = atand((zmax-zmin)/slope_range);
fprintf('Wedge angle = %.4f degrees \n',wedge_angle)

% bottom displacement (see Buckingham 1987)
f = 25.0   % frequency
omega = 2 * pi * f;
c0 = 1500;
c1 = 1700;
k1 = omega / c0;
alpha_c = acos( c0 / c1 )

delta = pi / ( 2 * k1 * sin( alpha_c ) )
delta_z = delta * cosd( wedge_angle )

dx = 1.0;
dy = 0.5;


xctr = 0;
xmin = -29.7;
xmax = +29.7;
x = (xmin:dx:xmax)+xctr;
nx = length(x);

yctr = dy/2;
ymin = -29.7;
ymax = 29.7;
y = (ymin:dy:ymax)+yctr;
ny = length(y);

z = zeros( ny, nx );
z0 = 200;   % Depth at (0,0)
zy = z0-y*1000*tand(wedge_angle);
zy = max(zy,zmin);
zy = min(zy,zmax);
% zy = fliplr(zy);
z = repmat(zy(:),1,nx);
z = z + delta_z;

fid = fopen( btyfil, 'w' );
fprintf( fid, '''%c'' \n', interp_type );

fprintf( fid, '%i \r\n', nx );
fprintf( fid, '%f %f /', xmin, xmax );
fprintf( fid, '\r\n');

fprintf( fid, '%i \r\n', ny );
fprintf( fid, '%f %f /', ymin, ymax );
fprintf( fid, '\r\n');

for iy = 1 : ny
   fprintf( fid, '%f ', z( iy, : ) );
   fprintf( fid, '\r\n');

end

fclose( fid );

plotbdry3d(btyfil)

