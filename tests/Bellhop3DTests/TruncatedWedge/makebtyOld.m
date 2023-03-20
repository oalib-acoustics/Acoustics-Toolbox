clear
btyfil = 'Wedge3d.bty';
interp_type = 'R';

zmax = 380;
zmin = 20;
slope_range = 7200;
wedge_angle = atand((zmax-zmin)/slope_range);
fprintf('Wedge angle = %.4f degrees \n',wedge_angle)

xctr = 0;
yctr = 0;
z0 = 200;   % Depth at (0,0)

xmin = -29.7;
xmax = +29.7;
dx = 0.1;
x = (xmin:dx:xmax);
nx = length(x);

ymin = -29.7;
ymax = 29.7;
dy = 0.1;
y = (ymin:dy:ymax);
ny = length(y);

z = zeros( ny, nx );

zy = z0-y*1000*tand(wedge_angle);
zy = max(zy,zmin);
zy = min(zy,zmax);
% zy = fliplr(zy);
z = repmat(zy(:),1,nx);

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

