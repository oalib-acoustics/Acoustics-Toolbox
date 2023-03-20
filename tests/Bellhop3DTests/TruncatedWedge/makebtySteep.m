btyfil = 'steep.bty';
interp_type = 'R';

zmax = 4000;
zmin = 20;

slope_range = 7200;
wedge_angle = atand((zmax-zmin)/slope_range);
fprintf('Wedge angle = %.4f degrees \n',wedge_angle)

dx = 10.;
dy = 10.;

xctr = 0;
xmin = -30;
xmax = +30;
x = (xmin:dx:xmax)+xctr;
nx = length(x);

yctr = 0; % dy/2;
ymin = -30;
ymax = +30;
y = (ymin:dy:ymax)+yctr;
ny = length(y);

z = zeros( ny, nx );
z0 = 200;   % Depth at (0,0)
zy = z0 - y * 1000 * tand( wedge_angle );

% truncate the wedge
%zy = max(zy,zmin);
%zy = min(zy,zmax);

% zy = fliplr(zy);
z = repmat(zy(:),1,nx);

fid = fopen( btyfil, 'w' );
fprintf( fid, '''%c'' \n', interp_type );

fprintf( fid, '%i \r\n', nx );
fprintf( fid, '%f %f /', min( x ), max( x ) );
fprintf( fid, '\r\n');

fprintf( fid, '%i \r\n', ny );
fprintf( fid, '%f %f /', min( y ), max( y ) );
fprintf( fid, '\r\n');

for iy = 1 : ny
   fprintf( fid, '%f ', z( iy, : ) );
   fprintf( fid, '\r\n');
end

fclose( fid );
