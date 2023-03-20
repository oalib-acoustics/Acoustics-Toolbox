btyfil = 'pwedge3d.bty';
interp_type = 'R';

th0 = 1.5;  % Wedge angle (deg)

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
ydepth = [ -y*1000*tand(th0)];

z = repmat(ydepth(:),1, length(x));

fid = fopen( btyfil, 'w' );
fprintf( fid, '''%c'' \n', interp_type );
fprintf( fid, '%i \r\n', ny );
fprintf( fid, '%f %f /', ymin, ymax );
fprintf( fid, '\r\n');

fprintf( fid, '%i \r\n', nx );
fprintf( fid, '%f %f /', xmin, xmax );
fprintf( fid, '\r\n');

for ix = 1 : nx
   fprintf( fid, '%f ', z( :, ix ) );
   fprintf( fid, '\r\n');
end

fclose( fid );

figure
plotbdry3d( btyfil )
