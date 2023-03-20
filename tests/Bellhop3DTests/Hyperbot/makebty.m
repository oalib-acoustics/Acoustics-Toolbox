% generate hyperbolic bottom bathymetry
% mbp 30 July 2019

btyfil = 'HyperBot.bty';
interp_type = 'C';

xmin = -20;
xmax = 20;
ymin = -20;
ymax = 20;

nx = 510;
ny = 510;
x = linspace( xmin, xmax, nx );
y = linspace( ymin, ymax, ny );

[ X, Y ] = meshgrid( x, y );
R = 1000 * sqrt( X.^2 + Y.^2 );   % range in meters

c = 500;
z = sqrt( c^2 + ( R / 10 ).^2 );

% bathymetry

fid = fopen( btyfil, 'w' );
fprintf( fid, '''%c'' \n', interp_type );

fprintf( fid, '%i \r\n', nx );
fprintf( fid, '%f ', x );
%fprintf( fid, '%f %f /', xmin, xmax );
fprintf( fid, '\r\n');

fprintf( fid, '%i \r\n', ny );
fprintf( fid, '%f ', y );
%fprintf( fid, '%f %f /', ymin, ymax );
fprintf( fid, '\r\n');

for iy = 1 : ny
   fprintf( fid, '%f ', z( iy, : ) );
   fprintf( fid, '\r\n');
end

fclose( fid );
