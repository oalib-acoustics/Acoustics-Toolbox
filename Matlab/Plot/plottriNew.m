function plottri( filename )

% plot the triangulation for a 3D KRAKEN run
%
% There are certain restrictions relative to the documentation for the flp file
% Commas used to separate fields should be removed.

loadssp	% this is used to read in surface SSP (optional)

fid = fopen( [ filename '.flp' ] );

title = fgetl( fid );
opt   = fgetl( fid );
m     = fscanf( fid, '%i', 1 );
temp  = fgetl( fid );

xs  = fscanf( fid, '%f', 1 );
ys  = fscanf( fid, '%f', 1 );
temp  = fgetl( fid );

% receiver ranges

nsd	= fscanf( fid, '%i', 1 );
sd		= fscanf( fid, '%f', 1 );
temp = fgetl( fid );

nrd	= fscanf( fid, '%i', 1 );
rd		= fscanf( fid, '%f', 1 );
temp = fgetl( fid );

nr		= fscanf( fid, '%f', 1 );
temp = fgetl( fid );
temp	= fscanf( fid, '%f', 2 );
rmin	= temp( 1 );
rmax	= temp( 2 );

temp = fgetl( fid );

% angles

ntheta	= fscanf( fid, '%f', 1 );
temp		= fscanf( fid, '%f', 2 );
theta = linspace( temp( 1 ), temp( 2 ), ntheta );
temp = fgetl( fid );

% nodes
nnodes = fscanf( fid, '%i', 1 )
temp = fgetl( fid );

x = zeros( nnodes, 1 );
y = zeros( nnodes, 1 );

for inode = 1 : nnodes
   temp = fscanf( fid, '%f', 2 );
   x( inode ) = temp( 1 );
   y( inode ) = temp( 2 );
   nodename( inode, : ) = fgetl( fid );
   
   % following can be commented out if SSP plot is omitted
   ii = findstr( nodename( inode, 8 ), 'abcde' );
   jj = str2num( nodename( inode, 9:10 ) );
   ctop( inode ) = c( ii, jj );

end

%figure
hold on

% elements

nelts = fscanf( fid, '%i', 1 );
temp  = fgetl( fid );

for ielt = 1:nelts
   temp = fscanf( fid, '%i', 3 );
   xv = [ x( temp( 1 ) ) x( temp( 2 ) ) x( temp( 3 ) ) ];
   yv = [ y( temp( 1 ) ) y( temp( 2 ) ) y( temp( 3 ) ) ];
   plot( xv, yv )
   
   % following can be commented out if SSP plot is omitted
   cv = [ ctop( temp( 1 ) ) ctop( temp( 2 ) ) ctop( temp( 3 ) ) ];
   fill( xv', yv', cv' )
   temp = fgetl( fid );
end

% Plot the bearing lines

clear xv yv

deg2rad = pi / 180;

for ith = 1:ntheta
   xv( 1 ) = xs + rmin * cos( deg2rad * theta( ith ) );
   yv( 1 ) = ys + rmin * sin( deg2rad * theta( ith ) );
   xv( 2 ) = xs + rmax * cos( deg2rad * theta( ith ) );
   yv( 2 ) = ys + rmax * sin( deg2rad * theta( ith ) );
   bearing = plot( xv, yv, 'w' );
   set( bearing, 'LineWidth', 2 )

end

axis equal
xlabel( 'x (km)' );
ylabel( 'y (km)' );
colorbar
axis( [ 100 600 200 600 ] )
view( 0, -90 )

