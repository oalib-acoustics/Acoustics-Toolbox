btyfil = 'Seamount.bty';

alpha0 = 15;

%%
rmax   = 15;

xmin = -rmax;
xmax =  rmax;

ymin = -rmax;
ymax =  rmax;

%%
dxy = 0.5;

x  = ( xmin : dxy : xmax );
nx = length( x );

y  = ( ymin : dxy : ymax );
ny = length( y );

%%
% different construction ...
nx = 500;
ny = 500;

%nx = 30;
%ny = 30;

x = linspace( xmin, xmax, nx );
y = linspace( ymin, ymax, ny );
%%
z = zeros( ny, nx );
[ X, Y ] = meshgrid( x, y );

R = sqrt( X.^2 + Y.^2 );
z = R * tand( alpha0 ) * 1000;

%%
Bdry.X = x;
Bdry.Y = y;
Bdry.depth = z;

interp_type = 'C';
writebdry3d( btyfil, interp_type, Bdry )

plotbdry3d( btyfil )