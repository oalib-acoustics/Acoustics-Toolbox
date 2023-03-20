
seamount_parameters

btyfil = 'DoubleSeamount3D.bty';
interp_type = 'R';
interp_type = 'C';

%%
rmax   = 8000;

dxy = 0.5;
xmin = -rmax;  xmax = rmax;
ymin = -rmax;  ymax = rmax;

x  = ( xmin : dxy : xmax );
nx = length( x );
y  = ( ymin : dxy : ymax );
ny = length( y );

%% different construction ...

nx = 1000;
ny = 1000;
x = linspace( xmin, xmax, nx );
y = linspace( ymin, ymax, ny );

%%
z = zeros( ny, nx );
[ X, Y ] = meshgrid( x, y );

z = sub_seamount( seamount, X, Y );

Bathy.X = x / 1000;
Bathy.Y = y / 1000;
Bathy.depth = z;

writebdry3d( btyfil, interp_type, Bathy )
