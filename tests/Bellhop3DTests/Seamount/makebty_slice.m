% creates a BELLHOP bathymetry from a single slice

btyfil = 'Seamount_slice.bty';

alpha0 = 15;   % cone angle
x_source = 3;
y_source = 0;

theta = 135;   % bearing line

%%
rmax   = 15;
nr = 500;

r = linspace( 0, rmax, nr );
x = x_source + r * cosd( theta );
y = y_source + r * sind( theta );


%%

R = sqrt( x.^2 + y.^2 );
z = R * tand( alpha0 ) * 1000;
Bathy = [ r' z' ];

%%
writebdry( btyfil, 'L', Bathy )

figure
plotbty( btyfil )