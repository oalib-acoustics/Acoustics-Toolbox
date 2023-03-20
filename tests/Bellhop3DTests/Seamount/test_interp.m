makebty_slice

figure( 10 )
plot( r, z )
hold on

xq = x;
yq = y;
rint = r;

figure
plotbdry3d( 'Seamount2D_side.bty' )

load   % needs a save statement added to plotbty3d
zint = interp2( X, Y, z, xq'*1000, yq'*1000 );

figure( 10 )
plot( rint, zint, 'r' )
