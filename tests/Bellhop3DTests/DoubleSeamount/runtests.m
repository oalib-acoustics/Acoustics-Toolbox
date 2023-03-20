% run the double seamount test case (created by YT Lin)

global units
units = 'km';

%%
makebty

figure
plotbdry3d( btyfil )
axis( [ 0 6 -4 4 ] )

%print -depsc2 Bathy

% ray trace

bellhop3d DoubleSeamount3D_ray

hold on
plotray3d( 'DoubleSeamount3D_ray.ray' )
axis( [ 0 6 -4 4 0 900 ] )
print -depsc2 DoubleSeamount3D_ray

%%
% default beam spacing becomes about 300 beams in each direction

bellhop3d DoubleSeamount3D


%%
figure
plotshdpol( 'DoubleSeamount3D.shd',  0, 0, 400 )
caxis( [ 40 80 ] )
axis equal
axis( [ 0 6 -4 4 ] )

print -depsc2 3DGaussian
%%
% default beam spacing becomes about 300 beams in each direction

copyfile( 'DoubleSeamount3D.bty', 'DoubleSeamount3DHat.bty' )

bellhop3d DoubleSeamount3DHat

%%
figure
plotshdpol( 'DoubleSeamount3DHat.shd',  0, 0, 400 )
caxis( [ 40 80 ] )
axis equal
axis( [ 0 6 -4 4 ] )

print -dpng 3DHat

%%
delete DoubleSeamount3D.bty
delete DoubleSeamount3DHat.bty
