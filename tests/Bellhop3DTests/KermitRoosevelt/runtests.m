
global units
units = 'km';

% bellhop3d KR2d
% 
% figure
% plotshdpol( 'KR2d.shd', 0, 0, 30 )
% caxisrev( [ 60 90 ] )
% axis image
% axis( [ 0 25 -4 4 ] )
% 
% 
% print -djpeg KR2D
% print -depsc2 KR2D


%%
copyfile( 'KR.bty', 'KR3dCart.bty' )   % copy over the bathymetry file
bellhop3d KR3dCart

figure
plotshdpol( 'KR3dCart.shd', 0, 0, 100 )
caxisrev( [ 40 75 ] )
axis image
axis( [ 0 10 -4 4 ] )

delete KR3dCart.bty
%print -depsc2 KR3DGeoHat
%print -djpeg KR3DGeoHat

%%
copyfile( 'KR.bty', 'KR3dGaussian.bty' )   % copy over the bathymetry file
bellhop3d KR3dGaussian

figure
plotshdpol( 'KR3dGaussian.shd', 0, 0, 100 )
caxisrev( [ 40 75 ] )
axis image
axis( [ 0 10 -4 4 ] )

delete KR3dGaussian.bty

%print -depsc2 KR3DGeoGaussian
%print -djpeg KR3DGeoGaussian
%%

% ray trace
copyfile( 'KR.bty', 'KR3d_ray.bty' )   % copy over the bathymetry file
bellhop3d KR3d_ray     % run BELLHOP3D on the KR3d test case

figure
plotray3d KR3d_ray.ray

delete KR3d_ray.bty
