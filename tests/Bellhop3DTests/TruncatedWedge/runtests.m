
global units
units = 'km';

makebty

%%
copyfile( 'wedge.bty', 'wedge2d.bty' );
bellhop3d wedge2d

figure
plotshdpol( 'wedge2d.shd', 0, 0, 30 )
caxisrev( [ 60 90 ] )
axis image
axis( [ 0 25 -4 4 ] )

delete( 'wedge2d.bty' )

print -depsc2 TWedge2D

%%
copyfile( 'wedge.bty', 'wedge3dHatRayCen.bty' );
bellhop3d wedge3dHatRayCen

%%
figure
plotshdpol( 'wedge3dHatRayCen.shd', 0, 0, 30 )
caxisrev( [ 60 90 ] )
axis image
axis( [ 0 25 -4 4 ] )

print -depsc2 TWedge3DHatRayCen

delete( 'wedge3dHatRayCen.bty' )

%%

copyfile( 'wedge.bty', 'wedge3dCart.bty' );
bellhop3d wedge3dCart

figure
plotshdpol( 'wedge3dCart.shd', 0, 0, 30 )
caxisrev( [ 60 90 ] )
axis image
axis( [ 0 25 -4 4 ] )

delete( 'wedge3dCart.bty' )

print -depsc2 TWedge3DHatCart

%%
copyfile( 'wedge.bty', 'wedge3dGaussian.bty' );
bellhop3d wedge3dGaussian

figure
plotshdpol( 'wedge3dGaussian.shd', 0, 0, 30 )
caxisrev( [ 60 90 ] )
axis image
axis( [ 0 25 -4 4 ] )

%delete( 'wedge3dGaussian.bty' )

print -dpng TWedge3DGeoGaussian


% set( gcf, 'Units', 'centimeters' )
% pos = get( gcf, 'Position' );
% set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [ pos(3), pos(4) ] )
% print -dpdf Twedge3dGaussian

%%
copyfile( 'wedge.bty', 'wedge3dGaussian_side.bty' );
bellhop3d wedge3dGaussian_side

figure
plotshd( 'wedge3dGaussian_side.shd' )
caxisrev( [ 60 90 ] )
%%
copyfile( 'wedge.bty', 'wedge3dHat_side.bty' );
bellhop3d wedge3dHat_side

figure
plotshd( 'wedge3dHat_side.shd' )
caxisrev( [ 60 90 ] )

%%

% ray trace
copyfile( 'wedge.bty', 'wedge3d_ray.bty' )   % copy over the bathymetry file
bellhop3d wedge3d_ray     % run BELLHOP3D on the wedge3d test case

figure
plotbdry3d wedge.bty
hold on
plotray3d wedge3d_ray.ray

delete( 'wedge3d_ray.bty' )
delete( 'wedge.bty' )
