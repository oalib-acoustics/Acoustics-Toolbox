% Penetrable Wedge test case

global units
units = 'km';

%% GeoHat beams, Ray centered coord.
bellhop3d pwedge3dRayCen

% polar plot of the TL
figure
plotshdpol( 'pwedge3dRayCen.shd', 0, -20, 100 )
caxisrev( [ 80 100 ] )
axis image
axis( [ 0 70 -50 0 ] )

print -dpng pwedge3dRayCen
%% GeoHat beams, Cart. coord.
bellhop3d pwedge3dCart

% polar plot of the TL
figure
plotshdpol( 'pwedge3dCart.shd', 0, -20, 100 )
caxisrev( [ 80 100 ] )
axis image
axis( [ 0 70 -50 0 ] )

print -dpng pwedge3dCart

%% GeoHat beams 2D
bellhop3d pwedge2d

% polar plot of the TL
figure
plotshdpol( 'pwedge2d.shd', 0, -20, 100 )
caxisrev( [ 80 100 ] )
axis image
axis( [ 0 70 -50 0 ] )

print -dpng pwedge2d

%% GeoHat beams 2D rotated
bellhop3d pwedge2d_rot

% polar plot of the TL
figure
plotshdpol( 'pwedge2d_rot.shd', -20, 0, 100 )
caxisrev( [ 80 100 ] )
axis image
axis( [ -50 0 -70 0 ] )

print -dpng pwedge2d_rot

%% GeoHat beams 3D rotated
bellhop3d pwedge3d_rot

% polar plot of the TL
figure
plotshdpol( 'pwedge3d_rot.shd', -20, 0, 100 )
caxisrev( [ 80 100 ] )
axis image
axis( [ -50 0 -70 0 ] )

print -dpng pwedge3d_rot

%% Side view
bellhop3d slice2d

% polar plot of the TL
figure
plotshd( 'slice2d.shd' )
caxisrev( [ 60 100 ] )
% axis( [ -50 0 -70 0 ] )
