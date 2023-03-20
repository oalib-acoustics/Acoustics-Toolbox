% run the seamount test case
% p517 computational acoustics
% You must use the curvilinear fit to the bathymetry to get a smooth plot
% The piecewise linear fit produces the disco-ball effect with a very
% ragged TL

% For a receiver depth at 100 m, you also get a noise picture even with the
% curvilinear bottom fit. This is the usual BELLHOP issue when a receiver
% gets close to a boundary because of BELLHOP does not compute an influence
% from both the incident beam and its mirror image at the same time.

global units
units = 'km';
%%
makebty              % make the bathymetry

figure
plotbdry3d Seamount.bty
shading flat

%% ray trace
copyfile( 'Seamount.bty', 'Seamount_ray.bty' )   % copy over the bathymetry file
bellhop3d Seamount_ray

hold on
plotray3d Seamount_ray.ray

delete Seamount_ray.bty

%%

% when I used 101 beams in declination angle over +/- 90 degrees there was
% a discontinuity in the field at a circle of radius 3.1 km or so. This is
% presumably because the beams are wide and the receiver is near the
% surface.
% This disappears with 1201 beams

copyfile( 'Seamount.bty', 'Seamount2D.bty' )   % copy over the bathymetry file

bellhop3d Seamount2D

% polar plot of the TL
figure
plotshdpol( 'Seamount2D.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

delete Seamount2D.bty

%%
% 3d run (GeoHat Cartesian)

copyfile( 'Seamount.bty', 'Seamount3DHatcart.bty' )   % copy over the bathymetry file

bellhop3d Seamount3DHatcart

% polar plot of the TL
figure
plotshdpol( 'Seamount3DHatcart.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

delete Seamount3DHatcart.bty
%%
% 3d run (GeoHat Ray centered)

copyfile( 'Seamount.bty', 'Seamount3DHatRaycen.bty' )   % copy over the bathymetry file

bellhop3d Seamount3DHatRaycen

% polar plot of the TL
figure
plotshdpol( 'Seamount3DHatRaycen.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

delete Seamount3DHatRaycen.bty

%%
% 3d run (GeoGaussian)

'This is way over-sampled'
'Can reduce rays, receivers by a factor of 2 in both directions'

copyfile( 'Seamount.bty', 'Seamount3DGaussian.bty' )   % copy over the bathymetry file

bellhop3d Seamount3DGaussian

% polar plot of the TL
figure
plotshdpol( 'Seamount3DGaussian.shd', 3, 0, 100 )
caxisrev( [ 40 100 ] )
axis( [ -6 3 0 9 ] )

%print -dpng Seamount3DGaussian

% set( gcf, 'Units', 'centimeters' )
% pos = get( gcf, 'Position' );
% set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [ pos(3), pos(4) ] )
% print -dpdf Seamount3DGaussian

delete Seamount3DGaussian.bty
%%
% sideview on the 135 degree bearing from BELLHOP3D

bellhop3d Seamount2D_side

figure
plotshd( 'Seamount2D_side.shd' )
caxisrev( [ 40 100 ] )

hold on
plotbty Seamount_slice

%%
% sideview on the 135 degree bearing from BELLHOP

bellhop Seamount_slice

figure
plotshd( 'Seamount_slice.shd' )
caxisrev( [ 40 100 ] )

hold on
plotbty Seamount_slice


