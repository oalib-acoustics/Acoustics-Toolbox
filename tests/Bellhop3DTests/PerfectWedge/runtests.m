% run the Wedge3D test case
% This is for the 'pefect' wedge with vacumm boundary conditions all around
% xaz, Jan 2012
% need to go beyong -90 to 90 degrees in azimuth to get the field on the
% y-axis

global units
units = 'km';

%%
makebty             % make the bathymetry
axis( [ 0 20 -40 -5 0 1000 ] )

%%
% Horizontal refraction disabled
% You get a horible mess in the halfspace above the source because of the
% interference between the energy going up the wedge and then getting fully
% reversed. This is a poor example in the sense that you wouldn't generally
% see that.

% There are further artifacts in shallower water since the receiver is then
% close to the bottom.

% The Gaussian beam option uses very fat beams because the 10-Hz
% source has a 150-m wavelength. The beams fill the water column in many
% areas. Stints need to be turned off inside the Fortran code

% The rays or beams are reversed when going upslope. Because of the
% cylindrical symmetry assumption, the beams refocus back at the source.
% A more realistic option would be to use a translationally invariant
% environment.

copyfile( 'wedge2d.bty', 'wedge2dray.bty' )
bellhop3d wedge2dray

hold on
plotray3d( 'wedge2dray.ray' )

%%
% With the 'g' option (geometric beams in ray-centered coordinates),
% you get several bands (artifacts) in the upper halfplane.
% These are zones near mode cutoff, i.e. where the rays go vertical
% just before getting dropped.

% With the 'G' option (geometric beams in cartesian coordinates),
% You get a picture with lots of static in the UHP due to rays
% that are reflected back to the source

bellhop3d wedge2d     % run BELLHOP3D on the wedge2d test case

% polar plot of the TL
figure
plotshdpol( 'wedge2d.shd', 0, -19.1, 80 )
caxisrev( [ 60 85 ] )
axis image
axis( [ 0 20 -40 -5 ] )

print -dpng   wedge2d

%%
% Hat-shaped beams
copyfile( 'wedge2d.bty', 'wedge3dhat.bty' )
bellhop3d wedge3dHat     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'wedge3dHat.shd', 0, -19.1, 80 )
caxisrev( [ 60 85 ] )
axis image
axis( [ 0 20 -40 -5 ] )

print -dpng   wedge3dHat

%%
% Hat-shaped beams, ray-centered
copyfile( 'wedge2d.bty', 'wedge3dHatRayCen.bty' )
bellhop3d wedge3dHatRayCen     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'wedge3dHatRayCen.shd', 0, -19.1, 80 )
caxisrev( [ 60 85 ] )
axis image
axis( [ 0 20 -40 -5 ] )

print -dpng   wedge3dHatRayCen

%%
% Gaussian beams

copyfile( 'wedge2d.bty', 'wedge3dGaussian.bty' )
bellhop3d wedge3dGaussian     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'wedge3dGaussian.shd', 0, -19.1, 80 )
caxisrev( [ 60 85 ] )
axis image
axis( [ 0 20 -40 -5 ] )

print -dpng wedge3dGaussian

%%
% Gaussian beams, Ray Cen

copyfile( 'wedge2d.bty', 'wedge3dGaussianRayCen.bty' )
bellhop3d wedge3dGaussianRayCen     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'wedge3dGaussianRayCen.shd', 0, -19.1, 80 )
caxisrev( [ 60 85 ] )
axis image
axis( [ 0 20 -40 -5 ] )

print -dpng wedge3dGaussianRayCen

% set( gca, 'ActivePositionProperty', 'Position', 'Units', 'centimeters' )
% set( gcf, 'Units', 'centimeters' )
% set( gcf, 'PaperPositionMode', 'auto');   % this is important; default is 6x8 inch page
%set( gca, 'Position', [ 2    2                       14.0       7.0 ] )
%set( gcf, 'Units', 'centimeters' )
%set( gcf, 'Position', [ 3 15 19.0 10.0 ] )
%set( gcf, 'PaperSize', [ 6 6 ] )

% set( gcf, 'Units', 'centimeters' )
% pos = get( gcf, 'Position' );
% set( gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [ pos(3), pos(4) ] )
% print -dpdf Perfectwedge3dGaussian


%%

% ray trace
copyfile( 'wedge2d.bty', 'wedge3d_ray.bty' )   % copy over the bathymetry file
bellhop3d wedge3d_ray     % run BELLHOP3D on the wedge3d test case

figure
plotbdry3d wedge3d_ray.bty
hold on
plotray3d wedge3d_ray.ray
