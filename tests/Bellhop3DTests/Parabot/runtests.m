% run the ParaBot test case
% p517 computational acoustics

% note that beams traced with negative azimuthal angles can exit the box
% where the boundaries are defined

global units
units = 'km';
%%
makebdry              % make the bathymetry

figure
plotbdry3d ParaBot.bty
shading flat

% ray trace
copyfile( 'ParaBot.bty', 'ParaBot_ray.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot_ray.ati' )   % copy over the bathymetry file

bellhop3d ParaBot_ray

hold on
plotray3d ParaBot_ray.ray

delete( 'ParaBot_ray.bty' )
delete( 'ParaBot_ray.ati' )

%%
copyfile( 'ParaBot.bty', 'ParaBot2D.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot2D.ati' )   % copy over the bathymetry file

bellhop3d ParaBot2D

% plot of the TL
figure
plotshd( 'ParaBot2D.shd' )
caxisrev( [ 60 100 ] )

delete( 'ParaBot2D.bty' )
delete( 'ParaBot2D.ati' )

%%
% 3d run polar (GeoHat Cartesian)
% looking for a perfect cylindrical symmetry without any ripple
% You get that with enough beams and using the analytic formulas for the
% curvature
% The curvilinear interpolation does not produce that

copyfile( 'ParaBot.bty', 'ParaBot3DHatCartPolar.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DHatCartPolar.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DHatCartPolar

% polar plot of the TL
figure
plotshdpol( 'ParaBot3DHatCartPolar.shd', 0, 0, 1000 )
caxisrev( [ 65 75 ] )

delete( 'ParaBot3DHatCartPolar.bty' )
delete( 'ParaBot3DHatCartPolar.ati' )
%%
% 3d run (GeoHat Cartesian)

copyfile( 'ParaBot.bty', 'ParaBot3DHatCart.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DHatCart.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DHatCart

% polar plot of the TL
figure
plotshd( 'ParaBot3DHatCart.shd' )
caxisrev( [ 60 100 ] )

delete( 'ParaBot3DHatCart.bty' )
delete( 'ParaBot3DHatCart.ati' )
%%
% 3d run, 90 degree bearing (GeoHat Cartesian)

copyfile( 'ParaBot.bty', 'ParaBot3DHatCart_theta90.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DHatCart_theta90.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DHatCart_theta90

% plot of the TL
figure
plotshd( 'ParaBot3DHatCart_theta90.shd' )
caxisrev( [ 60 100 ] )

delete( 'ParaBot3DHatCart_theta90.bty' )
delete( 'ParaBot3DHatCart_theta90.ati' )
%%
% 3d run (GeoHat Ray centered)

copyfile( 'ParaBot.bty', 'ParaBot3DHatRaycen.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DHatRaycen.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DHatRaycen

% plot of the TL
figure
plotshd( 'ParaBot3DHatRaycen.shd' )
caxisrev( [ 60 100 ] )

delete( 'ParaBot3DHatRaycen.bty' )
delete( 'ParaBot3DHatRaycen.ati' )


%%
% 3d run (GeoGaussian)

copyfile( 'ParaBot.bty', 'ParaBot3DGaussian.bty' )   % copy over the bathymetry file
copyfile( 'ParaBot.ati', 'ParaBot3DGaussian.ati' )   % copy over the bathymetry file

bellhop3d ParaBot3DGaussian

% plot of the TL
figure
plotshd( 'ParaBot3DGaussian.shd' )
caxisrev( [ 60 100 ] )

delete( 'ParaBot3DGaussian.bty' )
delete( 'ParaBot3DGaussian.ati' )

%delete( 'ParaBot.bty' )
%delete( 'ParaBot.ati' )
%print -depsc2 ParaBot3DGaussian
%print -djpeg ParaBot3DGaussian

delete( 'ParaBot3D.bty' )
delete( 'ParaBot3D.ati' )

%%
delete( 'ParaBot.bty' )
delete( 'ParaBot.ati' )
