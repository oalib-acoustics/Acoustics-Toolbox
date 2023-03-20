% run the HyperBot test case
% p517 computational acoustics

global units
units = 'km';
%%
makebty              % make the bathymetry
figure
plotbdry3d HyperBot
shading flat

% ray trace
copyfile( 'HyperBot.bty', 'HyperBot_ray.bty' )   % copy over the bathymetry file

bellhop3d HyperBot_ray     % run BELLHOP3D on the wedge3d test case

hold on
plotray3d HyperBot_ray.ray

delete( 'HyperBot_ray.bty' )

%%
copyfile( 'HyperBot.bty', 'HyperBot2D.bty' )   % copy over the bathymetry file

bellhop3d HyperBot2D     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'HyperBot2D.shd' )
caxisrev( [ 60 100 ] )

delete( 'HyperBot2D.bty' )

%%
% 3d run (GeoHat Cartesian)

copyfile( 'HyperBot.bty', 'HyperBot3DHatCart.bty' )   % copy over the bathymetry file

bellhop3d HyperBot3DHatCart     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshd( 'HyperBot3DHatCart.shd' )
caxisrev( [ 60 100 ] )

delete( 'HyperBot3DHatCart.bty' )

%%
% 3d run (GeoHat Cartesian)

copyfile( 'HyperBot.bty', 'HyperBot3DPolar.bty' )   % copy over the bathymetry file

bellhop3d HyperBot3DPolar     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol( 'HyperBot3DPolar.shd', 0, 0, 500 )
caxisrev( [ 60 100 ] )

delete( 'HyperBot3DPolar.bty' )
