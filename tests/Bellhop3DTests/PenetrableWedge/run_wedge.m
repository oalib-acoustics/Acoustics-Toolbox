% run the Wedge3D test case
% xaz, Jan 2012

global units
units = 'km';

makebty             % make the bathymetry

bellhop3d Wedge3d     % run BELLHOP3D on the wedge3d test case

% polar plot of the TL
figure
plotshdpol 'Wedge3d_2.shd'
%caxisrev( [ 60 120 ] )

if 0
    % ray trace
    copyfile( 'Wedge3d.bty', 'Wedge3d_ray.bty' )   % copy over the bathymetry file
    bellhop3d Wedge3d_ray     % run BELLHOP3D on the wedge3d test case
    
    figure
    plotray3d Wedge3d_ray.ray
end
