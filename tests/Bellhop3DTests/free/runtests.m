% Run tests to verify the surface and bottom reflection coefficents are
% incorporated correctly

% free space

global units
units = 'km';

scooter( 'freeSPoint' )
plotshd( 'freeSPoint.shd.mat', 5, 1, 1 );
caxisrev( [ 60 80 ] )

bellhop3d( 'freeBhat' )
plotshd( 'freeBhat.shd', 5, 1, 2 );
caxisrev( [ 60 80 ] )

bellhop3d( 'freeBhat_raycen' )
plotshd( 'freeBhat_raycen.shd', 5, 1, 3 );
caxisrev( [ 60 80 ] )

bellhop3d( 'freeBgaussian' )
plotshd( 'freeBgaussian.shd', 5, 1, 4 );
caxisrev( [ 60 80 ] )

bellhop3d( 'freeBgaussian_raycen' )
plotshd( 'freeBgaussian_raycen.shd', 5, 1, 5 );
caxisrev( [ 60 80 ] )

%%
% polar plots

bellhop3d( 'freeBhatpolar' )
figure
plotshdpol( 'freeBhatpolar.shd', 0.0, 0.0, 3000 );
caxisrev( [ 60 80 ] )

%%

bellhop3d( 'freeBgaussianpolar' )
figure
plotshdpol( 'freeBgaussianpolar.shd', 0.0, 0.0, 3000 );
caxisrev( [ 60 80 ] )
