% Run tests to verify the surface and bottom reflection coefficients are
% incorporated correctly. These tests compare BELLHOP3D in full 3D mode
% to BELLHOP.

bellhop3d( 'lower_halfB3' )
plotshd( 'lower_halfB3.shd', 2, 2, 1 );
caxisrev( [ 50 80 ] )

bellhop( 'lower_halfB' )
plotshd( 'lower_halfB.shd', 2, 2, 2 );
caxisrev( [ 50 80 ] )

bellhop3d( 'upper_halfB3' )
plotshd( 'upper_halfB3.shd', 2, 2, 3 );
caxisrev( [ 50 80 ] )

bellhop( 'upper_halfB' )
plotshd( 'upper_halfB.shd', 2, 2, 4 );
caxisrev( [ 50 80 ] )

% Run tests to verify the reflection coefficients are incorporated
% correctly for the bottom and a nearly vertical wall. This test compares
% BELLHOP3D in full 3D mode to BELLHOP. There are some small differences
% at ranges less than 1.5 km or so. This is to be expected as the BELLHOP3D
% wall is flat, while the BELLHOP wall is actually cylindrical which
% produces some focusing of the reflected rays.

bellhop3d( 'lower_half_wallB3' )
plotshd( 'lower_half_wallB3.shd', 1, 2, 1 );
caxisrev( [ 50 80 ] )

bellhop( 'lower_half_wallB' )
plotshd( 'lower_half_wallB.shd', 1, 2, 2 );
caxisrev( [ 50 80 ] )

%
% End of runtests.m
