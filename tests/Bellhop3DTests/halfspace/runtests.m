% Run tests to verify the surface and bottom reflection coefficents are
% incorporated correctly

% note that upper and lower cases are not symmetric because the source is
% not in the middle

% halfspaces

scooter( 'lower_halfS' )
plotshd( 'lower_halfS.shd.mat', 2, 2, 1 );
caxisrev( [ 60 80 ] )

bellhop3d( 'lower_halfB' )
plotshd( 'lower_halfB.shd', 2, 2, 2 );
caxisrev( [ 60 80 ] )

scooter( 'upper_halfS' )
plotshd( 'upper_halfS.shd.mat', 2, 2, 3 );
caxisrev( [ 60 80 ] )

bellhop3d( 'upper_halfB' )
plotshd( 'upper_halfB.shd', 2, 2, 4 );
caxisrev( [ 60 80 ] )

print -depsc Half