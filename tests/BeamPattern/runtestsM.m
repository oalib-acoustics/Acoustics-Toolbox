% Run tests to verify the source beam pattern is incorporated correctly in BELLHOP

global units
units = 'm';

bellhopM( 'omni' )
plotshd( 'omni.shd.mat', 3, 1, 1 )
caxisrev( [ 50 100 ] )
axis( 'image' )

bellhopM( 'shaded' )
plotshd( 'shaded.shd.mat', 3, 1, 2 )
caxisrev( [ 50 100 ] )
axis( 'image' )

scooterM shadedS
plotshd( 'shadedS.shd.mat', 3, 1, 3 )
caxisrev( [ 50 100 ] )
axis( 'image' )

% Munk with beam pattern
% scooter runs

scooterM MunkS
plotshd( 'MunkS.shd.mat', 2, 1, 1 )
caxisrev( [ 50 100 ] )

% kraken runs

krakenc MunkKleaky
%field MunkKleaky
plotshd( 'MunkKleaky.shd.mat', 2, 1, 2 )
caxisrev( [ 50 100 ] )
