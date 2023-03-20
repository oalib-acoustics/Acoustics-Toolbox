% Run tests to verify the capability to use a source beam pattern in
% BELLHOP
global units
units = 'm';

bellhop( 'omni' )
plotshd( 'omni.shd', 3, 1, 1 );
caxisrev( [ 50 100 ] )
axis( 'image' )

bellhop( 'shaded' )
plotshd( 'shaded.shd', 3, 1, 2 );
caxisrev( [ 50 100 ] )
axis( 'image' )

scooter shadedS
plotshd( 'shadedS.shd.mat', 3, 1, 3 )
caxisrev( [ 50 100 ] )
axis( 'image' )

%%
% Munk case

scooter MunkS
plotshd( 'MunkS.shd.mat', 2, 1, 1 )
caxisrev( [ 50 100 ] )

% kraken runs

krakenc MunkKleaky
%field MunkKleaky
plotshd( 'MunkKleaky.shd.mat', 2, 1, 2 )
caxisrev( [ 50 100 ] )
