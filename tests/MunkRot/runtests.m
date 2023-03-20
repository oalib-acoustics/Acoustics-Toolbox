% run the 3D Munk case

global units
units = 'km';

bellhop Munk

figure
plotshd( 'Munk.shd' )
caxisrev( [ 50 100 ] )

%%
bellhop MunkRot

units = 'm';

figure
plotshd( 'MunkRot.shd' )
caxisrev( [ 50 100 ] )
