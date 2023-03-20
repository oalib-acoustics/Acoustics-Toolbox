% Munk profile test cases
% mbp

global units jkpsflag
units = 'km';

%%
bellhop( 'arcticB' )
plotshd( 'arcticB.shd', 2, 2, 1 )
caxisrev( [ 50 100 ] )

bellhop( 'arcticB_gb' )
plotshd( 'arcticB_gb.shd', 2, 2, 2 )
caxisrev( [ 50 100 ] )

kraken( 'arcticK' )
plotshd( 'arcticK.shd.mat', 2, 2, 3 )
caxisrev( [ 50 100 ] )

scooter( 'arcticS' )
plotshd( 'arcticS.shd.mat', 2, 2, 4 )
caxisrev( [ 50 100 ] )

simplePE arcticB
figure
plotshd( 'arcticB.shd.mat' )
caxisrev( [ 50 100 ] )
