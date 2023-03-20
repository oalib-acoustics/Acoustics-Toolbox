% Surface duct test cases
% mbp

global units jkpsflag

units = 'km';

%%
bellhop( 'sductB' )
plotshd( 'sductB.shd', 2, 2, 1 )
caxisrev( [ 50 100 ] )

%%
bellhop( 'sductB_gb' )
plotshd( 'sductB_gb.shd', 2, 2, 2 )
caxisrev( [ 50 100 ] )

%%
kraken( 'sductK' )
plotshd( 'sductK.shd.mat', 2, 2, 3 )
caxisrev( [ 50 100 ] )

%%
scooter( 'sductS' )
plotshd( 'sductS.shd.mat', 2, 2, 4 )
caxisrev( [ 50 100 ] )

%%
simplePE sductB

figure
plotshd( 'sductB.shd.mat' )
caxisrev( [ 50 100 ] )
