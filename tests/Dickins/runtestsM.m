% Dickins seamount test
% mbp

global units jkpsflag

units = 'km';

bellhopM( 'DickinsB' )
plotshd( 'DickinsB.shd.mat', 2, 2, 1 )
caxisrev( [ 70 120 ] )
% superimpose a bathymetry plot
plotbty 'DickinsB'

bellhopM( 'DickinsFlatB' )
plotshd( 'DickinsFlatB.shd.mat', 2, 2, 2 )
caxisrev( [ 70 120 ] )

kraken( 'DickinsFlatK' )
plotshd( 'DickinsFlatK.shd.mat', 2, 2, 3 )
caxisrev( [ 70 120 ] )

scooterM( 'DickinsFlatS' )
plotshd( 'DickinsFlatS.shd.mat', 2, 2, 4 )
caxisrev( [ 70 120 ] )

