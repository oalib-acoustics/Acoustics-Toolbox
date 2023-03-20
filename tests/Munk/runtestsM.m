% Munk profile test cases
% mbp

figure
bellhopM( 'MunkB_ray' )
%figure
%plotray( 'MunkB_ray' )

%bellhop( 'MunkB_eigenray' )
%figure
%plotray( 'MunkB_eigenray' )

bellhopM( 'MunkB_Coh' )
plotshd( 'MunkB_Coh.shd.mat', 2, 2, 1 )
caxisrev( [ 50 100 ] )

%%
bellhopM( 'MunkB_gb' )
plotshd( 'MunkB_gb.shd.mat', 2, 2, 2 )
caxisrev( [ 50 100 ] )

kraken( 'MunkK' )
plotshd( 'MunkK.shd.mat', 2, 2, 3 )
caxisrev( [ 50 100 ] )

scooterM( 'MunkS' )
plotshd( 'MunkS.shd.mat', 2, 2, 4 )
caxisrev( [ 50 100 ] )

%%
% tests of incoherent, semi-coherent options

bellhopM( 'MunkB_Coh' )
plotshd( 'MunkB_Coh.shd.mat', 3, 2, 1 )
caxisrev( [ 50 100 ] )

bellhopM( 'MunkB_Coh_gb' )
plotshd( 'MunkB_Coh_gb.shd.mat', 3, 2, 2 )
caxisrev( [ 50 100 ] )

bellhopM( 'MunkB_Semi' )
plotshd( 'MunkB_Semi.shd.mat', 3, 2, 3 )
caxisrev( [ 50 100 ] )

bellhopM( 'MunkB_Semi_gb' )
plotshd( 'MunkB_Semi_gb.shd.mat', 3, 2, 4 )
caxisrev( [ 50 100 ] )

bellhopM( 'MunkB_Inc' )
plotshd( 'MunkB_Inc.shd.mat', 3, 2, 5 )
caxisrev( [ 50 100 ] )

bellhopM( 'MunkB_Inc_gb' )
plotshd( 'MunkB_Inc_gb.shd.mat', 3, 2, 6 )
caxisrev( [ 50 100 ] )

% test of Green's function plotting
%figure;
%plotgrn( 'MunkS.grn' )
