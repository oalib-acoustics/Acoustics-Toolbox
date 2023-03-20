function runtests
% Dickins seamount test
% mbp
global units jkpsflag
units = 'km';

%%
bellhop( 'DickinsBray' )
figure
plotray( 'DickinsBray' )
plotbty 'DickinsB'

bellhop( 'DickinsB_oneBeam' )

figure
plotshd( 'DickinsB_oneBeam.shd' )
caxisrev( [ 70 120 ] )

plotbty 'DickinsB'

%%

bellhop( 'DickinsB' )
plotshd( 'DickinsB.shd', 2, 2, 1 )
caxisrev( [ 70 120 ] )

% superimpose a bathymetry plot
plotbty 'DickinsB'

bellhop( 'DickinsFlatB' )
plotshd( 'DickinsFlatB.shd', 2, 2, 2 )
caxisrev( [ 70 120 ] )

copyfile( 'field_ri.flp', 'field.flp' )
kraken(  'DickinsFlatK' )
plotshd( 'DickinsFlatK.shd.mat', 2, 2, 3 )
caxisrev( [ 70 120 ] )

scooter( 'DickinsFlatS' )
plotshd( 'DickinsFlatS.shd.mat', 2, 2, 4 )
caxisrev( [ 70 120 ] )

%%
% Coupled mode runs for range-dependent case

% run Dickins case, using the wedge ocean

% cd Precalc
% runtests
% cd ..
% 
% % do it a second way using a simple marching
% 

cd March
runtests
cd ..

% Note: can see differences between deltar = 1 and 10 m
% ram
% 
% plotshd( 'DickinsB.shd', 2, 1, 1 )
% caxisrev( [ 70 120 ] )
% plotbty 'DickinsB'
% 
% plotshd( 'RAM.shd.mat', 2, 1, 2 )
% caxisrev( [ 70 120 ] )
% plotbty 'DickinsB'

%%

plotshd( 'DickinsFlatS.shd.mat', 2, 1, 1 )
caxisrev( [ 70 120 ] )

simplePE DickinsFlatB

plotshd( 'DickinsFlatB.shd.mat', 2, 1, 2 )
caxisrev( [ 70 120 ] )

simplePE DickinsB

figure
plotshd( 'DickinsB.shd.mat' )
caxisrev( [ 70 120 ] )


%%
% plots for jkps

if ( jkpsflag )
    
    figure
    plotray( 'DickinsBray' )
    plotbty 'DickinsB'
    axis( [ 0 100 0 3000 ] )
    print -depsc2 DickinsBray.eps
    
    figure
    plotshd( 'DickinsB.shd' )
    caxisrev( [ 70 120 ] )
    
    plotbty 'DickinsB'
    print -depsc2 DickinsB
end