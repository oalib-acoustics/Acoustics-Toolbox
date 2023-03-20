% Parabolic Bottom test
% mbp

global units jkpsflag
units = 'km';

% linear boundary interpolation

make_bdry( 'L' )

%%
% the rays:
figure
bellhopM ParaBot
%figure; plotray RAYFIL
hold on
plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot
axis( [ 0 20000 -5000 5000 ] )

% TL: Geometric ray theory
copyfile( 'ParaBot.bty', 'ParaBotTLGeom.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGeom.ati' );
bellhopM ParaBotTLGeom
plotshd( 'ParaBotTLGeom.shd.mat', 2, 1, 1 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

%%
% TL: Gaussian beams
copyfile( 'ParaBot.bty', 'ParaBotTLGB.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGB.ati' );
bellhopM ParaBotTLGB
plotshd( 'ParaBotTLGB.shd.mat', 2, 1, 2 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot
%%

% curvilinear boundary interpolation

make_bdry( 'C' )

%%
% the rays:
figure
bellhopM ParaBot
%figure; plotray RAYFIL
hold on
plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Geometric ray theory
copyfile( 'ParaBot.bty', 'ParaBotTLGeom.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGeom.ati' );
bellhopM ParaBotTLGeom
plotshd( 'ParaBotTLGeom.shd.mat', 2, 1, 1 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

%%
% TL: Gaussian beams
copyfile( 'ParaBot.bty', 'ParaBotTLGB.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGB.ati' );
bellhopM ParaBotTLGB
plotshd( 'ParaBotTLGB.shd.mat', 2, 1, 2 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

