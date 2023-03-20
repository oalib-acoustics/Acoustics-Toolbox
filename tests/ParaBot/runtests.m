% Parabolic Bottom test
% mbp

% Ray 118 hits the bottom extremely close to a node
% (Same thing happens for another ray hitting the surface)
% This allowed the ray to exit the domain
% A adjustment was made to BELLHOP to reduce the overstep from
% 1e-3 * deltas to
% 1e-4 * deltas, which fixed the problem in this case
% However, this is a fundmental flaw that can still happen in other cases

global units jkpsflag
units = 'km';

% *************************************************
% linear boundary interpolation

make_bdry( 'L' )

%%
% the rays:
bellhop ParaBot
figure; plotray ParaBot
axis( [ 0 20 -5000 5000 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Geometric ray theory
copyfile( 'ParaBot.bty', 'ParaBotTLGeom.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGeom.ati' );

bellhop ParaBotTLGeom

plotshd( 'ParaBotTLGeom.shd', 2, 1, 1 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

%%
% TL: Gaussian beams
copyfile( 'ParaBot.bty', 'ParaBotTLGB.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGB.ati' );

bellhop ParaBotTLGB

plotshd( 'ParaBotTLGB.shd', 2, 1, 2 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

%print -dpng ParabotLTL
%%

% *************************************************
% curvilinear boundary interpolation

make_bdry( 'C' )

%%
% the rays:
bellhop ParaBot
figure; plotray ParaBot
axis( [ 0 20 -5000 5000 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Geometric ray theory
copyfile( 'ParaBot.bty', 'ParaBotTLGeom.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGeom.ati' );
bellhop ParaBotTLGeom
plotshd( 'ParaBotTLGeom.shd', 2, 1, 1 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot

% TL: Gaussian beams
copyfile( 'ParaBot.bty', 'ParaBotTLGB.bty' );
copyfile( 'ParaBot.ati', 'ParaBotTLGB.ati' );
bellhop ParaBotTLGB
plotshd( 'ParaBotTLGB.shd', 2, 1, 2 )
caxisrev( [ 60 100 ] )

plotati 'ParaBot'   % superimpose an altimetry plot
plotbty 'ParaBot'   % superimpose a bathymetry plot
