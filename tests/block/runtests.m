% Effect of a block (forward and back-scatter)
% mbp

% the rays:
bellhop blockB_ray
figure; plotray blockB_ray
plotbty 'blockB_ray'   % superimpose a bathymetry plot

%print -depsc blockB_ray

%%

% TL: Geometric ray theory
bellhop blockB_geo
plotshd( 'blockB_geo.shd', 2, 1, 1 )

plotbty 'blockB_geo'   % superimpose a bathymetry plot

% TL: Gaussian beams
bellhop blockB_gb
plotshd( 'blockB_gb.shd', 2, 1, 2 )

plotbty 'blockB_gb'   % superimpose a bathymetry plot

%%

% print
% figure
% plotshd( 'blockB_gb.shd' )
% plotbty 'blockB_gb'   % superimpose a bathymetry plot
%print -depsc blockB_gb
