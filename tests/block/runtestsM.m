% Effect of a block (forward and back-scatter)
% mbp

% the rays:
figure
bellhopM blockB_ray
%figure; plotray blockB_ray
plotbty 'blockB_ray'   % superimpose a bathymetry plot

%print -depsc blockB_ray

%%

% TL: Geometric ray theory
bellhopM blockB_geo
plotshd( 'blockB_geo.shd.mat', 2, 1, 1 )

plotbty 'blockB_geo'   % superimpose a bathymetry plot

% TL: Gaussian beams
bellhopM blockB_gb
plotshd( 'blockB_gb.shd.mat', 2, 1, 2 )

plotbty 'blockB_gb'   % superimpose a bathymetry plot
