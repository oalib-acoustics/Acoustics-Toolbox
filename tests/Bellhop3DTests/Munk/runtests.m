% run the Munk3D test cases
% mbp, June 2011

% note: the earlier version of this test case had the seamount in a
% different place

global units
units = 'km';

%%
makebdry             % make the bathymetry

figure
plotbdry3d munk3d.bty

%%
bellhop3d munk3d     % run BELLHOP3D on the munk3d test case

% polar plot of the TL
figure
plotshdpol( 'munk3d.shd', 0, 0, 1000 )
caxisrev( [ 60 120 ] )

%%
% ray trace
copyfile( 'munk3d.bty', 'munk3d_ray.bty' )   % copy over the bathymetry file
bellhop3d munk3d_ray     % run BELLHOP3D on the munk3d test case

figure
plotray3d munk3d_ray.ray
%%
% sideviews
% the Gaussian, Nx2D is smoothed because of the stent
% the Gaussian  3D   has the stent turned off

bellhop3d slice2dHatCart
plotshd( 'slice2dHatCart.shd', 2, 2, 1 )
caxisrev( [ 50 100 ] )

bellhop3d slice3dHatCart
plotshd( 'slice3dHatCart.shd', 2, 2, 2 )
caxisrev( [ 50 100 ] )

bellhop3d slice2dGaussian
plotshd( 'slice2dGaussian.shd', 2, 2, 3 )
caxisrev( [ 50 100 ] )

bellhop3d slice3dGaussian
plotshd( 'slice3dGaussian.shd', 2, 2, 4 )
caxisrev( [ 50 100 ] )

%%
% Tests of different coherence options
% Hat in Ray-centered coordinates

% Hat

bellhop3d slice3dHatRaycen
plotshd( 'slice3dHatRaycen.shd', 3, 1, 1 )
caxisrev( [ 50 100 ] )

bellhop3d slice3dHatRaycenSemi
plotshd( 'slice3dHatRaycenSemi.shd', 3, 1, 2 )
caxisrev( [ 50 100 ] )

bellhop3d slice3dHatRaycenInc
plotshd( 'slice3dHatRaycenInc.shd', 3, 1, 3 )
caxisrev( [ 50 100 ] )

%%
% Hat in Cartesian coordinates

bellhop3d slice3dHatCart   % this is redundant if you ran the whole script
plotshd( 'slice3dHatCart.shd', 3, 1, 1 )
caxisrev( [ 50 100 ] )

bellhop3d slice3dHatCartSemi
plotshd( 'slice3dHatCartSemi.shd', 3, 1, 2 )
caxisrev( [ 50 100 ] )

bellhop3d slice3dHatCartInc
plotshd( 'slice3dHatCartInc.shd', 3, 1, 3 )
caxisrev( [ 50 100 ] )

%%
% Gaussian
bellhop3d slice3dGaussian   % this is redundant if you ran the whole script
plotshd( 'slice3dGaussian.shd', 3, 1, 1 )
caxisrev( [ 50 100 ] )

bellhop3d slice3dGaussianSemi
plotshd( 'slice3dGaussianSemi.shd', 3, 1, 2 )
caxisrev( [ 50 100 ] )

bellhop3d slice3dGaussianInc
plotshd( 'slice3dGaussianInc.shd', 3, 1, 3 )
caxisrev( [ 50 100 ] )
