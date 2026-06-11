% run the provinces3D test cases
% mbp, September 2024
global units
units = 'km';

%%
makebdry             % make the bathymetry

%%
figure
plotbdry3d provinces.bty

%%
bellhop3d provinces     % run BELLHOP3D on the provinces3d test case

% polar plot of the TL
figure
plotshdpol( 'provinces.shd', 0, 0, 50 )
axis( [ -10 10 -10 10 ] )
caxisrev( [ 60 80 ] )

%%
bellhop3d provinces_side2D
plotshd( 'provinces_side2D.shd', 2, 1, 1 )

bellhop3d provinces_side3D
plotshd( 'provinces_side3D.shd', 2, 1, 2 )


%%
bellhop provinces_2D
plotshd( 'provinces_2D.shd', 2, 2, 1 )

scooter provinces_2DS
plotshd( 'provinces_2DS.shd.mat', 2, 2, 2 )

bellhop provinces_2Dshear
plotshd( 'provinces_2Dshear.shd', 2, 2, 3 )

scooter provinces_2DshearS
plotshd( 'provinces_2DshearS.shd.mat', 2, 2, 4 )
