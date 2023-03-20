% Testing the power law attenuation in sediments

% The calib case compares TL in shallow water with 3 sediment attenuation
% laws:
%    1) no sediment loss
%    2) 0.6 dB/wavelength
%    3) quadratic power law with 0.1 dB at the base frequency of 250 Hz
% TL increases progressively with those changes except at the base
% frequency of 250 Hz. The power law gives the same result as the 0.6
% dB/wavelength case because it has the same attenuation at 250 Hz
%
% mbp

global units
units = 'km';

%%
% Calibration case

% reflection loss (compare to plot in JKPS)

bounce calibBounce
plotbrc 'calibBounce.brc'

%%

scooter calibS_noloss

plotshd( 'calibS_noloss.shd.mat', 250, 3, 1, 1 )
caxis( [ 40 80 ] )

plotshd( 'calibS_noloss.shd.mat', 500, 3, 1, 2 )
caxis( [ 40 80 ] )

plotshd( 'calibS_noloss.shd.mat', 1000, 3, 1, 3 )
caxis( [ 40 80 ] )

%%
scooter calibS_0.6dB

plotshd( 'calibS_0.6dB.shd.mat', 250, 3, 1, 1 )
caxis( [ 40 80 ] )

plotshd( 'calibS_0.6dB.shd.mat', 500, 3, 1, 2 )
caxis( [ 40 80 ] )

plotshd( 'calibS_0.6dB.shd.mat', 1000, 3, 1, 3 )
caxis( [ 40 80 ] )

%%
scooter calibS_PowLaw

plotshd( 'calibS_PowLaw.shd.mat', 250, 3, 1, 1 )
caxis( [ 40 80 ] )

plotshd( 'calibS_PowLaw.shd.mat', 500, 3, 1, 2 )
caxis( [ 40 80 ] )

plotshd( 'calibS_PowLaw.shd.mat', 1000, 3, 1, 3 )
caxis( [ 40 80 ] )
