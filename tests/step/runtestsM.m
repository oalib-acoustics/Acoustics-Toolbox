% Based on the Dickins seamount, but with a simple step up in the
% bathymetry
% mbp

% Oct. 2020: I'm getting sporadic segmentation faults; however, the code
% seems to successfully run to completion anyway
% It happens with gfortran but not the Intel compiler
% Probably there is a bug somewhere ...

global units jkpsflag

units = 'km';

Nprof = 4;  % change field.flp to also
% r     = linspace( 0.0, 30000, Nprof );
r = [ 0 50000 51000 100000 ];

kraken_rd( 'stepK', r )


%%
% Matlab version of field
field( 'stepK_rd' )

figure
plotshd( 'stepK_rd.shd.mat' )
caxisrev( [ 60 110 ] )

%%
% simplePE run

simplePE step
figure
plotshd step.shd.mat
caxisrev( [ 60 110 ] )

