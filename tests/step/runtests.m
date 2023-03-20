% Based on the Dickins seamount, but with a simple step up in the
% bathymetry
% mbp

global units jkpsflag

units = 'km';

Nprof = 4;  % change field.flp to also
% r     = linspace( 0.0, 30000, Nprof );
r = [ 0 50000 51000 100000 ];

kraken_rd( 'stepK', r )

%%
% Fortran version of field
%delete field.prt
%delete fort.6
runfield  = which( 'field.exe' );
eval( [ '! "' runfield '" stepK_rd' ] );

%%
figure
plotshd( 'stepK_rd.shd' )
caxisrev( [ 60 110 ] )

%%
% Matlab version of field
% field( 'DickinsK_rd' )
% 
% figure
% plotshd( 'DickinsK_rd.mat' )
% caxisrev( [ 70 120 ] )

%%
% simplePE run

simplePE step
figure
plotshd step.shd.mat
caxisrev( [ 60 110 ] )

