% Dickins seamount test
% mbp
global units jkpsflag

units = 'km';

Nprof = 601;  % change field.flp also
r     = linspace( 0.0, 30000, Nprof );

kraken_rd( 'DickinsK', r )

%%
% Fortran version of field

tic
runfield  = which( 'field.exe' );
eval( [ '! "' runfield '" DickinsK_rd' ] );
toc

figure
plotshd( 'DickinsK_rd.shd' )
caxisrev( [ 70 120 ] )

%%
% Matlab version of field

tic
field( 'DickinsK_rd' )
toc

figure
plotshd( 'DickinsK_rd.shd.mat' )
caxisrev( [ 70 120 ] )
