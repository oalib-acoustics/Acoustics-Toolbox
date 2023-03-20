function plotts( filename )
% plot a time series from the given file

% open the file
fid = fopen( filename, 'r' );

% read

pltitl = fgetl( fid );
nrd    = fscanf( fid, '%f', 1 );

rd     = fscanf( fid, '%f', nrd );
temp   = fscanf( fid, '%f', [ nrd + 1, inf ] );
fclose( fid );

% extract rts

t = temp(1, :)';
nt = length( t );

rts = temp(2:nrd+1, :)';

% plot

figure
orient tall

take = 1:nt;
title( pltitl )

for ird = 1:nrd
   subplot( nrd, 1, ird );
   plot( t( take ), rts( take, ird ) );
   ylabel( [ 'Rd = ', num2str( rd( ird ) ) ] );
end

xlabel( 'Time (s)' );
ylabel( [ 'Rd = ', num2str( rd( nrd ) ) ] );

set(1,'PaperPosition', [ 0.25 0.00 5.5 7.0 ] )
%print -deps bellhop.ps
