function plotts( filename )
% plot a time series from the given file

% open the file
fid = fopen( filename, 'r' );
if ( fid == -1 )
    error( 'No timeseries file with that name exists' );
end

% read the time series
[ PlotTitle, Pos, tout, RTS ] = read_ts( filename );
t   = tout;
nt  = length( t );
nrd = length( Pos.r.z );
rd  = Pos.r.z;

% plot
figure
orient tall
take = 1 : nt;
title( PlotTitle )
set( gca, 'YDir', 'Reverse' )

% scale all time series so that max of the absolute value is unity
scale  = max( max( abs( RTS( take, : ) ) ) );
RTS    = RTS / scale * rd( nrd ) / nrd;
offset = linspace( rd( 1 ), rd( nrd ), nrd );

hold on
threshold = 0;

for ird = 1 : nrd
   % ii = find( RTS( take, ird ) >  threshold );   % above threshold
   % jj = find( RTS( take, ird ) <= threshold );   % below threshold
   % area( t( take( ii ) ), RTS( take( ii ), ird ) + offset( ird ), offset( ird ) ); % plot part above threshold with shading under line
   % ylabel( [ 'Rd = ', num2str( rd( ird ) ) ] );
   % plot( t( take( ii ) ), rts( take( ii ), ird ) + offset( ird ) ); % plot part above threshold just as a line
   % area( t( take( jj ) ), RTS( take( jj ), ird ) + offset( ird ), offset( ird ) ); % plot part below threshold just as a line
   plot( t, RTS( :, ird ) + offset( ird ) ); % plot as a line
end

xlabel( 'Time (s)' );
%ylabel( [ 'Rd = ', num2str( rd( nrd ) ) ] );

%set(1,'PaperPosition', [ 0.25 0.00 5.5 7.0 ] )
%print -deps bellhop.ps

%%
% shaded color plot

%figure
%imagesc( rts )
%colorbar