function plotts( filename )
% plot a time series from the given .mat file

% open the file
load( [ filename '.rts.mat' ] )
rz  = Pos.r.z;
nrz = length( rz );
t   = tout;
RTS = RTS';

% plot
figure
orient tall
title( PlotTitle )
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature

% scale all time series so that max is unity
scale  = max( max( RTS( :, : ) ) );
RTS    = RTS / scale * rz( nrz )/nrz;
offset = linspace( rz( 1 ), rz( nrz ), nrz );

hold on
threshold = -1e30;

for irz = 1 : nrz
   ii = find( RTS( :, irz ) >  threshold );
   jj = find( RTS( :, irz ) <= threshold );
   h  = area( t( ii ), RTS( ii, irz ) + offset( irz ) ); % pos. part shading under line
   %ylabel( [ 'Rd = ', num2str( rd( ird ) ) ] );
   set( h, 'BaseValue', offset( irz ) );
   plot( t( jj ), RTS( jj, irz ) + offset( irz ) ); % negative part just line
end

xlabel( 'Time (s)' );
%ylabel( [ 'Rz = ', num2str( rz( nrz ) ) ] );

%set(1,'PaperPosition', [ 0.25 0.00 5.5 7.0 ] )
