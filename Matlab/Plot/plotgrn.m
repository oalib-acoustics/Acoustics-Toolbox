function plotgrn( filename )
% plot the Green's function calculated by SCOOTER
%
% usage:
% plotgrn( filename, m, n, p )
% (m, n, p) optional subplot spec
%
% mbp, Dec. 2003

isd = 1;

% read

[ PlotTitle, ~, freqVec, ~, ~, Pos, G ] = read_shd( filename );
k  = 2 * pi * freqVec( 1 ) ./ Pos.r.r;    % k values are stored in the range-vector

Gmag = abs( squeeze( G( 1, isd, :, : ) ) );
%Gmag = abs( squeeze( G( 1, isd, 251, : ) ) ).';
%Gmag = abs( squeeze( G( 1, isd, :, 4000 ) ) ).';

if ( size( Gmag, 1 ) > 1 && size( Gmag , 2 ) > 1 )
    pcolor( k, Pos.r.z, Gmag );  ...
        shading flat; colormap( jet );
    caxis( [ 0, max( max( Gmag ) ) ] ); colorbar( 'horiz' );
    set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
    xlabel( 'Wavenumber (1/m)' )
    ylabel( 'Depth (m)' );
    title( deblank( PlotTitle ) )
else
    if ( size( Gmag, 1 ) == 1 )
        plot( k, Gmag )
        %plot( Pos.r.z, Gmag )
        xlabel( 'Wavenumber (1/m)' )
        ylabel( '|G|' )
    end
end
