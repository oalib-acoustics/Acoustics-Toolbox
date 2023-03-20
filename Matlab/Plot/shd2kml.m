function varargout = plotshd( filename )

% convert a shd file to a single TL surface in dB, in kml format (for
% Google Earth)
%
% usage:
% shd2kml( filename )
% (m, n, p) optional subplot spec
% '.shd' is the fault file extension if not specified
%
% mbp

'just started on this; however, we don't have a package for displaying a range-depth plot in Google Earth'

global units jkpsflag

% read

%disp( 'PlotShd uses the first bearing and source depth in the shade file; check OK' )
itheta = 1;
isd    = 1;

[ PlotTitle, ~, freq, ~, Pos, pressure ] = read_shd( filename );

pressure = squeeze( pressure( itheta, isd, :, : ) );
zt       = Pos.r.z;
rt       = Pos.r.r;

% set labels in m or km
xlab     = 'Range (m)';
if ( strcmp( units, 'km' ) )
    rt      = rt / 1000.0;
    xlab    = 'Range (km)';
end


%%

tlt = double( abs( pressure ) );   % need 'double' because field.m produces a single precision result and pcolor needs double

tlt( isnan( tlt ) ) = 1e-6;   % remove NaNs
tlt( isinf( tlt ) ) = 1e-6;   % remove infinities

icount = find( tlt > 1e-37 );         % for stats, only these values count
tlt( tlt < 1e-37 ) = 1e-37;            % remove zeros
tlt = -20.0 * log10( tlt );          % so there's no error when we take the log

% compute some statistics to automatically set the color bar

tlmed = median( tlt( icount ) );    % median value
tlstd = std( tlt( icount ) );       % standard deviation
tlmax = tlmed + 0.75 * tlstd;       % max for colorbar
tlmax = 10 * round( tlmax / 10 );   % make sure the limits are round numbers
tlmin = tlmax - 50;                 % min for colorbar

% optionally remove cylindrical spreading:
% tlt = tlt + ones( nrd, 1 ) * 10.0 * log10( rt )';
%%

tej = flipud( jet( 256 ) );  % 'jet' colormap reversed

if ( size( tlt, 1 ) > 1 && size( tlt, 2 ) > 1 )
    % imagesc produces a better PostScript file, using PostScript fonts
    % however, it ignores the actual r, z, coordinates and assumes they're
    % equispaced
    %h = imagesc( rt, zt, tlt );
    h = pcolor( rt, zt, tlt );  ...
    shading flat
    colormap( tej )
    caxisrev( [ tlmin, tlmax ] )
    set( gca, 'YDir', 'Reverse' )
    xlabel( xlab )
    ylabel( 'Depth (m)' );
    title( { deblank( PlotTitle ); [ 'Freq = ' num2str( freq ) ' Hz    Sd = ' num2str( Pos.s.z ) ' m' ] } )
else   % line plots
    if ( size( Pos.r.r, 1 ) > 1 )   % TL vs. range
        h = plot( rt, tlt );
        xlabel( xlab );
        ylabel( 'TL (dB)' )
        set( gca, 'YDir', 'Reverse' )
        title( deblank( PlotTitle ) )
    else
        % TL vs. depth
        h = plot( tlt', zt );
        set( gca, 'YDir', 'Reverse' )
        set( gca, 'Xdir', 'Reverse' )
        xlabel( 'TL (dB)' )
        ylabel( 'Depth (m)' );
        title( deblank( PlotTitle ) )
    end
end

%text( 0.98 * max( rt ), min( zt ), '(a)' );

drawnow

if ( nargout == 1 )
    varargout( 1 ) = { h };   % return a handle to the figure
end

% fixed size for publications
if ( jkpsflag )
    set( gca, 'ActivePositionProperty', 'Position', 'Units', 'centimeters' )
    set( gcf, 'Units', 'centimeters' )
    set( gcf, 'PaperPositionMode', 'auto');   % this is important; default is 6x8 inch page

    if ( exist( 'm' ) )
        set( gca, 'Position', [ 2    2 + ( m - p ) * 9.0     14.0       7.0 ] )
        set( gcf, 'Position', [ 3                   15.0     19.0  m * 10.0 ] )
    else
        set( gca, 'Position', [ 2    2                       14.0       7.0 ] )
        set( gcf, 'Units', 'centimeters' )
        set( gcf, 'Position', [ 3 15 19.0 10.0 ] )
    end
    
    %     set( gcf, 'Units', 'centimeters' )
    %     set( gcf, 'PaperPositionMode', 'manual' );
    %     set( gcf, 'PaperPosition', [ 3 3 15.0 10.0 ] )
    
end
