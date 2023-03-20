function plotshd_multi( filename, m, n, p )

% Like plotshd.m, but works on 3D shade files. For example, if you
% generate TL images for multiple sources, or TL images for
% multiple time slices (using sparc), the pressure read from the
% shade file will be 3D. A 3D shade file will break plotshd.
%
% Say you have a 3D SHDFIL with 6 TL images. Here is how you would
% use plotshd_multi to display all 6 TL images in a single Matlab
% figure:
% for ii=1:6,
%    plotshd_new('SHDFIL',2,3,ii);
% end;
%
% Parameters m, n, p are OPTIONAL - if they are not specified, 
% plotshd_multi behaves like plotshd.
%
% Adapted from plotshd.m, Paul Hursky, July 24, 2002.
% usage:
% plotshd_new( filename, m, n, p )

% read

[ pltitl, ~, ~, ~, ~, nrr, ~, rd, rr, pressure ] = read_shd( filename );

zt = rd;
taker = 1:nrr;
rt = rr( taker );

rkm = rt / 1000.0;

if ( nargin == 1 )
  figure
else
  if ( p == 1 )
    figure; 
  end
  subplot( m, n, p )
end

tlt = abs( pressure(:,:,p) );
ii = find( tlt == 0 );                      % find places where pressure vanishes
icount = find( tlt > 1e-6 );                % for stats, only these values count
tlt( ii ) = max( max( tlt ) )/ 1e10; % and set it to a small number
tlt = -20.0 * log10( tlt );          % so there's no error when we take the log

% compute some statistics to automatically set the color bar

tlmed = median( tlt( icount ) );    % median value
tlstd = std( tlt( icount ) );       % standard deviation
tlmax = tlmed + 0.75 * tlstd;       % max for colorbar
tlmax = 10 * round( tlmax / 10 );   % make sure the limits are round numbers
tlmin = tlmax - 50;                 % min for colorbar

% optionally remove cylindrical spreading:
% tlt = tlt + ones( nrd, 1 ) * 10.0 * log10( rt )';

tej = flipud( jet );  % 'jet' colormap reversed
if ( size( tlt, 1 ) > 1 && size( tlt, 2 ) > 1 )
    pcolor( rt, zt, tlt );  ...
        shading flat; colormap( tej );
    caxis( [ tlmin, tlmax ] ); colorbar( 'horiz' ); view( 0, -90 );
    xlabel( 'Range (m)' ); ylabel( 'Depth (m)' );
    title( deblank( pltitl ) )
    hold on
else
    if ( size( tlt, 1 ) == 1 )
        plot( rt, tlt )
        xlabel( 'Range (m)' ); ylabel( 'TL (dB)' ); view( 0, -90 )
    else
        plot( tlt, zt )
        view( 0, -90 );
        xlabel( 'TL (dB)' ); ylabel( 'Depth (m)' );
    end
end

set( gcf, 'PaperPosition', [ 0.25 0.25 6.0 3.0 ] )
