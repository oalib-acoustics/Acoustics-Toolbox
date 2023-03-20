function p = HTS( kmin, deltak, stabil, r, G, option )

% Hankel transform of SCOOTER Green's function to produce pressure
% mbp, Dec. 2003

% G      = input k-space spectra, G( nrd, nt ) (destroyed by HTS)
% p      = output transmission loss, p( nrd, nr ).
% kmin   = minimum wavenumber
% deltak = spacing between wavenumbers
% stabil = amount integration has been moved off of the real axis.
% r      = vector of ranges for TL (must be evenly spaced)

if ndims( G ) > 2
   error( 'Hankel transform called with an n-dimensional array, n>3' )
end

if size( G, 2 ) == 1   % if G is a vector, expand it into a matrix with one row
    G = reshape( G, 1, length( G ) );
end

nk  = size( G, 2 );  % number of points in transform
nt2 = 2 * nk - 2;
k   = 0: deltak : ( nk - 1 ) * deltak;
r   = linspace( r( 1 ), r( end ), nt2 );
nrd = size( G, 1 );
ck  = stabil - 1i * kmin;
p   = zeros( nt2, nrd );

% scaling ...
if ( kmin > 0.5 * deltak )
   G( :, 1 ) =                        deltak   / sqrt( 2*pi ) * sqrt( kmin +             1i * stabil ) * G( :, 1 );
else
   G( :, 1 ) =  0.5 * ( 1i * stabil + deltak ) / sqrt( 2*pi ) * sqrt(                    1i * stabil ) * G( :, 1 );
end
G( :, 2:nk ) = scalecol( G( :, 2:nk ), deltak  / sqrt( 2*pi ) * sqrt( kmin + k( 2:nk ) + 1i * stabil ) );
G( :, 1    ) = real( G( :, 1 ) );

fftw( 'planner','patient' );   % set option to optimize fft on the first call

switch option( 2: 2 )
    case 'P'   % positive spectrum
        p = scalecol( G, exp( -1i * r( 1 ) * k + 1i * pi/4 ) ).'; % transposed so that FFT done down columns
        p = [ p; zeros( nk - 2, size( G, 1 ) ) ];
        p = fft( p, nt2 ).';    % exp (-IKX) transform; transposed so that each row is pressure vs. range
        p = scalecol( p, exp ( ck * r ) );
    case 'N'   % negative spectrum
        p = scalecol( G, exp( +1i * r( 1 ) * k - 1i * pi / 4 ) ).';
        p = [ zeros( size( p ) ); flipud( p( 2:nk-1, : ) ) ];
        p = fft( p, nt2 ).';  % exp (+IKX) transform with normalization factor, nt
        p = scalecol( p, exp ( -ck * r ) );
    case 'B'   % both positive and negative spectrum
        p = scalecol( G, exp( -1i * r( 1 ) * k + 1i * pi / 4 ) ).'; % transposed so that FFT done down columns
        p = [ p; zeros( nk - 2, size( G, 1 ) ) ];
        p = fft( p, nt2 ).';    % exp (-IKX) transform; transposed so that each row is pressure vs. range
        p = scalecol( p, exp ( ck * r ) );
        
        p2 = scalecol( G, exp( +1i * r( 1 ) * k - 1i * pi / 4 ) ).';
        p2 = [ zeros( size( p2 ) ); flipud( p2( 2:nk-1, : ) ) ];
        p2 = fft( p2, nt2 ).';  % exp (+IKX) transform with normalization factor, nt
        p = p + scalecol( p2, exp ( -ck * r ) );    
end

% cylindrical spreading
if ( option(1:1) == 'R' )
    %ii = find( r < eps( max( abs( r ) ) ) ); % look for zeros (or near-zeros) in the range vector
    %r( ii ) = eps( max( abs( r ) ) );        % get rid of them
    r( r < eps( max( abs( r ) ) ) ) = eps( max( abs( r ) ) ); % get rid of zeros in the range vector
    p = scalecol( p, 1 ./ sqrt( abs( r ) ) );
end
