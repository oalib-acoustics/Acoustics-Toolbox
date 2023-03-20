function p = HTS( kmin, deltak, stabil, r, G, option )

% Hankel transform of SCOOTER Green's function to produce pressure
% mbp, Dec. 2003

% G      = input k-space spectra, G( nrd, nt ) (destroyed by HTS)
% p      = output transmission loss, p( nrd, nr ).
% kmin   = minimum wavenumber
% deltak = spacing between wavenumbers
% stabil = amount integration has been moved off of the real axis.
% r      = vector of ranges for TL (must be evenly spaced)

% note: may be worth looking at Matlab tuning for fft (see fftw command)

if ndims( G ) > 2
   error( 'Hankel transform called with an n-dimensional array, n>3' )
end

if size( G, 2 ) == 1   % if G is a vector, expand it into a matrix with one row
    G = reshape( G, 1, length( G ) );
end

nt = size( G, 2 );  % number of points in transform
k  = 0: deltak : ( nt - 1 ) * deltak;
ck = stabil - 1i * kmin;
p  = zeros( size( G ) );

% scaling ...
if ( kmin > 0.5 * deltak )
    G( :, 1 ) =                      deltak   / sqrt( 2*pi ) * sqrt( kmin +             1i*stabil ) * G( :, 1 );
else
    G( :, 1 ) =   0.5 * ( 1i*stabil + deltak ) / sqrt( 2*pi ) * sqrt(                    1i*stabil ) * G( :, 1 );
end
G( :, 2:nt ) = scalecol( G( :, 2:nt ), deltak / sqrt( 2*pi ) * sqrt( kmin + k( 2:nt ) + 1i*stabil ) );

% *** First term, lower-half of cosine transform ***
if ( option(2:2) == 'P' || option(2:2) == 'B' )   % (P)ositive or (B)oth positive and negative spectrum
    Y = scalecol( G, exp( -1i * r( 1 ) * k + 1i*pi/4 ) ).';   % Y transposed so that FFT can be done down columns
    Y = fft( Y, nt ).';    % exp (-IKX) transform;        Y transposed so that each row is pressure vs. range
    p = scalecol( Y, exp ( ck * r ) );
end

% *** Second term, second-half of cosine transform ***
if ( option(2:2) == 'N' || option(2:2) == 'B' )   % (N)egative or (B)oth positive and negative spectrum
    Y = scalecol( G, exp( +1i * r( 1 ) * k - 1i*pi/4 ) ).';
    Y = nt * ifft( Y, nt ).';  % exp (+IKX) transform with normalization factor, nt
    p = p + scalecol( Y, exp ( -ck * r ) );
end

% cylindrical spreading
if ( option(1:1) == 'R' )
    ii = find( r < eps( max( abs( r ) ) ) ); % look for zeros (or near-zeros) in the range vector
    r( ii ) = eps( max( abs( r ) ) );        % get rid of them
    p = scalecol( p, 1 ./ sqrt( abs( r ) ) );
end
