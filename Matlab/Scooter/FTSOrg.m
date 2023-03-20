function p = FTS( kmin, deltak, stabil, r, G, option )

% Fourier transform of SCOOTER Green's function to produce pressure
% mbp, Dec. 2003

% G      = input k-space spectra, G( nrd, nt ) (destroyed by HTS)
% p      = output transmission loss, p( nrd, nr ).
% kmin   = minimum wavenumber
% deltak = spacing between wavenumbers
% stabil = amount integration has been moved off of the real axis.
% r      = vector of ranges for TL (must be evenly spaced)

if ndims( G ) > 2
   error( 'Fourier transform called with an n-dimensional array, n>3' )
end

if size( G, 2 ) == 1   % if G is a vector, expand it into a matrix with one row
    G = reshape( G, 1, length( G ) );
end

nt = size( G, 2 );  % number of points in transform
k  = 0 : deltak : ( nt - 1 ) * deltak;
ck = -kmin - 1i * stabil;
p = zeros( size( G ) );

% scaling ...
%if ( kmin > 0.5 * deltak )
%   G( :, 1 ) =                         deltak   / sqrt( 2*pi ) * G( :, 1 );
%else
%   G( :, 1 ) =      0.5 * ( i*stabil + deltak ) / sqrt( 2*pi ) * G( :, 1 );
%end
G( :, 1:nt ) = scalecol( G( :, 1:nt ), deltak   / sqrt( 2*pi ) );

% *** First term, lower-half of cosine transform ***
if ( option(2:2) == 'P' || option(2:2) == 'B' )   % (P)ositive or (B)oth positive and negative spectrum
    Y = scalecol( G, exp( -1i * r( 1 ) * k ) ).';     % Y transposed so that FFT can be done down columns
    Y = fft( Y, nt ).';    % exp (-IKX) transform; Y transposed so that each row is pressure vs. range
    p = scalecol( Y, exp ( 1i * ck * r ) );
end

% *** Second term, second-half of cosine transform ***
if ( option(2:2) == 'N' || option(2:2) == 'B' )   % (N)egative or (B)oth positive and negative spectrum
    Y = scalecol( G, exp( +1i * r( 1 ) * k ) ).';
    Y = nt * ifft( Y, nt ).';  % exp (+IKX) transform with normalization factor, nt
    p = p + scalecol( Y, exp ( -1i * ck * r ) );
end