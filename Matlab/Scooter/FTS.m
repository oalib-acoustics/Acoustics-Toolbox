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

nt  = size( G, 2 );  % number of points in transform
nt2 = 2 * nt - 2;

r = linspace( r( 1 ), r( end ), nt2 );
nr = size( G, 1 );

k  = 0 : deltak : ( nt - 1 ) * deltak;
ck = -kmin - 1i * stabil;

p  = zeros( nt2, nr );

% scaling ...
%if ( kmin > 0.5 * deltak )
%   G( :, 1 ) =                         deltak   / sqrt( 2*pi ) * G( :, 1 );
%else
%   G( :, 1 ) =      0.5 * ( i*stabil + deltak ) / sqrt( 2*pi ) * G( :, 1 );
%end
G( :, 1:nt ) = scalecol( G( :, 1 : nt ), deltak   / sqrt( 2*pi ) );
G( :, 1 )    = real( G( :, 1 ) );

switch option( 2: 2 )
    case 'P'
        Y = scalecol( G, exp( -1i * r( 1 ) * k ) ).';     % Y transposed so that FFT can be done down columns
        Y = [ Y; zeros( nt - 2, size( G, 1 ) ) ];
        Y = fft( Y, nt2 ).';    % exp (-IKX) transform; Y transposed so that each row is pressure vs. range
        p = scalecol( Y, exp ( 1i * ck * r ) );
    case 'N'
        Y = scalecol( G, exp( -1i * r( 1 ) * k ) ).';
        Y = [ zeros( size( Y ) ); flipud( conj( Y( 2:nt-1, : ) ) ) ];
        Y = fft( Y, nt2 ).';  % exp (+IKX) transform with normalization factor, nt
        p = scalecol( Y, exp ( -1i * ck * r ) );
    case 'B'
        Y = scalecol( G, exp( -1i * r( 1 ) * k ) ).';     % Y transposed so that FFT can be done down columns
        Y = [ Y; zeros( nt - 2, size( G, 1 ) ) ];
        Y = fft( Y, nt2 ).';    % exp (-IKX) transform; Y transposed so that each row is pressure vs. range
        p = scalecol( Y, exp ( 1i * ck * r ) );
 
        Y = scalecol( G, exp( -1i * r( 1 ) * k ) ).';
        Y = [ zeros( size( Y ) ); flipud( Y( 2:nt-1, : ) ) ];
        Y = fft( Y, nt2 ).';  % exp (+IKX) transform with normalization factor, nt
        p = p + scalecol( Y, exp ( 1i * ck * r ) );
        
        %Y = [ Y; flipud( conj( Y( 2 : nt - 1, : ) ) ) ];
        %G( :, nt+2, 2*nt ) = fliplr( G( :, 1:nt ) ).';
        %Y = fft( Y, nt2 ).';    % exp (-IKX) transform; Y transposed so that each row is pressure vs. range
        %p = scalecol( Y, exp ( i * ck * r ) );
        %p = Y;
end
