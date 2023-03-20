function [ s ] = mseq( m )

% generate an m-sequence
%
% Michael B. Porter April 2000
% formulas from Proakis, Digital Communications

if ( m < 2 || m > 15 || m ~= fix( m ) ); return; end

% set up the coefficients ( used in the recursion: m_l = \sum c_l m_l )

switch m
  
case  2, c = [ 1 1 ];
case  3, c = [ 1 0 1 ];
case  4, c = [ 1 0 0 1 ];
case  5, c = [ 1 0 0 1 0 ];
case  6, c = [ 1 0 0 0 0 1 ];
case  7, c = [ 1 0 0 0 0 0 1 ];
case  8, c = [ 1 0 0 0 1 1 1 0 ];
case  9, c = [ 1 0 0 0 0 1 0 0 0 ];
case 10, c = [ 1 0 0 0 0 0 0 1 0 0 ];
case 11, c = [ 1 0 0 0 0 0 0 0 0 1 0];
case 12, c = [ 1 0 0 0 0 0 1 0 1 0 0 1];
case 13, c = [ 1 0 0 0 0 0 0 0 0 1 1 0 1];
case 14, c = [ 1 0 0 0 1 0 0 0 1 0 0 0 0 1 ];
case 15, c = [ 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ];
end

length = 2^m - 1;

%successive shifts with feedback following Proakis p. 433

seed = [ 1 zeros( 1, m - 1 ) ];	% all zero seed except for one
s    = zeros( 1, length );

for ii = 1 : length
  out( 1 : m - 1 ) = seed( 2 : m );
  out( m )         = mod( c * seed' , 2 );	% addition mod 2
  seed             = out;
  s( ii )          = out( 1 );
end

% convert to a +/- sequence

%ii = find( s == 0 );
%s( ii ) = -1;
s( s == 0 ) = -1;

% check autocorrelation properties

%shat = fft( s );
%scorr = real( ifft( shat .* conj( shat ) ) );
%plot( scorr )
