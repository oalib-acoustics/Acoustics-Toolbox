function s = bpsk( s_bipolar, fc, fs, chips_per_sec )

% take the binary source sequence and encodes it as a BPSK signal

%fc = 12000;	carrier frequency
%fs = 48000;	sample frequency
%chips_per_sec = 3000;

% Michael B. Porter April 2000

samples_per_chip = fs / chips_per_sec;

if ( samples_per_chip ~= fix( samples_per_chip ) )
  error( 'samples_per_chip is not an integer' )
end

deltat = 1 / fs;
tsin   = 0 : deltat : ( samples_per_chip - 1 ) * deltat;

sinwave = sin( 2 * pi * fc * tsin )';	% sine wave as column vector
s = sinwave * s_bipolar;	% column * row is a matrix

s = reshape( s, 1, size( s, 1 ) * size( s, 2 ) );	% unravel it into a row vector


%Tfinal = length( s ) * deltat;
%t = 0 : deltat : Tfinal - deltat;
%plot( t(1:200), s(1:200 ) )
