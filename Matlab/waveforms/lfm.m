function [ s, time ] = lfm( fmin, fmax, T, sample_rate )

% generate an LFM pulse (chirp)
%
% fmin        = min frequency
% fmax        = max frequency
% T           = duration of time-series
% sample_rate = samples per second

N          = T * sample_rate;
deltat     = T / N;
time       = linspace( 0.0, T - deltat, N );
feffective = fmin + ( fmax - fmin ) .* time / ( 2 * T );

s = sin( 2.0 * pi * ( feffective .* time  ) );

% figure; specgram( s, 256, sample_rate, [], [] );

% above is more simply done with the voltage controlled oscillator:
%
% x = linspace( -1, 1, T * samplerate );
% s = vco( x, [fmin, fmax], samplerate );
%
% or matlab chirp command

