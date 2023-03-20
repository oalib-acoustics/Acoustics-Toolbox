function [ s, time ] = hfm( fmin, fmax, T, sample_rate )

% generate an HFM pulse (hyperbolic frequency modulation, a.k.a. linear period modulation)
%
% fmin       = min frequency
% fmax       = max frequency
% T          = duration of time-series
% sample_rate = samples per second

N      = T * sample_rate;
deltat = T / N;

time = linspace( 0.0, T - deltat, N );

b = ( fmin - fmax )/( fmin * fmax * T );
P1 = 1 / fmin;
s = sin( ( 2 * pi / b ) * log( 1 + b * time / P1 ) );

%figure
%specgram( s, 256, sample_rate, [], [] );
