function y = gaussian( time, delay, duration )

% generate a Gaussian pulse
% usage: y = gaussian( time, delay, duration )
%
% time is a vector of sample times
% delay is time of the pulse (peak location)
% duration is, well, the duration
% time, delay, duration should all be in the same units, e.g. seconds
%
% mbp 2001

y = exp( -( ( time - delay ) / duration ) .^2 );
