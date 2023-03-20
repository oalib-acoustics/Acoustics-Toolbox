function y = whiten( x, ~, ~ )

% usage: y = whiten( x )
%
% x is an input  time series
% y is an output time series
%
% mbp 10/00, Barletta, Italy
%
N = 2^13;				% basic transform size
Nframes_avg = 20;		% # of frames to average for spectrum

% break the vector into frames

group_len = Nframes_avg * N;	% group length in samples

Npts		= length( x );
Ngroups 	= floor( Npts / group_len ) ;		% # of groups (a group is a set of frames, pre-whitened using the same estimated spectrum)
Nframes	= Ngroups * Nframes_avg;
Nused		= Npts * group_len;

xfft = fft( reshape( x( 1 : N * Nframes ), [ N Nframes ] ) );
clear x

% FFT and form the average spectrum

xmag = abs( xfft );
xmag = reshape( xmag, [ N Nframes_avg Ngroups ] );

mean_spec = squeeze( mean( xmag, 2 ) );	% mean across Nframes_avg (2nd dimension)
clear xmag

% whiten each group of frames, using the mean spectrum for that group

for igroup = 1 : Ngroups
  iframe1 = 1 + ( igroup - 1 ) * Nframes_avg
  iframe2 = iframe1 + Nframes_avg - 1;
  for iframe = iframe1 : iframe2
    xfft( :, iframe ) = xfft( :, iframe ) ./ mean_spec( :, igroup );
  end
end

y = real( ifft( xfft ) );					% go back to the time domain
clear xfft

y = reshape( y, [ N * Nframes 1 ] );	% reshape frames into a long vector

 