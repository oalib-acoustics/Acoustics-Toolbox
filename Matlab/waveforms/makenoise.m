function nts = makenoise( fc, BW, T, fs )

% make a waveform consisting of bandpass-filtered, Gaussian random noise
%
% useage:
% sts = make_noise( f, T, fs )
%
% input:
% f  = center frequency of band
% BW = Bandwidth
% A  = associated amplitude
% T  = duration
% fs = sample rate
%
% output:
%    nts is a column vector with the noise timeseries
%
% mbp 27 September 2007 for PN07 experiment

deltat = 1 / fs;
N      = T / deltat;             % number of samples to generate
time   = 0:deltat:T - deltat;
time   = time';               % make time a column vector

deltat2 = 1 / BW;
N2      = T / deltat2;           % number of samples to generate

nts = randn( N2, 1 );                       % generate basebanded waveform
nts = resample( nts, N, N2 );               % resample to fs rate
nts = sin( 2 * pi * fc * time ) .* nts;     % hetorodyne with carrier frequency
