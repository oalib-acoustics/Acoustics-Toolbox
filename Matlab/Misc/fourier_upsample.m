function upsampled_ts = fourier_upsample( time_series, ratio )

% up sample a time series by zero padding in freq domain
%
% Usage: upsampled_ts = fourier_upsample(time_series, ratio)
%
%  Input arguments:
%   time_series  - real valued vector containing the time series
%   ratio        - integer valued scalar containing the upsample ratio 
%
%  Output argument:
%   upsampled_ts - real valued vector containing the upsampled time series
%
% John Peterson 2012?

% Check for missing input arguments

if nargin < 2
  error( [mfilename, ': one or more input arguments is missing.'] );
end

% Some sanity checks on the time series

if ndims(time_series) > 2
  error( [mfilename, ': time series must be a row or col vector.'] );
end

[nrow, ncol] = size(time_series);

if min([nrow, ncol]) > 1
  error( [mfilename, ': time series must be a row or col vector.'] );
end

len_ts = length(time_series);

% Some sanity checks on the upsample ratio

if ratio < 1
  error( [mfilename, ': upsample ratio must be a positive integer.'] );
end

if (ratio - floor(ratio)) > 0
  error( [mfilename, ': upsample ratio is not integer valued.'] );
end

% Compute the discrete fourier transform of the time series

fft_ts = fft(time_series);

% Allocate memory for the fourier transform of the upsampled time series

len_usts = ratio * len_ts;

fft_usts = complex(zeros(len_usts, 1), zeros(len_usts, 1));

% Copy the frequency bins from the transform of the original time series

fft_usts(1                       : 1 + len_ts/2) = fft_ts( 1            : 1 + len_ts/2);
fft_usts(len_usts - len_ts/2 + 2 : len_usts    ) = fft_ts( 2 + len_ts/2 : len_ts      );

% The handling of the Nyquist frequency is a little subtle. In the transform
% of the original time series, the amp of the Nyquist frequency is contained
% in a single bin. In the transform of the upsampled time series, the amp of
% the Nyquist frequency is contained in TWO bins. The orginal amp is divided
% by two and placed in the correct bins of the padded transform.

amp_nyquist = 0.5 * fft_usts(1 + len_ts/2);

fft_usts(1 + len_ts/2           ) = amp_nyquist;
fft_usts(len_usts - len_ts/2 + 1) = amp_nyquist;

% The MATLAB [sic] discrete fourier transform includes a weighting by the
% length of the transform, which must be adjusted to reflect the new length.

fft_usts = ratio * fft_usts;

% Compute the inverse fourier transform to get the upsampled time domain series

upsampled_ts = ifft(fft_usts);
