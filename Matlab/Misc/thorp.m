function alpha = thorp( f )

% Thorp volume attenuation (dB/km)
% Updated formula from JKPS Eq. 1.34
%
% f is a frequency vector in kHz
% alpha = attenuation   (dB/km)

alpha = 0.0033 + 0.11 * f.^2 ./ ( 1 + f.^2 ) + 44 * f.^2 ./ ( 4100 + f.^2 ) + 0.0003 * f.^2;
