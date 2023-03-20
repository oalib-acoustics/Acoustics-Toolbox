function noise = noise_maker( power, len )
% usage: noise = noise_maker( power,len )
% produces a Gaussian random noise vector with specified power and length
% power is given in dB
% the dB reference is a square wave of unit amplitude
%

noise = 10^( power / 20 ) * randn( len, 1 );     % add Gaussian noise
 