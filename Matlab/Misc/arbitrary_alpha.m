function s_alpha = arbitrary_alpha( s, alpha )

% function s_alpha = arbitrary_alpha( s, alpha )
%
% time scales the input signal with an arbitrary scaling
% factor, alpha:
%
%      s_alpha(t) = s( alpha*t )
%
% the algorithm should work well for values of alpha near 1
% a more general approach maybe required for cases when alpha
% is not near 1
%
% the modified chirp z-transform is used to compute
% the required frequency domain samples of the 
% time scaled waveform
%

% J.M.Tattersall NUWC 9/5/96

% note alpha > 1 means time series compressed, i.e. closing in range/ mbp
% s should be a row vector

nt = length( s );

if alpha == 1
  s_alpha = s;
  return
end

Nfft = next_radix2( length(s) );
s    = [ s zeros(1, Nfft - length( s ) ) ];   % zero pad to power of 2
ups  = 1;


if alpha < 1
  s = interp( s, 2);
  ups = 2;
end

delta_phi = 2*pi/( ups * Nfft * alpha );

phi_0 = 0;

S_alpha = modified_chirp_z( s, Nfft, delta_phi, phi_0 );
S_alpha = S_alpha .* [ ones( 1, Nfft / 2 + 1 ) zeros( 1, Nfft / 2 - 1 ) ];
s_alpha = real( ( 2 / ups ) * ( 1 / alpha ) * ifft( S_alpha, Nfft ) );

s_alpha = s_alpha( 1, 1 : nt );   % trim off zero-filled part
