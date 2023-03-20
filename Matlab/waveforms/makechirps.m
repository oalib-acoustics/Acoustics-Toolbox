function probe = makechirps( fmin, fmax, fs, T, shortcleartime, T_tot )
% makechirps( fmin, fmax, fs, T, shortcleartime, T_tot )
% [ fmin, fmax ] frequency band
% fs = sampling rate
% T = chirp duration
% T_tot = total duration for chirp sequence

% 2'' leader
% 232 lfm's
% comb
% trailer

% mbp 11/29/99

lead_time = 0.2;
chirp_duration = T + shortcleartime;
Nprobes    = floor( ( T_tot - lead_time ) / chirp_duration );

% lfm's

[ s, time ] = lfm( fmin, fmax, T, fs );

% assemble chirps
shortclear = zeros( 1, shortcleartime * fs );

probe = [];
for iprobe = 1:Nprobes
  probe = [ probe, s, shortclear ];
end

% normalize and zero fill

leader  = zeros( 1, lead_time * fs );
probe = [ leader 0.95 * probe / max( abs( probe ) ) ];  % normalize
n = length( probe );
if ( n < T_tot * fs )
  probe = [ probe, zeros( 1, T_tot * fs - n ) ];   % zero-fill to T_tot
end
