function probe = makemseq( fmin, fmax, fs, T_tot )
% make an m-sequence probe
%  makemseq( fmin, fmax, fs, T_tot )

% leader
% m-sequences
% trailer

% mbp 11/29/99

lead_time       = 0.2;

% m-sequence
fc = 0.5 * ( fmin + fmax );
chips_per_sec = 0.5 * ( fmax - fmin );

s_m = mseq( 10 );
s = bpsk( s_m, fc, fs, chips_per_sec );

Nreps = floor( T_tot * chips_per_sec / length( s_m ) );
probe = [];
for irep = 1 : Nreps
  probe = [ probe, s ];
end

leader  = zeros( 1, lead_time * fs );
probe = [ leader 0.95 * probe / max( abs( probe ) ) ];  % normalize

n = length( probe );
if ( n < T_tot * fs )
  probe = [ probe, zeros( 1, T_tot * fs - n ) ];   % zero-fill to 1'
end