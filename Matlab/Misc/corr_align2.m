
function r_align = align( r )

% Line up the time series for an ensemble of pings
% Each row of r should contain a timeseries.

% This version aligns each pulse relative to a previous one
% trying to maximize the correlation between the two pulses
% mbp 7/96

Npings = size( r, 1 );
N      = size( r, 2 );

r_align( 1, : ) = r( 1, : );

for ping = 2:Npings
   
   % Calculate shift from correlation peak
   
   % align against previous ping
   if ping == 2 || ( max( r_align( ping-1, : ) ) >0 && ping == 10 * floor( ping / 10 ) )
      oldfft = fft( r_align( ping-1, : ) );
   end
   
   % for ping 134 ship has moved back towards the testbed
   if ping == 134
      oldfft = fft( r_align( 32, : ) );
   end
   
   newfft = fft( r( ping, : ) );
   corr   = ifft( conj( oldfft ) .* newfft );
   [ ~, PeakIndex ] = max( corr );
   shift  = PeakIndex - 1;
   
   r_align( ping, : ) = [ r( ping, shift+1:N   ) r( ping, 1:shift   ) ];
   
end   % next ping

