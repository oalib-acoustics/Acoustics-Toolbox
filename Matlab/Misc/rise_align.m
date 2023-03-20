
function r_align = rise_align( r, Threshold, Offset )

% usage: r_align = rise_align( r, Threshold, Offset )
%
% Line up the time series for an ensemble of pings
% Each row of r should contain a timeseries.
% This version uses the leading edge to align things.
%
% r is the ensemble of time series.
% Threshold is used to define the leading edge.
%    (= 0.9 means that the leading edge is a 0.9 times the peak)
% Offset is the sample number where the peak is placed
%
% mbp 7/96

Npings = size( r, 1 );
N      = size( r, 2 );
Offset = fix( Offset );	% offset for leading edge on plot (in samples)

for ping = 1:Npings

    % Calculate shift from correlation peak

    temp = r( ping, : );
    [ Peak, ipeak ] = max(  abs( temp ) );
    Avg  = norm( temp ) / length( temp );
    J = find( abs( temp( 1: ipeak ) ) > Avg + Threshold * ( Peak - Avg ) );

    if isempty( J )   % this happens if time-series is all zeros
       shift = 0;
     else
       % note that a shift of -m is equivalent to a shift of N-m
       shift = mod( J( 1 ) - Offset, N );   % put the peak at Offset
    end

    r_align( ping, : ) = [ r( ping, shift+1:N   ) r( ping, 1:shift   ) ];

end   % next ping

