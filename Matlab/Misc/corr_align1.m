
function r_align = align( r )

% Line up the time series for an ensemble of pings
% Each row of r should contain a timeseries.
% mbp 7/96

Npings = size( r, 1 );
N      = size( r, 2 );

r_align( 1, : ) = r( 1, : );   % first ping is unshifted
ref = r( 1, : );               % this is the reference for alignment

shift_vec = zeros( Npings, 1 );

for ping = 2 : Npings

    % Calculate shift from correlation peak

    % ref = r_align( ping-1, : );

    temp = xcorr( r( ping, : ), ref );
    [ ~, shift_vec( ping ) ] = max( temp );   % locate peak of correlation
    shift_vec( ping ) = N - shift_vec( ping );

    shift = shift_vec( ping );

    if shift < 0
       shift = N + shift;
    end

    r_align( ping, : ) = [ r( ping, shift+1:N   ) r( ping, 1:shift   ) ];

    % ref = ref + r_align( ping, : );   % add aligned ping into reference ping

end   % next ping

