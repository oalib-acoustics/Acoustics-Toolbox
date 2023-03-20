function s = Nwave( time, F )

% make an N-wave
% peak at F, support [0, 4F], [0,3F] also OK
% F = nominal source frequency

omega = 2 * pi * F;

%IF ( T .LE. 1.0 / F )

s = sin( omega * time ) - 0.5 * sin( 2 * omega * time );

%ii = find( time > 1/F | time < 0 );
%s( ii ) = 0;

s( time > 1/F | time < 0 ) = 0;
