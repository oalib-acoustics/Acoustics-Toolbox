function rootz = PekerisRoot( z )

% Calculates the Pekeris branch of the square root
% mbp 4/2009

ii = find( real( z ) >= 0.0 );
rootz ( ii ) = sqrt( z( ii ) );

ii = find( real( z ) < 0.0 );
rootz( ii ) = 1i * sqrt( -z( ii ) );

rootz = rootz( : );   % make output a column vector

