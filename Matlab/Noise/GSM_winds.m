
function NL = GSM_winds( freq, W )

% GSM wind noise from Eller and Cavanaugh report
% result is a PSD (1 Hz band) dB re 1 microPascal^2 at 1 m
%
% freq is the frequency in Hz
% W is the wind speed in knots

ii = find( freq <  1000 );
jj = find( freq >= 1000 );

NL( ii ) = 44 + sqrt( 21 * W ) - 17 * ( log10( freq( ii ) ) - 3 ) .* ( log10( freq( ii ) ) - 2 );
NL( jj ) = 95 + sqrt( 21 * W ) - 17 *   log10( freq( jj ) );

