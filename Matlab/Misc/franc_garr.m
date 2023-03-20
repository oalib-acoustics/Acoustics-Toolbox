
function alpha = franc_garr( f, T, S, pH, z )

% Francois Garrison formulas for attenuation
% Script provided by D. Jackson APL-UW

% mbp Feb. 2019
% some reformatting and mod to allow f to be a vector
% Verified using F-G Table IV

% alpha = attenuation   (dB/km)
% f     = frequency     (kHz)
% T     = temperature   (deg C)
% S     = salinity      (psu)
% pH    = 7 for neutral water
% z     = depth         (m)

c = 1412 + 3.21 * T + 1.19 * S + 0.0167 * z;

% Boric acid contribution
A1 = 8.86 / c * 10^( 0.78 * pH - 5 );
P1 = 1;
f1 = 2.8 * sqrt( S / 35 ) * 10^( 4 - 1245 / ( T + 273 ) );

% Magnesium sulfate contribution
A2 = 21.44 * S / c * ( 1 + 0.025 * T );
P2 = 1 - 1.37 * 10^-4 * z + 6.2 * 10^-9 * z^2;
f2 = 8.17 * 10^( 8 - 1990 / ( T + 273 ) ) / ( 1 + 0.0018 * ( S - 35 ) );

% Viscosity
P3 = 1 - 3.83 * 10^-5 * z + 4.9 * 10^-10 * z^2;
if T < 20
   A3 = 4.937 * 10^-4 -2.59  * 10^-5 * T + 9.11 * 10^-7 * T^2 - 1.5 * 10^-8  * T^3;
else
   A3 = 3.964 * 10^-4 -1.146 * 10^-5 * T + 1.45 * 10^-7 * T^2 - 6.5 * 10^-10 * T^3;
end

alpha = A1 * P1 * ( f1 * f.^2 ) ./ ( f1^2 + f.^2 ) + A2 * P2 * ( f2 * f.^2 ) ./ ( f2^2 + f.^2 ) + ...
        A3 * P3 * f.^2;