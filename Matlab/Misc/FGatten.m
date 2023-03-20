function alpha = FGatten( f, OceanSea ) %  S, T, pH, z )

% Francois-Garrison/Ainslie-McColm attenuation as a function of
% f = frequency (kHz) (can be a vector)
% S = salinity (ppt)
% T = temperature (degrees Celsius)
% pH
% z = depth (km)
% See Ainslie and McColm, JASA 103(3):1671- 1998

switch OceanSea
    case 'Pacific'
        pH = 7.7;
        S  = 34;
        T  = 4;
        z  = 1;
    case 'Red'
        pH = 8.2;
        S  = 40;
        T  = 22;
        z  = 0.2;
    case 'Arctic'
        pH = 8.2;
        S  = 30;
        T  = -1.5;
        z  = 0;
    case 'Baltic'
        pH = 7.9;
        S  = 8;
        T  = 4;
        z  = 0;
    otherwise
        disp( 'Unknown ocean or sea' )
end

f1 = 0.78 * sqrt( S / 35 ) * exp( T / 26 );   % boron
f2 = 42 * exp( T / 17 );   % magnesium

% attenuation in dB/km:
alpha = 0.106 * f1 * f.^2 ./ ( f.^2 + f1^2 ) * exp( ( pH - 8 ) / 0.56 ) + ...
    + 0.52 * ( 1 + T / 43 ) * ( S / 35 ) * f2 * f.^2 ./ ( f.^2 + f2^2 ) * exp( -z / 6 ) + ...
    + 0.00049 * f.^2 * exp( - ( T /27 + z / 17 ) );
