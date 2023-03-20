function crci = crci( c, alpha, freq, AttenUnit )

%     Convert real wave speed and attenuation to a single complex wave speed

%     5 cases:
%     N for Nepers/meter
%     M for dB/meter      (M for Meters)
%     F for dB/m-kHZ      (F for frequency dependent)
%     W for dB/wavelength (W for Wavelength)
%     Q for Q

%     T for Thorpe

global T S pH z_bar

omega = 2.0 * pi * freq;

% *** Convert to Nepers/m ***

alphaT = 0.0;

switch ( AttenUnit( 1 : 1 ) )
   case ( 'N' )   % Nepers/m
      alphaT = alpha;
   case ( 'M' )   % dB/meter
      alphaT = alpha / 8.6858896;
   case ( 'F' )   % dB/m-kHZ
      alphaT = alpha * freq / 8685.8896;
   case ( 'W' )   % dB/wavelength
      if ( c ~= 0.0 )
         alphaT = alpha * freq / ( 8.6858896 * c );
      end
   case ( 'Q' )
      if( c * alpha ~= 0.0 )
         alphaT = omega / ( 2.0 * c * alpha );
      end
end

% added volume attenuation
switch ( AttenUnit( 2 : 2 ) )
   case ( 'T' )   % Thorp
      f2     = ( freq / 1000.0 )^2;
      
      % Original Thorp (1967) formula
      %    alphaT = 40.0 * f2 / ( 4100.0 + f2 ) + 0.1 * f2 / ( 1.0 + f2 );
      %    alphaT = alphaT / 914.4;     % dB / m
      %    alphaT = alphaT / 8.6858896; % Nepers / m
      
      % Updated formula from JKPS Eq. 1.34
      Thorp  = 3.3d-3 + 0.11 * f2 / ( 1.0 + f2 ) + 44.0 * f2 / ( 4100.0 + f2 ) + 3d-4 * f2;   % dB/km
      Thorp  = Thorp / 8685.8896; % Nepers / m
      alphaT = alphaT + Thorp;
   case ( 'F' )   % Francois-Garrison
      %T      = 20; % need to set these 4 values
      %S      = 35;
      %pH     = 8;
      %z_bar  = 0;
      %disp( [ T, S, pH, z_bar ] )
      FG     = franc_garr( freq / 1000, T, S, pH, z_bar );   % dB/km
      FG     = FG / 8685.8896; % Nepers / m
      alphaT = alphaT + FG;
end

%     *** Convert Nepers/m to equivalent imaginary sound speed ***

alphaT = alphaT * c^2 / omega;
crci   = c + 1i * alphaT;
