function U = scalep( Dalpha, c, r, U, RunType, TopOpt, freq )

% Scale the pressure field

% Compute scale factor
switch ( RunType(2:2) )
   case ( 'C' )
      const = -Dalpha * sqrt( freq ) / c;
   case ( 'R' )
      const = -Dalpha * sqrt( freq ) / c;
   otherwise
      const = -1.0;
end

% Volume attenuation
% This is no longer used since attenuation is now incorporated in the ray
% trace itself

% switch ( TopOpt( 4 : 4 ) )
%    case ( 'T' )   % Thorp
%       f2 = ( freq / 1000.0 )^2;   % frequency to kHz then squared
%       
%       % Original Thorp (1967) formula
%       %     alpha = 40.0 * f2 / ( 4100.0 + f2 ) + 0.1 * f2 / ( 1.0 + f2 );
%       %     alpha = alpha / 914.4;     % dB / m
%       %     alpha = alpha / 8.6858896; % Nepers / m
%       
%       % Updated formula from JKPS Eq. 1.34
%       alpha = 3.3d-3 + 0.11 * f2 / ( 1.0 + f2 ) + 44.0 * f2 / ( 4100.0 + f2 ) + 3d-4* f2;   % dB/km
%    case ( 'F' )   % Francois-Garrison (untested)
%       T  = 20
%       S  = 35
%       pH = 8
%       z  = 0
%       alpha = franc_garr( freq / 1000, T, S, pH, z )   % dB/km
%    otherwise
%       alpha = 0.0;
% end
% alpha = alpha / 8685.8896; % convert dB/km to Nepers / m

% For incoherent RunType, convert intensity to pressure
if ( RunType(1:1) ~= 'C' )
   U = sqrt( real( U ) );
end

% add in attenuation
for ir = 1 : length( r )
   if RunType(4:4) == 'X' % line source
      factor = -4 * sqrt( pi ) * const;
   else
      if ( r( ir ) ~= 0 )
         % factor = const * exp( -alpha * r( ir ) ) / sqrt( r( ir ) );
         factor = const / sqrt( r( ir ) );
      else
         factor = 0.0;
      end
   end

   U( :, ir ) = factor * U( :, ir );
end
