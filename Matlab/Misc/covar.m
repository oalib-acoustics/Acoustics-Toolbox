function covar

%     Computes synthesized data vectors with Gaussian random noise added
%     Averages data vector to from covariance matrix

clear read_modes_bin % to force rewind to beginning of mode file

filename = 'Envs_Like_MBPJOE+HB/000.5_-00.5.mod';
NoiseType = 'C';   % colored or white noise
SNRDB     = 200;
NTimes    = 30000;

Pos.r.z = [ 5 6 ];   % vector or receiver depths
NRD         = length( Pos.r.z );

SNR    = 10.0 ^ ( SNRDB / 10.0 );
dNoisy = zeros( NRD, 1 );

PhonePower = 1.0 / NRD;

% Zero out covariance matrix

Noise      = zeros( NRD, 1 );
Covariance = zeros( NRD, NRD );

% *** If making colored noise, need modal data ***

if ( NoiseType == 'C' )
   Pos.s.z = 6.0;  % noise sources are at 1 m
   
   [ Modes ] = read_modes( filename );
   M = Modes.M;
   
   % calculate mode values at source depth
   zs   = Pos.s.z;
   isd  = find( Modes.z >= zs );    % index of source depth
   isd  = isd( 1 );
   phiS = Modes.phi( isd, : ).';    % column vector with modal weights
   
   % calculate mode values at receiver depths
   irdvec = zeros( 1, length( Pos.r.z ) );
   for ii = 1 : length( Pos.r.z )
      zr  = Pos.r.z( ii );
      ird = find( Modes.z >= zr );	% index of source depth
      irdvec( ii ) = ird( 1 );
   end
   phiR = Modes.phi( irdvec, : ).';
   
   if ( any( imag( Modes.k( 1 : M ) ) >= 0.0 ) )
      error( 'Modes present with no loss: Noise field is infinite ' )
   end
   
   % Construct phi and Power
   
   for mode = 1 : M
      phiR( mode, : ) = phiS( mode ) * phiR( mode, : ) / ...
         sqrt( -real( Modes.k( mode ) ) * imag( Modes.k( mode ) ) );
      Power = norm( phiR( mode, : ) );
   end
end

% *** Main loop: generate realizations of noise field ***

% AverageSignal = 0.0;
AverageNoise  = 0.0;
sigma         = sqrt( PhonePower / SNR );   % Noise level on a phone

for iTime = 1 : NTimes
   
   switch NoiseType( 1 : 1 )
      case 'W'
         Noise = White( sigma, NRD );
      case'C'
         Noise = Color( phiR, M, NRD, sigma, Power );
      otherwise
         disp( 'Unknown noise type' )
   end
   
   dNoisy = Noise;
   
   % AverageSignal = sum( abs( d     ) ^ 2 );
   AverageNoise = sum( abs( Noise ) .^ 2 );
   
   Covariance = Covariance + dNoisy * dNoisy';   % cross-sensor correlation matrix R
   
end

% AverageSignal = AverageSignal / ( NTimes * NRD )
AverageNoise  = AverageNoise  / ( NTimes * NRD )

% disp( 'Average signal (dB) = ' )
% disp( 10.0 * log10( AverageSignal ) )
disp( 'Average noise  (dB) = ' )
disp( 10.0 * log10( AverageNoise  ) )
% disp( 'Average SNR    (dB) = ' )
% disp( 10.0 * log10( AverageSignal / AverageNoise ) )

% *** Analytic formula for covariance matrix in isovelocity ocean

CovarianceISO = zeros( NRD, NRD );

for I = 1 : NRD
   for J = 1 : NRD
      Depth = 100.0;
      for mode = 1 : M
         gamma = ( mode - .5 ) * pi / Depth;
         CovarianceISO( I, J ) = CovarianceISO( I, J ) + sin( gamma * Pos.s.z( 1 ) ) ^ 2 * ...
            sin( gamma * Pos.r.z( I ) ) * sin( gamma * Pos.r.z( J ) );
      end
   end
end

% *** normalize cross-sensor correlation matrix

Covariance = NRD * Covariance / ( ( 1.0 + sigma^2 ) * NTimes )
%  ABS( CovarianceISO( I, J ) * Covariance( 1, 1 ) / CovarianceISO( 1, 1 ) )

end

%----------------------------------------------------------------------

function Noise = White( sigma, NRD )

% Generates Gaussian random white noise

Noise = zeros( 1 : NRD, 1 );  % Zero out noise vector

X = rand( NRD, 1 );
Y = rand( NRD, 1 );

X = X + 0.000001;

for IRD = 1 : NRD
   Noise( IRD ) = sigma * sqrt( -log( X( IRD ) ) ) * exp( 1i * 2.0 * pi * Y( IRD ) );
end

% Note power = NRD * sigma ^ 2

end

%----------------------------------------------------------------------

function Noise = Color( phi, M, NRD, sigma, Power )

% Generates colored noise vector

Noise = zeros( NRD, 1 );  % Zero out noise vector

% compute the random coefficient
X = rand( M, 1 );
Y = rand( M, 1 );
X = X + 0.000001;

Coef = sqrt( -log( X ) ) .* exp( 1i * 2.0 * pi * Y );

% modal sum for each receiver
Noise = Noise + ( Coef.' * phi ).';

% Normalize so that power = NRD * sigma^2

Noise = sigma * Noise * sqrt( NRD / Power );

end
