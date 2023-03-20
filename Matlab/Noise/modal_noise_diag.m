function NL = modal_noise_diag( ModeFile, rho_SL_dB, sd, rd, freq, Rmax_km, Component )

% Computes synthesized data vectors with Gaussian random noise added
% Averages data vector to from covariance matrix
%
% rho_SL_dB source level density (dB expressed per m^2 of the noise sheet)
% ModeFile  name of the mode file
% sd        source depth (m)
% Rmax_km   radius (km) of the noise disc (can be a vector)
% rd        receiver depth (m) (can be a vector)
% Component 'P', 'V', or 'H' for pressure, vertical, or horizontal
%
% returns
% NL in dB with nrd rows and nrr columns
%
% mike porter, 2012

% This is based on Eq, 9.27 from JKPS
% That equation should be corrected when the eigenfunctions are complex
% Also, the assumption that off-diagonal terms can be dropped is not always valid

clear read_modes_bin % to force rewind to beginning of mode file

Rmax        = 1000. * Rmax_km;
Pos.s.z = sd;  % noise sources are at 1 m
Pos.r.z = rd;  % vector of receiver depths (need two values for dphi/dz

Nrd         = length( Pos.r.z );
Nrr         = length( Rmax );

rho_SL = 10^( rho_SL_dB / 10 ); % source level per m^2

% *** Read modal data ***

[ Modes ] = read_modes( ModeFile, freq );
M = Modes.M;

if ( M > 0 )
   % calculate mode values at source depth
   phiS = interp1_mbp( Modes.z, Modes.phi, Pos.s.z );
   
   % calculate mode values at receiver depths
   % Modes.phi is a matrix nrd x nmodes
   % phiR should be a matrix nmodes x nrd
   
   phiR = interp1_mbp( Modes.z, Modes.phi, Pos.r.z ).';
   
   % calculated k0 used as a scaling factor in Eq. 9.17 JKPS
   % Becomes irrelevant because it goes into q2, but q2 is divided by k0 in
   % the final noise level
   omega = 2 * pi * Modes.freqVec( 1 );
   rho   = 1000;       % s.i. units (kg/m^3)
   c     = 1480;       % nominal sound speed
   k0    = omega / c;  % nominal wavenumber at the source depth
   
   % Here's the factor to convert a SL density to q2
   q2 = k0^2 * 4 * pi * rho_SL;
   
   % compute scaled mode, phiScaled, depending on the component
   % it is either the mode or something involving a derivative
   
   phiScaled = zeros( size( phiR ) );
   
   switch Component
      case 'P' % pressure
         phiScaled = conj( phiR );
         
      case 'V' % vertical velocity formula (do a depth derivative)
         disp( 'Check dz used for depth derivative' )
         dz          = 10;   % used for the depth derivative when vertical veleocity is selected
         
         for mode = 1 : M
            % [ abs( ( phiR( mode, 2 : end ) - phiR( mode, 1 : end-1) ) / abs( phiR( mode, 1 ) ) ) / ...
            %   ( omega / c )  abs( Modes.k( mode ) ) / ( omega / c ) ];
            % add a zero on the end to preserve length of phiScaled
            phiScaled( mode, : ) = ( phiR( mode, 2 : end ) - phiR( mode, 1 : end-1 ) ) / ( dz * rho * omega );
         end
         
      case 'H' % horizontal velocity formula
         phiScaled = diag( Modes.k ) * phiR / ( 2^(1/4) * rho * omega );
         
   end
   
   % any growing modes?
   if ( any( imag( Modes.k( 1 : M ) ) >= 0.0 ) )
      disp( 'Modes present with no loss: Noise field is infinite ' )
      %find( imag( Modes.k( 1 : M ) ) >= 0.0 )
   end
   
   ii = find( imag( Modes.k( 1 : M ) ) < 0 )';   % these are the ones we keep
   
   % *** Construct phi and Power  (Eq. 9.22 and 9.25 jkps) ***
   
   phiSR     = zeros( M, Nrd );
   Intensity = zeros( Nrd, Nrr );
   
   for ir = 1 : length( Rmax )
      for mode = ii
         % kr   = real( Modes.k( mode ) );
         k    = Modes.k( mode );
         ki   = imag( k );
         kabs = abs(  k );
         
         % pressure formula using Appendix J of Carey and Evans
         % original
         %phiSR( mode, : ) = ( abs( phiS( mode ) * phiScaled( mode, : ) ) / ...
         %   sqrt( -kabs .* ki ) ).^2 * ( 1. - exp( 2 * ki * Rmax( ir ) ) );
         
         % equivalent form
         phiSR( mode, : ) = ( abs( phiS( mode ) * phiScaled( mode, : ) ) ).^2 / ...
            ( kabs * ki ) * ( 1. - exp( 2 * ki * Rmax( ir ) ) );
         
         % previous forms assume small ki; this one fixes that
         phiSR( mode, : ) = -4 * 1i * ( abs( phiS( mode ) * phiScaled( mode, : ) ) ).^2 / ...
            ( k^2 - conj( k )^2 ) * ( 1. - exp( 2 * ki * Rmax( ir ) ) );
         % this one works
         phiSR( mode, : ) = ( phiS( mode ) * phiScaled( mode, : ) ).^2 / ...
            ( -kabs .* ki ) * ( 1. - exp( 2 * ki * Rmax( ir ) ) );
         
         % this one also works
         % phiSR( mode, : ) = ( abs( phiS( mode ) ) * phiScaled( mode, : ) ) .^2 / ...
         %   ( -kabs .* ki ) * ( 1. - exp( 2 * ki * Rmax( ir ) ) );
      end
      
      % modal sum for each receiver
      dNoisy = squeeze( sum( phiSR, 1 ) );
      dNoisy = dNoisy.';
      
      % I removed dNoisy' and put a square in pressure formula above ...
      Covariance         = ( ( pi / 2 ) * q2 / k0^2 ) * dNoisy; % * dNoisy'   % cross-sensor correlation matrix R
      Intensity( :, ir ) = abs( Covariance( :, 1 ) );
   end
   
   switch Component
      case 'P'
         NL = Intensity;
      case 'V'  % vertical   velocity
         NL = ( rho * c )^2 * Intensity;
      case 'H'  % horizontal velocity
         NL = ( rho * c )^2 * Intensity;
   end
   
   NL = 10 * log10( NL );   % convert to dB
   
else
   NL = NaN;   % return a NaN if there are no modes
end
% *** Analytic formula for covariance matrix in isovelocity ocean

% CovarianceISO = zeros( NRD, NRD );
%
% for I = 1 : NRD
%    for J = 1 : NRD
%       Depth = 100.0;
%       for mode = 1 : M
%          gamma = ( mode - .5 ) * pi / Depth;
%          CovarianceISO( I, J ) = CovarianceISO( I, J ) + sin( gamma * Pos.s.z( 1 ) ) ^ 2 * ...
%             sin( gamma * Pos.r.z( I ) ) * sin( gamma * Pos.r.z( J ) );
%       end
%    end
% end
