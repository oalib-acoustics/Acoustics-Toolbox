function NL = modal_noise_full( ModeFile, rho_SL_dB, sd, rd, freq, Rmax_km, Component )

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
% Cov        covariance matrix for the noise with nrd rows and columns
%
% mike porter, 2012

% This is based on Eq, 9.27 from JKPS
% That equation should be corrected when the eigenfunctions are complex
% See also Carey and Evans appendix

% The Component selection in this code is NOT implemented correctly
% but I'm leaving the code in as it is a good basis for that

clear read_modes_bin % to force rewind to beginning of mode file

Rmax        = 1000. * Rmax_km;
Pos.s.z = sd;  % noise sources are at 1 m
Pos.r.z = rd;  % vector of receiver depths (need two values for dphi/dz

dz = 10;
Nrd         = length( Pos.r.z );
Nrr         = length( Rmax );

rho_SL = 10^( rho_SL_dB / 10 ); % source level per m^2

% *** Read modal data ***

[ Modes ] = read_modes( ModeFile, freq );
M = Modes.M;

if ( M > 0 )
   
   % calculate mode values at source depth
   phiS = interp1( Modes.z, Modes.phi, Pos.s.z );
   phiS = phiS.';   % make it a column vector
   
   % calculate mode values at receiver depths
   % Modes.phi is a matrix nrd x nmodes
   % phiR should be a matrix nmodes x nrd
   
   phiR = interp1( Modes.z, Modes.phi, Pos.r.z ).';
   
   % calculated k0 used as a scaling factor in Eq. 9.17 JKPS
   % Becomes irrelevant because it goes into q2, but q2 is divided by k0 in
   % the final noise level
   omega = 2 * pi * Modes.freqVec( 1 );
   rho   = 1000;       % s.i. units (kg/m^3)
   c     = 1480;       % nominal sound speed
   k0    = omega / c;  % nominal wavenumber at the source depth
   
   % Here's the factor to convert a SL density to q2
   q2 = k0^2 * 4 * pi * rho_SL;
   
   % compute scaled mode, phiScale, depending on the component
   % it is either the mode or something involving a derivative
   
   switch Component
      case 'P' % pressure
         phiR_Scaled = phiR;
         
      case 'V' % vertical velocity formula (do a depth derivative)
         for mm = 1 : M
            phiR_Scaled( mm, : ) = [ abs( ( phiR( mm, 2 : end ) - phiR( mm, 1 : end - 1) ) / abs( phiR( mm, 1 ) ) ) / ...
               ( omega / c )  abs( Modes.k( mm ) ) / ( omega / c ) ];
            % add a zero on the end to preserve length of phiScaled
            phiR_Scaled( mm, : ) = ( phiR( mm, 2 : end ) - phiR( mm, 1 : end-1 ) ) / ( dz * rho * omega );
         end
         
      case 'H' % horizontal velocity formula
         phiR_Scaled = diag( Modes.k ) * phiR / ( 2^(1/4) * rho * omega );
         
   end
   
   % any growing modes?
   if ( any( imag( Modes.k( 1 : M ) ) >= 0.0 ) )
      disp( 'Modes present with no loss: Noise field is infinite ' )
      %find( imag( Modes.k( 1 : M ) ) >= 0.0 )
   end
   
   %%% ii = find( imag( Modes.k( 1 : M ) ) < 0 )';   % these are the ones we keep
   
   % *** Construct phi and Power  (Eq. 9.22 and 9.25 jkps) ***
   
   C         = zeros( M, Nrd );
   Intensity = zeros( Nrd, Nrr );
   
   % A = zeros( M, M );
   
   % loop over radii of noise discs
   
   kmat = zeros( M, M );
   for ii = 1 : M
      kmat( :, ii ) = conj( Modes.k( ii )^2 - conj( Modes.k( : ) ).^2 );
   end
   
   for ir = 1 : length( Rmax )
      % pressure formula using Appendix I of Carey and Evans
      
      for mm = 1 : Nrd
         for nn = 1 : Nrd
            
            % set up the A matrix
            %for ii = 1 : M
            %for jj = 1 : M;
            %   A( ii, jj ) = 4 * phiR_Scaled( ii, mm ) * conj( phiR_Scaled( jj, nn ) ) / ...
            %       ( Modes.k( ii )^2 - conj( Modes.k( jj ) )^2 ); % * ( 1. - exp( ( Modes.k( ii ) - conj( Modes.k( jj ) ) ) * Rmax( ir ) ) );
            %end
            %A( :, ii ) = 4 * conj( phiR_Scaled( ii, mm ) ) * phiR_Scaled( :, nn ) ./ ...
            %       kmat( ii, : ).'; % * ( 1. - exp( ( Modes.k( ii ) - conj( Modes.k( : ) ) ) * Rmax( ir ) ) );
            %end
            A = 4 * phiR_Scaled( :, nn ) * ( phiR_Scaled( :, mm ) )' ./ kmat;
            C( mm, nn ) = phiS' * A * phiS;
         end
      end
      
      dNoisy = diag( C );
      
      % I removed dNoisy' and put a square in pressure formula above ...
      Covariance         = ( ( pi / 2 ) * q2 / k0^2 ) * dNoisy;
      Intensity( :, ir ) = abs( Covariance );
   end
   
   switch Component
      case 'P'
         NL = Intensity;
      case 'V' % vertical velocity
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
