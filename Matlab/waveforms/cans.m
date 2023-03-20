function [ S, PulseTitle ] = cans( t, omega, Pulse )

% Computes the source time series
%
%     t      is the time (can be a scalar or a vector)
%     omega  is some frequency characterizing the pulse
%     Pulse  is the letter indicating which pulse
%     S      is the time series (takes same shape as t)
%
% note: if you're going to form a Hilbert transform of these,
% need to evaluate pulses at negative time to avoid artifacts

% mbp, based on the original 1988 Fortran version from SPARC

S          = zeros( size( t ) );   % initialize to zero
PulseTitle = Pulse; % must return this, even if T<=0
F          = omega / ( 2.0 * pi );

%     Evaluate selected Pulse

switch ( Pulse(1:1 ) )
   case ( 'P' )   % Pseudo gaussian
      % (peak at 0, support [0, 3F]
      ii = find( t > 0 & t <= 1 / F );
      T = t( ii );
      
      S( ii ) = 0.75 - cos( omega * T ) + 0.25 * cos( 2.0 * omega * T );
      PulseTitle = 'Pseudo gaussian';
      
   case ( 'R' )   % Ricker wavelet
      % (peak at F, support [0, 2F]
      ii = find( t > 0 );
      T = t( ii );
      
      U = omega * T - 5.0;
      S( ii ) = 0.5 * ( 0.25 * U .* U - 0.5 ) * sqrt( pi ) .* exp( -0.25 * U .* U );
      PulseTitle = 'Ricker wavelet';
      
   case ( 'A' )    % Approximate Ricker wavelet
      % (peak at F, support [0, 2.5F]
      TC = 1.55 / F;
      ii = find( t > 0 & t <= TC );
      T = t( ii );
      
      S( ii ) = ...
         +0.48829 *     cos( 2.0 * pi * T / TC ) ...
         -0.14128 * 4 * cos( 4.0 * pi * T / TC ) ...
         +0.01168 * 9 * cos( 6.0 * pi * T / TC );
      
      PulseTitle = 'Approximate Ricker wavelet';
      
   case ( 'S' )   % Single sine
      % (peak at F, support [0, infinity], nulls at nF
      ii = find( t > 0 & t <= 1 / F );
      T = t( ii );
      
      S( ii ) = sin( omega * T );
      
      PulseTitle = 'Single sine';
      
   case ( 'H' )   % Hanning weighted four sine
      % (peak at F, support [0, infinity], first null at about 1.5F
      ii = find( t > 0 & t <= 4 / F );
      T = t( ii );
      
      S( ii ) = 0.5 * sin( omega * T ) .* ( 1 - cos( omega * T / 4.0 ) );
      
      PulseTitle = 'Hanning weighted four sine';
      
   case ( 'N' )   % N-wave
      % (peak at 1 F, support [0, 4F], [0,3F] OK
      ii = find( t > 0 & t <= 1 / F );
      T = t( ii );
      
      S( ii ) = sin( omega * T ) - 0.5 * sin( 2.0 * omega * T );
      
      PulseTitle = 'N-wave';
      
   case ( 'M' )   % Miracle wave (has a Hilbert transform?)
      % (peak at 0, support [0, infinity]
      ii = find( t > 0 );
      T = t( ii );
      A  = 1.0 / ( 6.0 * F );
      T0 = 6.0 * A;
      TS = ( T - T0 ) / A;
      S( ii )  = 1.0 ./ ( 1.0 + TS .* TS );
      % HS is the Hilbert transform of the time series
      % HS = TS / ( 1.0 + TS * TS );
      PulseTitle = 'Miracle wave';
      
   case ( 'G' )   % Gaussian
      % (peak at 0, support [0, infinity]
      % T0 is the peak time
      % A is the 3dB down time
      % NOTE S(0) = exp( -NSIG ** 2 )
      ii = find( t > 0 );
      T = t( ii );
      
      NSIG = 3;
      A  = 1.0 / F / ( 2.0 * NSIG );
      T0 = NSIG * A;
      S( ii )  = exp( -( ( T - T0 ) / A ) .^ 2 );
      PulseTitle = 'Gaussian';
   case ( 'T' )   % Tone burst (gated sinewave)
      % (peak at F, support [0, infinity]
      ii = find( t > 0 & t <= 0.4 );
      T = t( ii );
      
      S( ii ) = sin( omega * T );
      
      PulseTitle = 'Tone';
   case ( 'C' )   % Sinc
      % ( uniform spectrum from [0, F]
      
      T = t( ii );
      S( ii ) = sin( omega * T ) / ( omega * T );
      PulseTitle = 'Sinc';
   case( 'E' )   % One-sided exponential
      ii = find( t > 0 );
      T = t( ii );
      S( ii ) = exp( -omega * T );
      PulseTitle = 'One-sided exponential';
end
