% Pekeris Turboelectroencabulator

% Main calling routine to evaluate the eigenvalues, k, for a Pekeris waveguide
% This version is used to test different estimates of imag( k )

% Mike Porter August 2013

clear all

global freq D c1 c2 rho1 rho2 realflag

freq = 250;   % frequency in Hz
D    = 100;    % depth in meters
c1   = 1500;   % sound speed in water
c2   = 1800;   % sound speed in lower halfspace
rho1 = 1.0;    % density in water
rho2 = 1.8;    % density in lower halfspace

omega  = 2 * pi * freq;
omega2 = omega^ 2;

xMax   = omega2 / c1^2 - 100 * eps;   % upper wavenumber limit (squared) (reduced to avoid branch pt.)
xMin   = omega2 / c2^2;   % lower
xMin = 0

%%
realflag = 1;
x_noloss = bisect_by_count( xMin, xMax );   % use bisection to get the roots
k_noloss = sqrt( x_noloss );   % convert them to wavenumbers

Iterations = 10;
AtUnit     = 'W';
c2_hold    = c2;
x          = x_noloss;
realflag   = 0;

figure

for alpha = 0.5 % : .5 : 5 % (loss in dB / wavelength )
   % c2 = c2 + alpha * 50
   c2 = crci( c2_hold, alpha, freq, AtUnit )   % make complex c based on attenuation
   
   x = Newton( x, @pekeris, Iterations );
   k = sqrt( x );
   
   plot( k, 'b-*' );
   hold all
   
   % Loop length formula
   gamma = sqrt( omega2 / c1^2 - k_noloss.^2 );
   Dloop = 2 * k_noloss * D ./ gamma;
   
   % Reflection coefficient for a homogeneous halfspace
   gamma1SQ = omega2 / c1^ 2 - x_noloss;
   gamma2SQ = omega2 / c2^ 2 - x_noloss;
   gamma1   = sqrt( -gamma1SQ );
   gamma2   = sqrt( -gamma2SQ );
   
   Refl = ( rho2 * gamma1 - gamma2 ) ./ ( rho2 * gamma1 + gamma2 );
   
   % figure( 3 ); plot( k_noloss, abs( Refl ) )
   
   dk = log( abs( Refl ) ) ./ Dloop;
   
   plot( k_noloss + i * dk, 'r-' ); hold on
   
   % perturbation formula
   dk = -2 * ( 1 - abs( Refl ) ) ./ ( 1 + abs( Refl ) ) ./ Dloop;

   plot( k_noloss + i * dk, 'k-' ); hold on
   
   % loop formula 2
   u      = sqrt( 2 / D ) *          sin( gamma * D );
   uprime = sqrt( 2 / D ) * gamma .* cos( gamma * D );

   Dloop = 4 * gamma .* k_noloss ./ ( ( gamma .* u ) .^2 + uprime .^2 );
   Dloop = 4 * gamma .* k_noloss ./ ( ( gamma .* u ) .^2 );

   dk = log( abs( Refl ) ) ./ Dloop;
   
   plot( k_noloss + i * dk, 'g-' ); hold on
   pause( 1 )
end

%%
% Following is a two stage version where the roots are first bracketted or
% isolated, then a second rood finder is used to refine them
% This is slower because bisect does not vectorize as easily
%
% [ xL, xR ] = bisect( xMin, xMax );
% [ sqrt( xL' ) sqrt( xR' ) ];
% x = bisect_final( xL, xR );
%


% call Matlab's root finder to refine the root
% This one is also slower because fzero is not vectorized
%
% for mode = 1 : length( xL )
%    x0 = [ xL( mode ), xR( mode ) ];   % bracketting interval
%    [ xx( mode ), fval, exitflag ] = fzero( @pekeris, x0 );
%    % disp( [ mode sqrt( xx( mode ) ) ] )
% end
