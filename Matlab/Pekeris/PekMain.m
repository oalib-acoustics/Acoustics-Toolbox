% Pekeris Turboelectroencabulator

% Main calling routine to evaluate the eigenvalues, k, for a Pekeris waveguide

% Mike Porter 11/2009

global freq D c1 c2 rho1 rho2

freq = 250;   % frequency in Hz
D    = 100;    % depth in meters
c1   = 1500;   % sound speed in water
c2   = 1700;   % sound speed in lower halfspace
rho1 = 1.0;    % density in water
rho2 = 2.0;    % density in lower halfspace


freq = 10;   % frequency in Hz
D    = 600;    % depth in meters
c1   = 1500;   % sound speed in water
c2   = 1553.9;   % sound speed in lower halfspace
rho1 = 1.0;    % density in water
rho2 = 2.0;    % density in lower halfspace

omega2 = ( 2 * pi * freq ) ^ 2;
xMax   = omega2 / c1^2;   % upper wavenumber limit (squared)
xMin   = omega2 / c2^2;   % lower

xMax = xMax - 100 * eps;   % reduce xMax so that it doesn't lie on the branch point
xMin = 0

x    = bisect_by_count( xMin, xMax );   % use bisection to get the roots
k    = sqrt( x );   % convert them to wavenumbers

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
