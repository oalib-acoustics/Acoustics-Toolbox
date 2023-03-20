function [ f, ModeCount ] = pekeris( k2 )

% returns the characteristic function for the Pekeris waveguide
% k is the horizontal wavenumber (row or column vector)
% returns f, the characteristic function

% note: characteristic function vanishes at k = omega / c1

% Mike Porter 11/2009

global freq D c1 c2 rho1 rho2 realflag

omega2 = ( 2 * pi * freq ) ^ 2;

kz1 = sqrt( omega2 / c1^2 - k2 );
% kz2 = sqrt( k2 - omega2 / c2^2 );
kz2 = PekerisRoot( k2 - omega2 / c2^2 ).';

% Here is the Pekeris characteristic function
% multiplied by 1i and rearranged a bit

f = rho2 * kz1 .* cos( kz1 * D ) + rho1 * kz2 .* sin( kz1 * D );   % formula in JKPS, Ch. 2

if realflag
   f = real( f );
end

% calculate the mode count for the trial eigenvalue, k2
ModeCount = floor( kz1 * D / pi );

ii = find( sign( f .* sin( kz1 * D ) ) < 0 );
ModeCount( ii ) = ModeCount( ii ) + 1;

%%
% plotting ...

% fz1 = [ real( kz1 ) ; imag( kz1 ) ];
% fz2 = [ real( kz2 ) ; imag( kz2 ) ];
% ff = [ real( f ) ; imag( f ) ];
%
% figure; plot( k, fz1 )
% figure; plot( k, fz2 )
% figure; plot( k, ff  )
%
% axis( [ min( k ), max( k ), -10, 10 ] )