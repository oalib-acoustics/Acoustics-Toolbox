function C = noisecorr( zprime, R, z1, z2, filename )

% useage:
%   C = noisecorr( R, z1, z2, filename )
%
% program to compute correlation of the noise field using normal modes
%
% zprime depth of the sheet of noise sources
% z1        depth of the first phone
% z2        depth of the second phone (can be a vector)
% R         separation in range between the two receivers
% filename  mode file created by KRAKEN
% modes     optional vector of modes to include (not implemented yet)

% mike porter 12/2009 based on Eq. 9.25 in JKPS

rho = 1;   % density
q   = 1;   % noise source strength

if nargin == 5
  [ Modes ] = read_modes( filename );
else
  [ Modes ] = read_modes( filename, modes );
end
M = length( Modes.k );
k = 2 * pi * Modes.freq / 1500.0;   % should really be the wavenumber at z' !!!

% weight modes by mode excitation
% warning: should be doing linear interpolation !!!

idzprime = find( Modes.z >= zprime );		% index of first phone
phiprime = Modes.phi( idzprime( 1 ), : ).';   % column vector with modal weights

id = find( Modes.z >= z1 );		% index of first phone
phi1   = Modes.phi( id( 1 ), : ).';   % column vector with modal weights

% mode values at second phone
idvec = zeros( 1, length( z2 ) );
for ii = 1 : length( z2 )
   id = find( Modes.z >= z2( ii ) );		% index of second phone
   idvec( ii ) = id( 1 );
end
phi2 = Modes.phi( idvec, : );

% loop over modes to compute the correlation

% ok; here I'm not allowing z2 to be a vector !!!
C = 0;
for mode = 1 : M
   kappa = real( Modes.k( mode ) );
   alpha = -imag( Modes.k( mode ) );
   C = C + phiprime( mode ) ^2 * phi1( mode ) * phi2( 1, mode ) * besselj( 0, kappa * R ) / ( alpha * kappa );
end

C = pi * q^2 / ( 2 * rho^2 * k^2 ) * C;   % scale factor




