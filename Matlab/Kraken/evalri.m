function p = evalri( phiS, phi, R, k, Opt, Comp )

% conversion from Fortran eval.f90
% mbp 4/2009
%
% Computes pressure field from modes
% Normalized to pressure of point source at 1 meter
% Returns the pressure matrix on a rectangular grid for a single source depth
%
% Opt = X     Cartesian   (X, z) coordinates
% Opt = R     Cylindrical (R, z) coordinates
% Opt = S     Scaled cylindrical (same as cylindrical but without the sqrt( r ) decay in range)
%
% R is the range in meters

% Receiver range offsets (for array tilt) have not been implemented

% For vertical component of displacement field, take finite difference in depth
% This is done using a backward difference.
% A centered difference or an FFT formula would be better

if ( Comp == 'V' )
   phidiff = diff( phi );   % should divide by diff( z ) also
   phi( 1: end - 1, : ) = phidiff;
   phi( end,        : ) = zeros( 1, size( phi, 2 ) ); % no derivative for last row
end

%phi = phi * diag( phiS, 0 );	% scale modes by phiS
phi = scalecol( double( phi ), phiS );
% form pressure field

% avoid singularity at R=0 by replacing it by a small number
Rt = R;
Rt( Rt == 0.0 ) = 1e-9;

m = length( k );
n = length( Rt );

switch Opt( 1 : 1 )
   case 'R'
      phase = spdiags( 1.0 ./ sqrt( k ), 0 , m, m ) * exp( -1i * k * Rt' ) * spdiags( realsqrt( 2 * pi ./ Rt ), 0, n, n );
   case 'X'
      k0 = real( k( 1 ) );   % this should be k0, the wavenumber at the source depth
      k0 = 2 * pi * 1000 / 1500
      factor = realsqrt( 2 * pi * k0 );   % using far-field approximation to Hankel function
      factor = 2 / besselh( 0, k0 );      % exact Hankel function at 1 m

      phase = spdiags( 1.0 ./       k  , 0 , m, m ) * exp( -1i * k * Rt' ) * factor;
   case 'S'
      phase = spdiags( 1.0 ./ sqrt( k ), 0 , m, m ) * exp( -1i * k * Rt' ) * realsqrt( 2 * pi );
end

% for horizontal component take derivative in range direction
% The following formula approximates the derivative, assuming the phase term
% (e^(-i k rr )) dominates
if ( Comp == 'H' )
   phase = diag( -1i .* sqrt( k ) ) * exp( -1i * k * Rt' ) * diag( realsqrt( 2 * pi ./ Rt ) );
end

if Opt( 4 : 4 ) ~= 'I'   % coherent sum
  p = phi * phase;

else  % incoherent sum
  nz = size( phi, 1 );
  nr = length( Rt );
  p  = NaN( nz, nr );

  % following should be vectorized
  for jr = 1 : nr
    for jz = 1 : nz
      p( jz, jr ) = sqrt( sum( abs( phi( jz, : ) .* phase( :, jr ).' ).^2 ) );
    end
  end
end