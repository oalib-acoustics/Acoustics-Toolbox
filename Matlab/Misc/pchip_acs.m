function rv = pchip_acs(x, y, xx)

% PCHIP_ACS - a variant of monotone PCHIP developed at HLS Research, Inc.
%
% Usage: rv = pchip_acs(x, y, xx)
%
% Inputs:
%    x - array of strictly increasing abscissa data values
%    y - array of dependent values associates with the above
%   xx - optional array of points where the polynomial will be evaluated
%
% Output:
%   rv - return value:
%         a piecewise polynomial struct when argument xx is absent, otherwise
%         array with the values of the interpolant evaluated at xx
%
% The resulting piecewise polynomial respects the monotonicity of the
% input data (on any interval where the y are strictly increasing or
% decreasing, the interpolating piecewise polynomial will be too). As
% with all versions of monotone piecewise cubic Hermite interpolation,
% the first derivative is continuous everywhere.
%
% This is a new variant of monotone piecewise cubic Hermite interpolation
% that attempts to provide an additional degree of smoothness. First, the
% cubic spline that interpolates the data points is found using specified
% derivatives at the first and last point (estimated using a non-centered
% three point difference formula). At the knots where the derivatives from
% the cubic spline interpolant don't conflict with the monotonicity
% requirement, the 2-nd derivative will be continuous. When the derivatives
% from the cubic spline conflict with the monotonicity requirement, they are
% adjusted accordingly to provide monotonicity, but the 2-nd derivative will
% fail to be continuous at that data point. Given that the interpolant tries
% to be as smooth as the spline, but sometimes fails, we have called our new
% algorithm the monotone PCHIP ACS (almost a cubic spline) interpolant.
%
% $Id: pchip_acs.m,v 1.3 2019/08/21 22:32:07 jcp Exp $

% Check for missing input arguments

if nargin < 2
  error( [mfilename, ': one or more required input arguments is missing.'] );
end

if nargin > 3
  error( [mfilename, ': too many input arguments, expecting at most three.'] );
end

% Some further sanity checks on the input

npts = length(x);

if npts ~= length(y)
  error( [mfilename, ': lengths of x and y arrays are not equal.'] );
end

if npts < 2
  error( [mfilename, ': there must be 2 or more data points.'] );
end

if ~isreal(x)
  error( [mfilename, ': argument x must be of type real.'] );
end

if min(diff(x)) <= 0.0
  error( [mfilename, ': x values are not strictly monotone increasing.'] );
end

% Check for complex valued ordinate

if ~isreal(y)
  if nargin == 3
    rv_r = pchip_acs(x, real(y), xx);
    rv_i = pchip_acs(x, imag(y), xx);
    rv = complex(rv_r, rv_i);
  else
    [ breaks, coefs_r, ~, ~, dim ] = pchip_acs(x, real(y));
    [ ~     , coefs_i, ~, ~, ~   ] = pchip_acs(x, imag(y));
    rv = mkpp(breaks, complex(coefs_r, coefs_i), dim);
  end
  return
end

% Handle the special case of exactly two data points separately

if npts == 2

% Case of exactly two data points degenerates to linear intepolation % % % % % %

coefs(1, 1) = 0.0;
coefs(1, 2) = 0.0;
coefs(1, 3) = (y(2) - y(1)) / (x(2) - x(1));
coefs(1, 4) = y(1);

else

% Case of more than two data points  % % % % % % % % % % % % % % % % % % % % % %

% Allocate memory for the components of the piecewise polynomial

coefs = zeros(npts-1, 4);

fprime = zeros(npts, 1);

% Estimate the derivative at the first node, using a non-centered three
%   point difference formula (see the Appendix of UCRL-85104 for details)

h1 = x(2) - x(1);
h2 = x(3) - x(2);
del1 = (y(2) - y(1)) / h1;
del2 = (y(3) - y(2)) / h2;

f1prime = ( (2.0*h1+h2)*del1 - (h1)*del2 ) / (h1 + h2);

if sign(f1prime) ~= sign(del1)
  % set to zero if the signs don't agree
  f1prime = 0.0;
elseif (sign(del1) ~= sign(del2)) && (abs(f1prime) > abs(3.0*del1))
  % enforce monotonicity
  f1prime = 3.0*del1;
end

fprime(1) = f1prime;

% Estimate the derivative at the last node, using a non-centered three
%   point difference formula (see the Appendix of UCRL-85104 for details)

h1 = x(end-1) - x(end-2);
h2 = x(end)   - x(end-1);
del1 = (y(end-1) - y(end-2)) / h1;
del2 = (y(end)   - y(end-1)) / h2;

f2prime = ( (-h2)*del1 + (h1+2.0*h2)*del2 ) / (h1 + h2);

if sign(f2prime) ~= sign(del2)
  % set to zero if the signs don't agree
  f2prime = 0.0;
elseif (sign(del2) ~= sign(del1)) && (abs(f2prime) > abs(3.0*del2))
  % enforce monotonicity
  f2prime = 3.0*del2;
end

fprime(end) = f2prime;

% Compute cubic spline coefficients (to compute estimates of 1-st derivative)

pp_spline = spline(x(:), [fprime(1); y(:); fprime(end)]);
[ breaks_s, coefs_s, pieces_s, order_s, dim_s ] = unmkpp(pp_spline);
coefs_ydot = repmat(order_s-1:-1:1, pieces_s, 1) .* coefs_s(:, 1:order_s-1);
pp_dot = mkpp(breaks_s, coefs_ydot, dim_s); % reassemble the differentiated pp

% Determine estimates of the derivative at each interior node, knot

del_j = (y(2) - y(1)) / (x(2) - x(1));

for j = 2:npts-1

  % save the secant, slope of the previous interval

  del_jm1 = del_j;

  % compute the secant, slope of this interval

  del_j = (y(j+1) - y(j)) / (x(j+1) - x(j));

  % determine an estimate the derivative at this knot that ensures monotonicity

  fprime(j) = G(del_jm1, del_j, ppval(pp_dot, x(j)));

end

% Compute the coefficients of the interpolating polynomial on each interval

for j = 1:npts-1

  f1 = y(j);
  f2 = y(j+1);

  f1prime = fprime(j);
  f2prime = fprime(j+1);

  % compute coefficients of the standard Hermite polynomial

  dx = x(j+1) - x(j);

  coefs(j, 1) = (f1prime + f2prime)/dx^2 - 2.0*(f2 - f1)/dx^3;
  coefs(j, 2) = 3.0*(f2-f1)/dx^2 - (2.0*f1prime + f2prime)/dx;
  coefs(j, 3) = f1prime;
  coefs(j, 4) = f1;

end

end	% case of more than two data points

% Create a piecewise polynomial struct

pp = mkpp(x, coefs, 1);

% Check if the caller wants to evaluate the piecewise polynomial

if nargin == 3
  rv = ppval(pp, xx);
else
  rv = pp;
end

return

end

% Define the G helper function 

function dydx = G(del1, del2, fprime)

% Check if the derivative is within the trust region, project into it if not

if (del1*del2) > 0.0
  if sign(del1) > 0.0
    dydx = min(max(fprime, 0.0), 3.0*min(del1, del2));
  else
    dydx = max(min(fprime, 0.0), 3.0*max(del1, del2));
  end
else
  % force the interpolant to have an extrema here
  dydx = 0.0;
end

return

end

%
% End of pchip_acs.m
