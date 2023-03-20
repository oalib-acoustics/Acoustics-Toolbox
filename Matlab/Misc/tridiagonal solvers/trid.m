% MATH5315: High Performance Numerical Computing
% File:     trid.m
% Type:     Matlab function
% Requires: -
% Usage:    x = trid(a, b, c, d)
% Purpose:  Solution of tridiagonal linear system T*x = d where
%           T = diag(a,-1) + diag(b,0) + diag(c,1)
%           using  standard sequential factorization and substitution
% Example:  See Matlab script tridtst.m
function x = trid(a, b, c, d)

n = length(b);
x = zeros(n,1);

% Gaussian elimination
for k = 1:n-1

   p = a(k)/b(k);
   b(k+1) = b(k+1) - p*c(k);
   d(k+1) = d(k+1) - p*d(k);

end

% backsubstitution
x(n) = d(n)/b(n);
for k = n-1:-1:1

   x(k) = (d(k) - c(k)*x(k+1)) / b(k);

end

