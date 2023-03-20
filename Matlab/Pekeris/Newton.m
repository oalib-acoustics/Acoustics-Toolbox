function x = Newton( x0, f, Iterations )

% August 2013, Mike Porter
% Implements Newton's method for finding roots of a function
%
% x0 is an initial guess
% f  is the handle to the function
% x  is the approximation to the root
% Iterations is the number of iterations done

x = x0;

for ii = 1 : Iterations
   h      = 100 * eps * x;
   fprime = ( f( x + h ) - f( x ) ) ./ h;
   x = x - f( x ) ./ fprime;
end

