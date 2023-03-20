function [ Iteration, ErrorMessage ] = secant_method( x2, Tolerance, MAXIteration )

% Secant method

ErrorMessage = ' '
if ( Tolerance <= 0.0 )
    ErrorMessage = 'Non-positive tolerance specified'
    stop
end

x1 = x2 + 100.0 * Tolerance;
f1 = funct( x1 );

for Iteration = 1 : MAXIteration
    x0 = x1;
    f0 = f1;

    x1 = x2;
    f1 = funct( x1 );

    % ugly stuff to block overflows by forcing shift to be bounded
    cNum = f1 * ( x1 - x0 );
    cDen = f1 - f0;

    if ( abs( cNum ) >= abs( cDen * x1 ) )
        shift = 100.0 * Tolerance * abs( x1 );
    else
        shift = cNum / cDen;
    end

    x2 = x1 - shift;
    if ( abs( x2 - x1 ) + abs( x2 - x0 ) < Tolerance ) return
    end
end

ErrorMessage = ' *** FAILURE TO CONVERGE IN SECANT'
