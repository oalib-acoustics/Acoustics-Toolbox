function out = equally_spaced( x )

% Test whether the vector x is composed of roughly equally-spaced values

n     = length( x );
xtemp = linspace( x( 1 ), x( end ), n );

delta = abs( x( : ) - xtemp( : ) );   % the colon converts both to column vectors

if ( max( delta ) < 1e-9 )
    out = 1;
else
    out = 0;
end