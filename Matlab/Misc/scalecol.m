function B = scalecol( A, mults )

% Scales each column of A by the respective multiplier in mults
% Equivalent to B = A * diag( mults )
% except the full matrix doesn't need to be constructed
%
% mults can be either a row or a column vector
%
% usage: B = scalecol( A, mults )

mults = mults( : );  % make sure mults is a row vector
B     = zeros( size( A ) );

if ndims( A ) == 1
    B( : ) = mults .* A( : );
end

if ndims( A ) == 2
    Ncols = length( mults );
    B = full( A * spdiags( mults, 0, Ncols, Ncols ) );
    %B = A * diag( mults' );    % generates memory error, apparently trying to allocate diag( mults' );
    %for icol = 1: Ncols
    %    B( :, icol ) = mults( icol ) * A( :, icol );
    %end
end
