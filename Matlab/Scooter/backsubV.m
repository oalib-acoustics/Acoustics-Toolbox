function [x,b] = backsubV( N, mults, d, e, b )

% Performs back-substitution for a symmetric tridiagonal linear system
% [x,b] = backsub( N, mults, d, e, b )
%
% Based on TINVIT in EISPACK
% N is the order of the matrix
% d contains the diagonal of the upper triangular matrix
% e contains the upper diagonal (in its last N-1 positions)
% m contains the multipliers used during forward elimination
% b contains the right hand side (Ax=b)
% the answer is returned in vector b
%
% Michael B. Porter 7/1/85

% make sure input vectors are column vectors
%d = d(:);
%e = e(:);
%b = b(:);

M = size( b, 1 );

% Forward elimination
for I = 2 : N
    b( :, I ) = b( :, I ) - mults( :, I - 1 ) .* b( :, I - 1 );
end

x = zeros( M, N );  % preallocate and make sure result is a column vector

% Back-substitution (result in b)
x( :, N ) = b( :, N ) ./ d( :, N );
if ( N >= 2 )
    for I = N - 1 : -1 : 1
        x( :, I ) = ( b( :, I ) - e( :, I + 1 ) .* x( :, I + 1 ) ) ./ d( :, I );
    end
end
