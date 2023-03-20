function [ mults, dt, et ] = factortriV( N, d, e )

% [ mults, dt, et ] = factor( N, d, e )
% Gaussian elimination to factor a symmetric tridiagonal linear system
% 
% N is the order of the matrix
% d contains the diagonal elements of the input matrix
% e              subdiagonal in its last N-1 positions
% dt contains the diagonal after reduction
% et contains the upper diagonal
% mults contains the multipliers used during elimination
%
% d and e can be matrices in which case the factorization is done for each
% row in parallel

% LU decomposition without interchanges

M = size( d, 1 );
dt    = zeros( M, N );
mults = zeros( M, N );

dt( 1 ) = d( 1 );

if N >= 2
    for I = 1 : N-1
        mults( :, I ) = e( :, I+1 ) ./ dt( :, I ); % multiplier
        dt( :, I+1  ) = d( :, I+1 ) - mults( :, I ) .* e( :, I+1 );  % new diagonal
    end
end
et = e;
